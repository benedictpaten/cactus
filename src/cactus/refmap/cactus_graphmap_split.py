#!/usr/bin/env python3

"""This sciprt will take the cactus-graphmap output and input and split it up into chromosomes.  This is done after a whole-genome graphmap as we need the whole-genome minigraph alignments to do the chromosome splitting in the first place.  For each chromosome, it will make a PAF, GFA and Seqfile (pointing to chrosome fastas)
"""
import os, sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit
import shutil
import re

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import unzip_gz, write_s3
from cactus.shared.common import get_faidx_subpath_rename_cmd
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.version import cactus_commit
from cactus.preprocessor.fileMasking import get_mask_bed_from_fasta
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from cactus.refmap.cactus_graphmap import filter_paf
from cactus.refmap.cactus_minigraph import check_sample_names, minigraph_gfa_from_pansn
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory, getTempFile, catFiles

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("seqFile", help = "Seq file (gzipped fastas supported)")
    parser.add_argument("minigraphGFA", help = "Minigraph-compatible reference graph in GFA format (can be gzipped)")
    parser.add_argument("graphmapPAF", type=str, help = "Output pairwise alignment file in PAF format (can be gzipped)")
    parser.add_argument("--outDir", required=True, type=str, help = "Output directory")
    parser.add_argument("--refContigs", nargs="*", help = "Subset to these reference contigs (multiple allowed)", default=[])
    parser.add_argument("--refContigsFile", type=str, help = "Subset to (newline-separated) reference contigs in this file")
    parser.add_argument("--otherContig", type=str, help = "Lump all reference contigs unselected by above options into single one with this name")
    parser.add_argument("--reference", required=True, nargs='+', type=str, help = "Name of reference (in seqFile).  Ambiguity filters will not be applied to it")
    parser.add_argument("--maskFilter", type=int, help = "Ignore softmasked sequence intervals > Nbp")
    parser.add_argument("--minIdentity", type=float, help = "Ignore PAF lines with identity (column 10/11) < this (overrides minIdentity in <graphmap_split> in config)")
    parser.add_argument("--permissiveContigFilter", nargs='?', const='0.25', default=None, type=float, help = "If specified, override the configuration to accept contigs so long as they have at least given fraction of coverage (0.25 if no fraction specified). This can increase sensitivity of very small, fragmented and/or diverse assemblies.")

    parser.add_argument("--collapse", help = "Incorporate minimap2 self-alignments.", action='store_true', default=False)
    
    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    if options.outDir and not options.outDir.startswith('s3://'):
        if not os.path.isdir(options.outDir):
            os.makedirs(options.outDir)
        
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    # support but ignore multi reference
    if options.reference:
        options.reference = options.reference[0]

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()
    cactus_graphmap_split(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap-split has finished after {} seconds".format(run_time))

def cactus_graphmap_split(options):
    with Toil(options) as toil:
        importSingularityImage(options)

        #load cactus config
        config_node = ET.parse(options.configFile).getroot()
        config = ConfigWrapper(config_node)
        config.substituteAllPredefinedConstantsWithLiterals(options)

        if options.collapse:
            findRequiredNode(config_node, "graphmap").attrib["collapse"] = 'all'

        #Run the workflow
        if options.restart:
            wf_output = toil.restart()
        else:
            options.cactusDir = getTempDirectory()
            
            # load up the contigs if any
            ref_contigs = set(options.refContigs)
            # todo: use import?
            if options.refContigsFile:
                with open(options.refContigsFile, 'r') as rc_file:
                    for line in rc_file:
                        if len(line.strip()):
                            ref_contigs.add(line.strip().split()[0])
            options.refContigs = list(ref_contigs)

            if options.refContigs:
                max_ref_contigs = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "maxRefContigs", typeFn=int, default=128)
                if len(options.refContigs) > max_ref_contigs:
                    logger.warning('You specified {} refContigs, which is greater than the suggested limit of {}. This may cause issues downstream.'.format(len(options.refContigs), max_ref_contigs))
                

            if options.otherContig:
                assert options.otherContig not in options.refContigs

            # get the minigraph "virutal" assembly name
            graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")

            # load the seqfile
            seqFile = SeqFile(options.seqFile, defaultBranchLen=config.getDefaultBranchLen(pangenome=True))

            #import the graph
            logger.info("Importing {}".format(options.minigraphGFA))
            gfa_id = toil.importFile(makeURL(options.minigraphGFA))

            #import the paf
            logger.info("Importing {}".format(options.graphmapPAF))
            paf_id = toil.importFile(makeURL(options.graphmapPAF))

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            input_seq_id_map = {}
            input_name_map = {}
            leaves = set([seqFile.tree.getName(node) for node in seqFile.tree.getLeaves()])
            
            if graph_event not in leaves:
                raise RuntimeError("Minigraph name {} not found in seqfile".format(graph_event))
            if options.reference and options.reference not in leaves:
                raise RuntimeError("Name given with --reference {} not found in seqfile".format(options.reference))

            # validate the sample names
            check_sample_names(seqFile.pathMap.keys(), options.reference)
                        
            for genome, seq in seqFile.pathMap.items():
                if genome in leaves:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    logger.info("Importing {}".format(seq))
                    input_seq_id_map[genome] = toil.importFile(seq)
                    input_name_map[genome] = os.path.basename(seq)

            # run the workflow
            wf_output = toil.start(Job.wrapJobFn(graphmap_split_workflow, options, config, input_seq_id_map, input_name_map,
                                                 gfa_id, options.minigraphGFA,
                                                 paf_id, options.graphmapPAF))

        #export the split data
        export_split_data(toil, wf_output[0], wf_output[1], wf_output[2], wf_output[3], options.outDir, config)

def graphmap_split_workflow(job, options, config, seq_id_map, seq_name_map, gfa_id, gfa_path, paf_id, paf_path, sanitize=True,
                            pansn_gfa_input=True):

    root_job = Job()
    job.addChild(root_job)

    # can be a list coming in from cactus-pangenome, but we only need first item
    if type(options.reference) is list:
        options.reference = options.reference[0]

    #override the minIdentity
    if options.minIdentity is not None:
        findRequiredNode(config.xmlRoot, "graphmap").attrib["minIdentity"] = str(options.minIdentity)

    #override the contig filter parameters with a single value taken from this option
    #also, disabling the uniqueness threshold
    if options.permissiveContigFilter is not None:
        findRequiredNode(config.xmlRoot, "graphmap_split").attrib["minQueryCoverages"] = str(options.permissiveContigFilter)
        findRequiredNode(config.xmlRoot, "graphmap_split").attrib["minQueryCoverageThresholds"] = ""
        findRequiredNode(config.xmlRoot, "graphmap_split").attrib["minQueryUniqueness"] = "1"
    
    # get the sizes before we overwrite below
    gfa_size = gfa_id.size
    paf_size = paf_id.size
    genome_names = set(seq_id_map.keys())    

    # fix up the headers
    if sanitize:
        sanitize_job = root_job.addChildJobFn(sanitize_fasta_headers, seq_id_map, pangenome=True)
        seq_id_map = sanitize_job.rv()
    else:
        sanitize_job = Job()
        root_job.addChild(sanitize_job)

    # auto-set --refContigs
    if not options.refContigs:
        refcontig_job = sanitize_job.addFollowOnJobFn(detect_ref_contigs, config, options, seq_id_map)
        ref_contigs = refcontig_job.rv()
        options.otherContig = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "otherContigName", typeFn=str, default="chrOther")
        sanitize_job = refcontig_job
    else:
        # move from options to avoid promise confusion
        ref_contigs = options.refContigs

    # convert the GFA from PanSN to Cactus names
    if pansn_gfa_input:
        rename_gfa_job = root_job.addChildJobFn(minigraph_gfa_from_pansn, genome_names, gfa_path, gfa_id,
                                                disk=gfa_size*2)
        new_root_job = Job()
        root_job.addFollowOn(new_root_job)
        root_job = new_root_job
        gfa_id = rename_gfa_job.rv()
            
    # use file extension to sniff out compressed input
    if gfa_path.endswith(".gz"):
        gfa_id = root_job.addChildJobFn(unzip_gz, gfa_path, gfa_id, delete_original=False, disk=gfa_size * 10).rv()
        gfa_size *= 10
    if paf_path.endswith(".gz"):
        paf_id = root_job.addChildJobFn(unzip_gz, paf_path, paf_id, delete_original=False, disk=paf_id.size * 10).rv()
        paf_size *= 10

    # do some basic paf filtering
    paf_filter_mem = max(paf_id.size * 10, 2**32)
    paf_filter_job = root_job.addFollowOnJobFn(filter_paf, paf_id, config, reference=options.reference,
                                               skip_overlap_filter=True,
                                               disk = paf_id.size * 10, memory=cactus_clamp_memory(paf_filter_mem))
    paf_id = paf_filter_job.rv()
    root_job = paf_filter_job

    mask_bed_id = None
    if options.maskFilter:
        mask_bed_id = sanitize_job.addFollowOnJobFn(get_mask_bed, seq_id_map, options.maskFilter).rv()
        
    # use rgfa-split to split the gfa and paf up by contig
    split_gfa_job = root_job.addFollowOnJobFn(split_gfa, config, gfa_id, [paf_id], ref_contigs,
                                              options.otherContig, options.reference, mask_bed_id,
                                              disk=(gfa_size + paf_size) * 5,
                                              memory=cactus_clamp_memory((gfa_size + paf_size) * 3))

    # use the output of the above splitting to do the fasta splitting
    split_fas_job = split_gfa_job.addFollowOnJobFn(split_fas, seq_id_map, seq_name_map, split_gfa_job.rv(0))

    # gather everythign up into a table
    gather_fas_job = split_fas_job.addFollowOnJobFn(gather_fas, split_gfa_job.rv(0), split_fas_job.rv(0), split_fas_job.rv(1))

    # lump "other" contigs together into one file (to make fewer align jobs downstream)
    bin_other_job = gather_fas_job.addFollowOnJobFn(bin_other_contigs, config, ref_contigs, options.otherContig, gather_fas_job.rv(0),
                                                    disk=(gfa_size + paf_size) * 2)

    # return all the files, as well as the 2 split logs
    return (seq_name_map, bin_other_job.rv(), split_gfa_job.rv(1), gather_fas_job.rv(1))

def detect_ref_contigs(job, config, options, seq_id_map):
    """ automatically determine --refContigs from the data """
    work_dir = job.fileStore.getLocalTempDir()
    # it will be unzipped since we're after sanitize
    fa_path = os.path.join(work_dir, options.reference + '.fa')
    job.fileStore.readGlobalFile(seq_id_map[options.reference], fa_path)
    cactus_call(parameters=['samtools', 'faidx', fa_path])
    contigs = []
    with open(fa_path + '.fai', 'r') as fai_file:
        for line in fai_file:
            toks = line.split('\t')
            assert len(toks) > 2
            # clean out prefix that was added by sanitize
            assert toks[0].startswith('id={}|'.format(options.reference))
            clean_name = toks[0][len('id={}|'.format(options.reference)):]
            contigs.append((clean_name, float(toks[1])))

    # get the cutoff heuristics 
    max_ref_contigs = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "maxRefContigs", typeFn=int, default=128)
    ref_contig_dropoff = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "refContigDropoff", typeFn=float, default=10.0)
    ref_contig_regex = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "refContigRegex", typeFn=str, default=None)
    other_contig_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "otherContigName", typeFn=str, default="chrOther")

    # sort by decreasing size
    sorted_contigs = sorted(contigs, key=lambda x : x[1], reverse=True)
    
    ref_contigs = set()
    # pass 1: do the regex matches
    if ref_contig_regex:
        for contig, length in contigs:
            if re.fullmatch(ref_contig_regex, contig):
                ref_contigs.add(contig)
                if len(ref_contigs) >= max_ref_contigs:
                    break
        
    # pass 2: do the dropoff
    longest_length = sorted_contigs[0][1]
    for contig_len in sorted_contigs:
        if len(ref_contigs) >= max_ref_contigs:
            break
        if contig_len[0] not in ref_contigs:
            dropoff = longest_length / contig_len[1]
            if dropoff >= ref_contig_dropoff:
                break
            ref_contigs.add(contig_len[0])

    ref_contigs = sorted(list(ref_contigs))

    msg = "auto-detected --refContigs {}".format(' '.join(ref_contigs))
    if len(ref_contigs) < len(sorted_contigs):
        msg += " --otherContig {}".format(other_contig_name)

    RealtimeLogger.info(msg)

    return sorted(list(ref_contigs))
    

def get_mask_bed(job, seq_id_map, min_length):
    """ make a bed file from the fastas """
    beds = []
    for event in seq_id_map.keys():
        fa_id = seq_id_map[event]
        fa_path = '{}.fa'.format(event) 
        beds.append(job.addChildJobFn(get_mask_bed_from_fasta, event, fa_id, fa_path, min_length, disk=fa_id.size * 5).rv())
    return job.addFollowOnJobFn(cat_beds, beds).rv()

def cat_beds(job, bed_ids):
    in_beds = [job.fileStore.readGlobalFile(bed_id) for bed_id in bed_ids]
    out_bed = job.fileStore.getLocalTempFile()
    catFiles(in_beds, out_bed)
    return job.fileStore.writeGlobalFile(out_bed)
    
def split_gfa(job, config, gfa_id, paf_ids, ref_contigs, other_contig, reference_event, mask_bed_id):
    """ Use rgfa-split to divide a GFA and PAF into chromosomes.  The GFA must be in minigraph RGFA output using
    the desired reference. """

    if not paf_ids:
        # we can bypass when, ex, doing second pass on ambiguous sequences but not are present
        return [None, None]

    if not gfa_id and not getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "remap", typeFn=bool, default=False):
        # also bypass if remapping is off in the config (we know it's the second pass because gfa_id is None)
        return [None, None]
    
    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    paf_path = os.path.join(work_dir, "mg.paf")
    out_prefix = os.path.join(work_dir, "split_")
    bed_path = os.path.join(work_dir, "mask.bed")
    log_path = os.path.join(work_dir, "split.log")
    if (mask_bed_id):
        job.fileStore.readGlobalFile(mask_bed_id, bed_path)

    if gfa_id:
        job.fileStore.readGlobalFile(gfa_id, gfa_path)
        
    paf_paths = []
    for i, paf_id in enumerate(paf_ids):
        paf_paths.append('{}.{}'.format(paf_path, i) if len(paf_ids) > 1 else paf_path)
        job.fileStore.readGlobalFile(paf_id, paf_paths[-1])
    if len(paf_paths) > 1:
        catFiles(paf_paths, paf_path)

    # get the minigraph "virutal" assembly name
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    # and look up its unique id prefix.  this will be needed to pick its contigs out of the list
    mg_id = graph_event

    # get the specificity filters
    query_coverages = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQueryCoverages", default=None)
    query_coverage_thresholds = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQueryCoverageThresholds", default=None)
    coverage_opts = []
    if query_coverages is not None:
        try:
            query_cov_vals = [float(x) for x in query_coverages.split()]
            query_thresh_vals = [int(x) for x in query_coverage_thresholds.split()]
            assert len(query_cov_vals) == len(query_thresh_vals) + 1
            assert sorted(query_thresh_vals) == query_thresh_vals
            for cv in query_cov_vals:
                coverage_opts += ['-n', str(cv)]
            for tv in query_thresh_vals:
                coverage_opts += ['-T', str(tv)]
        except:
            raise RuntimeError("minQueryCoverages and / or minQueryCoverageThresholds malspecified in config")
    query_uniqueness = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQueryUniqueness", default="0")
    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")

    cmd = ['rgfa-split',
           '-p', paf_path,
           '-b', out_prefix,
           '-Q', query_uniqueness,
           '-a', amb_name,
           '-L', log_path]
    cmd += coverage_opts    
    if gfa_id:
        cmd += ['-g', gfa_path, '-G']
    if reference_event:
        cmd += ['-r', 'id={}|'.format(reference_event)]
    if mask_bed_id:
        cmd += ['-B', bed_path]
    min_mapq = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minMAPQ")
    if min_mapq:
        cmd += ['-A', min_mapq]
    # optional stuff added to second pass:
    if not gfa_id:
        remap_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "remapSplitOptions", default=None)
        if remap_opts:
            cmd += remap_opts.split(' ')
    if not other_contig:
        # we use this option only when subsetting to given ref contigs *without* using the "other"
        # option.  otherwise, we rely on bin_other_contigs here in python to do the job
        for contig in ref_contigs:
            cmd += ['-c', contig]
            # hack in to make sure we deal with gfas with/without prexixes
            # todo: once all merged, this needs to be simplified and cleaned
            if reference_event and not contig.startswith('id={}|'.format(reference_event)):
                cmd += ['-c', 'id={}|{}'.format(reference_event, contig)]
            if contig.startswith('id=') and contig.find('|') > 3:
                cmd += ['-c', contig[contig.find('|')+1:]]            

    cactus_call(parameters=cmd, work_dir=work_dir, job_memory=job.memory)

    output_id_map = {}
    for out_name in os.listdir(work_dir):
        file_name, ext = os.path.splitext(out_name)
        if file_name.startswith(os.path.basename(out_prefix)) and ext in [".gfa", ".paf", ".fa_contigs"] and \
           os.path.isfile(os.path.join(work_dir, file_name + ".fa_contigs")):
            name = file_name[len(os.path.basename(out_prefix)):]

            # don't leave unique identifier in names
            if name.startswith('id=') and name.find('|') > 3:
                name=name[name.find('|') + 1:]

            if name not in output_id_map:
                output_id_map[name] = {}
            if ext == '.paf':
                # apply the hacky naming correction so that subpaths have no special characterse in the hal (to make hubs happy)
                # this gets undone by hal2vg
                cmd = get_faidx_subpath_rename_cmd()
                cmd += ['-e', 's/ /\t/g', '-i', os.path.join(work_dir, out_name)]
                cactus_call(parameters=cmd)
            output_id_map[name][ext[1:]] = job.fileStore.writeGlobalFile(os.path.join(work_dir, out_name))
            
    return output_id_map, job.fileStore.writeGlobalFile(log_path)

def split_fas(job, seq_id_map, seq_name_map, split_id_map):
    """ Use samtools to split a bunch of fasta files into reference contigs, using the output of rgfa-split as a guide"""

    if (seq_id_map, split_id_map) == (None, None):
        return None

    root_job = Job()
    job.addChild(root_job)

    # map event name to dict of contgs.  ex fa_contigs["CHM13"]["chr13"] = file_id
    fa_contigs = {}
    # map event name to dict of contig sizes ex fa_contigs["CHM13"]["chr13"] = N
    # where N is total number of fasta bases
    fa_contig_sizes = {}
    
    # we do each fasta in parallel
    for event, fa_id in seq_id_map.items():
        fa_path = seq_name_map[event]
        if fa_id.size:
            split_job = root_job.addChildJobFn(split_fa_into_contigs, event, fa_id, fa_path, split_id_map,
                                               strip_prefix=False, 
                                               disk=fa_id.size * 3)
            fa_contigs[event] = split_job.rv(0)
            fa_contig_sizes[event] = split_job.rv(1)

    return fa_contigs, fa_contig_sizes

def split_fa_into_contigs(job, event, fa_id, fa_name, split_id_map, strip_prefix=False):
    """ Use samtools turn on fasta into one for each contig. this relies on the informatino in .fa_contigs
    files made by rgfa-split """

    # download the fasta
    work_dir = job.fileStore.getLocalTempDir()
    fa_path = os.path.join(work_dir, '{}.fa'.format(event))
    is_gz = fa_name.endswith(".gz")
    job.fileStore.readGlobalFile(fa_id, fa_path)

    unique_id = 'id={}|'.format(event)
                        
    contig_fa_dict = {}
    contig_size_dict = {}

    for ref_contig in split_id_map.keys():
        query_contig_list_id = split_id_map[ref_contig]['fa_contigs']
        list_path = os.path.join(work_dir, '{}.fa_contigs'.format(ref_contig))
        job.fileStore.readGlobalFile(query_contig_list_id, list_path)
        faidx_input_path = os.path.join(work_dir, '{}.fa_contigs.clean'.format(ref_contig))
        contig_count = 0
        with open(list_path, 'r') as list_file, open(faidx_input_path, 'w') as clean_file:
            for line in list_file:
                query_contig = line.strip()
                if query_contig.startswith(unique_id):
                    if strip_prefix:
                        query_contig = query_contig[len(unique_id):]                        
                    clean_file.write('{}\n'.format(query_contig))
                    contig_count += 1
        contig_fasta_path = os.path.join(work_dir, '{}_{}.fa'.format(event, ref_contig))
        if contig_count > 0:
            cmd = [['samtools', 'faidx', fa_path, '--region-file', faidx_input_path]]
            # transform chr1:10-15 (1-based inclusive) into chr1_sub_9_15 (0-based end open)
            # this is a format that contains no special characters in order to make assembly hubs
            # happy.  But it does require conversion going into vg which wants chr[9-15] and
            # hal2vg is updated to do this autmatically
            cmd.append(get_faidx_subpath_rename_cmd())
            if is_gz:
                cmd.append(['bgzip', '--threads', str(job.cores)])
            cactus_call(parameters=cmd, outfile=contig_fasta_path)

            # now get the total sequence size
            size_cmd = [['cat', contig_fasta_path], ['grep', '-v', '>'], ['wc'], ['awk', '{printf \"%f\", $3-$1}']]
            if is_gz:
                size_cmd[0] = ['bgzip', '-dc', contig_fasta_path, '--threads', str(job.cores)]
            num_bases = int(float(cactus_call(parameters=size_cmd, check_output=True).strip()))
        else:
            # TODO: review how cases like this are handled
            with open(contig_fasta_path, 'w') as empty_file:
                empty_file.write("")
            num_bases = 0
        contig_fa_dict[ref_contig] = job.fileStore.writeGlobalFile(contig_fasta_path)
        contig_size_dict[ref_contig] = num_bases
        
    return contig_fa_dict, contig_size_dict

def gather_fas(job, output_id_map, contig_fa_map, contig_size_map):
    """ take the split_fas output which has everything sorted by event, and move into the ref-contig-based table
    from split_gfa.  return the updated table, which can then be exported into the chromosome projects """

    if not contig_fa_map:
        return None

    events = set()
    for ref_contig in output_id_map.keys():
        output_id_map[ref_contig]['fa'] = {}
        for event, fa_id in contig_fa_map.items():
            if ref_contig in fa_id:
                output_id_map[ref_contig]['fa'][event] = fa_id[ref_contig]
            events.add(event)

    # go ahead and make the size table file here
    # download the fasta
    work_dir = job.fileStore.getLocalTempDir()
    size_table_path = os.path.join(work_dir, 'contig_sizes.tsv')
    events = sorted(events)
    with open(size_table_path, 'w') as size_table_file:
        #write the header
        size_table_file.write("Contig\t" + '\t'.join(events) + '\tmin\tmax\tavg\n')
        for ref_contig in output_id_map.keys():
            sizes = [contig_size_map[event][ref_contig] for event in events]
            #total size for each event
            size_table_file.write(ref_contig + '\t' + '\t'.join([str(s) for s in sizes]))
            #tack on 3 basic stats
            size_table_file.write('\t{}\t{}\t{}\n'.format(min(sizes), max(sizes), int(sum(sizes) / len(sizes)) if sizes else 0))
    contig_size_table_id = job.fileStore.writeGlobalFile(size_table_path)

    return output_id_map, contig_size_table_id

def bin_other_contigs(job, config, ref_contigs, other_contig, output_id_map):
    """ take all the other (ie non ref) contigs (in practice, unplaced bits of grch38) and merge them
    all up into a single "other" contig.  this avoids a 1000 align jobs getting created downstream. 
    Note we've moved this to the ned here (it used to be done usinc -c -o in rgfa-split) so as to
    avoid letting these contigs glom together into componenets: the final output will be identital
    to if they were kept in separate files. """
    work_dir = job.fileStore.getLocalTempDir()
    if not ref_contigs or not other_contig or not output_id_map:
        return output_id_map

    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")
    
    ref_contigs = set(ref_contigs)
    other_contigs = set()
    merged_output_path_map = {} # map[TYPE][EVENT] = local file path

    for ref_contig in output_id_map.keys():
        if ref_contig not in ref_contigs and ref_contig != amb_name:
            other_contigs.add(ref_contig)
            for type_key in output_id_map[ref_contig].keys():
                if type_key == 'fa':
                    for event in output_id_map[ref_contig][type_key].keys():
                        # we read the existing id
                        local_path = os.path.join(work_dir, '{}_{}.{}'.format(ref_contig, event, type_key))
                        job.fileStore.readGlobalFile(output_id_map[ref_contig][type_key][event], local_path)

                        # find its new home
                        other_path = os.path.join(work_dir, '{}_{}.{}'.format(other_contig, event, type_key))

                        # update the map
                        if type_key not in merged_output_path_map:
                            merged_output_path_map[type_key] = {}
                        if event not in merged_output_path_map[type_key]:
                            merged_output_path_map[type_key][event] = other_path

                        # append the file
                        with open(local_path, 'rb') as local_file, open(other_path, 'ab') as other_file:
                            shutil.copyfileobj(local_file, other_file)
                else:
                    # we read the existing id
                    local_path = os.path.join(work_dir, '{}.{}'.format(ref_contig, type_key))
                    job.fileStore.readGlobalFile(output_id_map[ref_contig][type_key], local_path)

                    # find its new home
                    other_path = os.path.join(work_dir, '{}.{}'.format(other_contig, type_key))

                    # update the map
                    if type_key not in merged_output_path_map:
                        merged_output_path_map[type_key] = other_path

                    # append the file
                    with open(local_path, 'rb') as local_file, open(other_path, 'ab') as other_file:
                        shutil.copyfileobj(local_file, other_file)
                    
    # and convert to ids
    for type_key in merged_output_path_map.keys():
        if type_key == 'fa':
            for event in merged_output_path_map[type_key].keys():
                merged_output_path_map[type_key][event] = job.fileStore.writeGlobalFile(merged_output_path_map[type_key][event])
        else:
            merged_output_path_map[type_key] = job.fileStore.writeGlobalFile(merged_output_path_map[type_key])

    # remove the other contigs from the map
    for contig in other_contigs:
        del output_id_map[contig]

    # drop in the other contig
    output_id_map[other_contig] = merged_output_path_map

    return output_id_map

def export_split_data(toil, input_name_map, output_id_map, split_log_id, contig_size_table_id, output_dir, config):
    """ download all the split data locally """

    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")
    
    chrom_file_map = {}

    # export the sizes
    contig_size_table_path = os.path.join(output_dir, 'contig_sizes.tsv')
    toil.exportFile(contig_size_table_id, makeURL(contig_size_table_path))

    # export the log
    toil.exportFile(split_log_id, makeURL(os.path.join(output_dir, 'minigraph.split.log')))

    # hack to filter out reference contigs where nothing maps to it (not even itself)
    # this is usually worked around by using --refContigs <chroms> --otherContig chrOther
    # but filtering them out here lets users who don't do that make graphs
    # 
    empty_contigs = set()
    
    for ref_contig in output_id_map.keys():        
        if output_id_map[ref_contig] is None or len(output_id_map[ref_contig]) == 0:
            # todo: check ambigous?
            continue
        ref_contig_path = os.path.join(output_dir, ref_contig)
        if not os.path.isdir(ref_contig_path) and not ref_contig_path.startswith('s3://'):
            os.makedirs(ref_contig_path)

        # GFA: <output_dir>/<contig>/<contig>.gfa
        if 'gfa' in output_id_map[ref_contig]:
            # we do this check because no gfa made for ambiguous sequences "contig"
            toil.exportFile(output_id_map[ref_contig]['gfa'], makeURL(os.path.join(ref_contig_path, '{}.gfa'.format(ref_contig))))

        # PAF: <output_dir>/<contig>/<contig>.paf
        paf_path = os.path.join(ref_contig_path, '{}.paf'.format(ref_contig))
        if 'paf' in output_id_map[ref_contig]:
            # rgfa-split doesn't write empty pafs: todo: should it?
            toil.exportFile(output_id_map[ref_contig]['paf'], makeURL(paf_path))

        # Fasta: <output_dir>/<contig>/fasta/<event>_<contig>.fa ..
        seq_file_map = {}
        for event, ref_contig_fa_id in output_id_map[ref_contig]['fa'].items():
            fa_base = os.path.join(ref_contig_path, 'fasta')
            if not os.path.isdir(fa_base) and not fa_base.startswith('s3://'):
                os.makedirs(fa_base)
            fa_path = makeURL(os.path.join(fa_base, '{}_{}.fa'.format(event, ref_contig)))
            if input_name_map[event].endswith('.gz'):
                fa_path += '.gz'
            seq_file_map[event] = fa_path
            toil.exportFile(ref_contig_fa_id, fa_path)

        seq_file_map_size = 0
        for event, fa_path in seq_file_map.items():
            if output_id_map[ref_contig]['fa'][event].size > 0:
                seq_file_map_size += 1
        if seq_file_map_size < 2:
            logger.warning("Omitting reference contig {} from the graph because it doesn't align well enough. If you absolutely want to include it, rerun with --refContigs and --otherContig specified".format(ref_contig))
            empty_contigs.add(ref_contig)
            continue
        
        # Seqfile: <output_dir>/seqfiles/<contig>.seqfile
        seq_file_path = os.path.join(output_dir, 'seqfiles', '{}.seqfile'.format(ref_contig))
        if seq_file_path.startswith('s3://'):
            seq_file_temp_path = getTempFile()
        else:
            seq_file_temp_path = seq_file_path
            if not os.path.isdir(os.path.dirname(seq_file_path)):
                os.makedirs(os.path.dirname(seq_file_path))
        with open(seq_file_temp_path, 'w') as seq_file:
            for event, fa_path in seq_file_map.items():
                # cactus can't handle empty fastas.  if there are no sequences for a sample for this
                # contig, just don't add it.
                if output_id_map[ref_contig]['fa'][event].size > 0:
                    seq_file.write('{}\t{}\n'.format(event, fa_path))
        if seq_file_path.startswith('s3://'):
            write_s3(seq_file_temp_path, seq_file_path)

        # Top-level seqfile
        chrom_file_map[ref_contig] = seq_file_path, paf_path

    # Chromfile : <coutput_dir>/chromfile.txt
    chrom_file_path = os.path.join(output_dir, 'chromfile.txt')
    if chrom_file_path.startswith('s3://'):
        chrom_file_temp_path = getTempFile()
    else:
        chrom_file_temp_path = chrom_file_path        
    with open(chrom_file_temp_path, 'w') as chromfile:
        for ref_contig, seqfile_paf in chrom_file_map.items():
            if ref_contig != amb_name and ref_contig not in empty_contigs:
                seqfile, paf = seqfile_paf[0], seqfile_paf[1]
                chromfile.write('{}\t{}\t{}\n'.format(ref_contig, seqfile, paf))
    if chrom_file_path.startswith('s3://'):
        write_s3(chrom_file_temp_path, chrom_file_path)
        
if __name__ == "__main__":
    main()
