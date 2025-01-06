#!/usr/bin/env python3

"""This script will take the output of cactus-align-batch (which was run on the output of cactus-graphmap-split).  This would be
a set of .vg files, one for each chromosome.  It will do the following in order to make a single output graph
- clip unmapped regions out of each graph
- join the ids of the clipped graph
- convert each clipped graph to gfa
- merge the gfas into a whole-genome graph
- convert the gfa into a gbwt/gbwt-graph pair
- export the gbwt to xg
- export the xg/gbwt to vcf

so the final output will be
- set of clipped vg graphs
- graph.gfa.gz
- graph.gbwt
- graph.gg
- graph.xg
- graph.snarls
- graph.vcf.gz
- graph.vcf.gz.tbi
"""
import os, sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit
import re
import subprocess

from operator import itemgetter
from collections import defaultdict

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
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import clean_jobstore_files
from cactus.shared.version import cactus_commit
from cactus.preprocessor.fileMasking import get_mask_bed_from_fasta
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from cactus.refmap.cactus_minigraph import check_sample_names
from cactus.progressive.cactus_prepare import human2bytesN

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory, getTempFile, catFiles
    
def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("--vg", required=True, nargs='+',  help = "Input vg files (PackedGraph or HashGraph format)")
    parser.add_argument("--hal", nargs='+', default = [], help = "Input hal files (for merging)")    
    parser.add_argument("--outDir", required=True, type=str, help = "Output directory")
    parser.add_argument("--outName", required=True, type=str, help = "Basename of all output files")
    parser.add_argument("--reference", required=True, nargs='+', type=str, help = "Reference event name(s). The first will be the \"true\" reference and will be left unclipped and uncollapsed. It also should have been used with --reference in all upstream commands. Other names will be promoted to reference paths in vg")

    graphmap_join_options(parser)
        
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

    #make sure our options make sense and fill in sensible defaults
    options = graphmap_join_validate_options(options)
        
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()
    graphmap_join(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap-join has finished after {} seconds".format(run_time))

def graphmap_join_options(parser):
    """ we share these options with cactus-pangenome """
    parser.add_argument("--clip", type=int, default=10000, help = "Generate clipped graph by removing anything longer than this amount that is unaligned to the underlying minigraph. Set to 0 to disable (must also set --filter 0 as well). [default=10000]")
    
    parser.add_argument("--filter", type=int, default=2, help = "Generate a frequency filtered graph (from the clipped graph) by removing any sequence present in fewer than this many sequences. Set to 0 to disable. [default=2]")

    parser.add_argument("--gfa", nargs='*', default=None, help = "Produce a GFA for given graph type(s) if specified. Valid types are 'full', 'clip', and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space. [--gfa clip assumed by default]")

    parser.add_argument("--unchopped-gfa", nargs='*', default=None, help = "Many vg tools require nodes have lengths of at most 1024bp, so by default all output graphs are chopped in this way. This option will produce GFA files where this limit is not imposed - their nodes are as long as the graph toplogy allows. These files will have the .unchopped.gfa.gz suffix and *not* be compatible with any other graph files or VCFs which are all based on the chopped vg coordinates. Valid types are 'full', 'clip', and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space. [--gfa clip assumed by default]")
    
    parser.add_argument("--gbz", nargs='*', default=None, help = "Generate GBZ/snarls indexes for the given graph type(s) if specified. Valid types are 'full', 'clip' and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space. --giraffe will also produce these (and other) indexes")

    parser.add_argument("--xg", nargs='*', default=None, help = "Generate XG index (from GBZ) for the given graph type(s) if specified. Valid types are 'full', 'clip' and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space.")    

    parser.add_argument("--odgi", nargs='*', default=None, help = "Generate ODGI (.og) versions for the given graph type(s) if specified. Valid types are 'full', 'clip' and 'filter'. If no type specified 'full' will be used. Multiple types can be provided separated by a space.")

    parser.add_argument("--viz", nargs='*', default=None, help = "Generate 1D ODGI visualizations of each chromosomal graph for the given graph type(s) if specified. Valid types are 'full', 'clip' (but NOT `filter`). If no type specified 'full' will be used. Multiple types can be provided separated by a space.")

    parser.add_argument("--draw", nargs='*', default=None, help = "Generate 2D ODGI visualizations of each chromosomal graph for the given graph type(s) if specified. WARNING: EXPERIMENTAL OPTION. More work needs to be done to figure out the best ODGI parameters to use, and it can be quite slow in the meantime. For larger graphs, use --chrom-og and run the drawing by hand in order to avoid having draw issues prevent you from getting the rest of your output. Valid types are 'full', 'clip' (but NOT `filter`). If no type specified 'full' will be used. Multiple types can be provided separated by a space.")
    
    parser.add_argument("--chrom-vg", nargs='*', default=None, help = "Produce a directory of chromosomal graphs is vg format for the graph type(s) specified. Valid typs are 'full', 'clip' and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space.  The output will be the <outDir>/<outName>.chroms/ directory")

    parser.add_argument("--chrom-og", nargs='*', default=None, help = "Produce a directory of chromosomal graphs is odgi format for the graph type(s) specified. Valid typs are 'full', 'clip' and 'filter'. If no type specified 'full' will be used. Multiple types can be provided separated by a space.  The output will be the <outDir>/<outName>.chroms/ directory")
    
    parser.add_argument("--vcf", nargs='*', default=None, help = "Generate a VCF from the given graph type(s). Valid types are 'full', 'clip' and 'filter'. If no type specified, 'clip' will be used ('full' used if clipping disabled). Multipe types can be provided separated by space")
    parser.add_argument("--vcfReference", nargs='+', default=None, help = "If multiple references were provided with --reference, this option can be used to specify a subset for vcf creation with --vcf. By default, --vcf will create VCFs for the first reference only")
    parser.add_argument("--vcfbub", type=int, default=100000, help = "Use vcfbub to flatten nested sites (sites with reference alleles > this will be replaced by their children)). Setting to 0 will disable, only prudcing full VCF [default=100000].")
    parser.add_argument("--vcfwave", action='store_true', default=False, help = "Create a vcfwave-normalized VCF. vcfwave realigns alt alleles to the reference, and can help correct messy regions in the VCF. This option will output an additional VCF with 'wave' in its filename, other VCF outputs will not be affected")
    parser.add_argument("--vcfwaveCores", type=int, help = "Number of cores for each vcfwave job [default=2].", default=2)
    parser.add_argument("--vcfwaveMemory", type=human2bytesN, help = "Memory for reach vcfwave job [default=32Gi].", default=32000000000)
    parser.add_argument("--snarlStats", nargs='*', help = "Write a list of snarl statistics for the graph type(s). Valid types are 'full', 'clip' and 'filter'. If no type specified, 'clip' will be used ('full' used if clipping disabled). Multipe types can be provided separated by space")
        
    parser.add_argument("--giraffe", nargs='*', default=None, help = "Generate Giraffe (.dist, .min) indexes for the given graph type(s). Valid types are 'full', 'clip' and 'filter'. If not type specified, 'filter' will be used (will fall back to 'clip' than full if filtering, clipping disabled, respectively). Multiple types can be provided seperated by a space. NOTE: do not use this option if you want to use haplotype sampling. Use --haplo instead.")

    parser.add_argument("--haplo", nargs='*', default=None, help = "Generate haplotype subsampling (.ri, .hapl) indexes for the given graph type(s). Haplotype subsampling is a new, better alternative filtering by allele frequency. Valid types are 'full' and 'clip'. If not type specified, 'clip' will be used ('full' will be used if clipping disabled). Multiple types can be provided seperated by a space. NOTE: you normally want to use this option without --giraffe!")
    
    parser.add_argument("--indexCores", type=int, default=None, help = "cores for general indexing and VCF constructions (defaults to the same as --maxCores)")

    
    parser.add_argument("--indexMemory", type=human2bytesN,
                        help="Memory in bytes for each indexing and vcf construction job (defaults to an estimate based on the input data size). If specified will also be used to upper-bound per-chromosome memory estimates -- ie no job will request more than this much memory."
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))", default=None)   
    
    parser.add_argument("--collapse", help = "Incorporate minimap2 self-alignments.", action='store_true', default=False)

def graphmap_join_validate_options(options):
    """ make sure the options make sense and fill in sensible defaults """
    if options.hal and len(options.hal) != len(options.vg):
        raise RuntimeError("If --hal and --vg should specify the same number of files")

    # apply cpu override
    if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
        if not options.indexCores:
            options.indexCores = sys.maxsize
        options.indexCores = min(options.indexCores, max(1, cactus_cpu_count() - 1), int(options.maxCores) if options.maxCores else sys.maxsize)
    else:
        if not options.indexCores:
            raise RuntimeError("--indexCores required run *not* running on single machine batch system")
    
    # sanity check the workflow options and apply defaults
    if options.filter and not options.clip:
        raise RuntimeError('--filter cannot be used without also disabling --clip.')

    # check the reference name suffix
    check_sample_names(options.reference, options.reference[0])

    # if no graph output specified, default to gfa.gz
    if not options.gfa and not options.unchopped_gfa and not options.xg and not options.gbz:
        options.gfa = ['clip'] if options.clip else ['full']
    options.gfa = list(set(options.gfa)) if options.gfa else []
    for gfa in options.gfa:
        if gfa not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --gfa: {}. Must be one of {{clip, filter, full}}'.format(gfa))
        if gfa == 'clip' and not options.clip:
            raise RuntimError('--gfa cannot be set to clip since clipping is disabled')
        if gfa == 'filter' and not options.filter:
            raise RuntimeError('--gfa cannot be set to filter since filtering is disabled')

    if options.unchopped_gfa == []:
        options.unchopped_gfa = ['clip'] if options.clip else ['full']
    options.unchopped_gfa = list(set(options.unchopped_gfa)) if options.unchopped_gfa else []
    for gfa in options.unchopped_gfa:
        if gfa not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --unchopped-gfa: {}. Must be one of {{clip, filter, full}}'.format(gfa))
        if gfa == 'clip' and not options.clip:
            raise RuntimError('--unchopped-gfa cannot be set to clip since clipping is disabled')
        if gfa == 'filter' and not options.filter:
            raise RuntimeError('--unchopped-gfa cannot be set to filter since filtering is disabled')
        
    if options.gbz == []:
        options.gbz = ['clip'] if options.clip else ['full']
    options.gbz = list(set(options.gbz)) if options.gbz else []
    for gbz in options.gbz:
        if gbz not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --gbz: {}. Must be one of {{clip, filter, full}}'.format(gbz))
        if gbz == 'clip' and not options.clip:
            raise RuntimError('--gbz cannot be set to clip since clipping is disabled')
        if gbz == 'filter' and not options.filter:
            raise RuntimeError('--gbz cannot be set to filter since filtering is disabled')

    if options.xg == []:
        options.xg = ['clip'] if options.clip else ['full']
    options.xg = list(set(options.xg)) if options.xg else []
    for xg in options.xg:
        if xg not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --xg: {}. Must be one of {{clip, filter, full}}'.format(xg))
        if xg == 'clip' and not options.clip:
            raise RuntimError('--xg cannot be set to clip since clipping is disabled')
        if xg == 'filter' and not options.filter:
            raise RuntimeError('--xg cannot be set to filter since filtering is disabled')
        if xg not in options.xg:
            raise RuntimeError('Specifying --xg {} requires also specifying --gbz {}'.format(xg, xg))

    if options.odgi == []:
        options.odgi = ['full']
    options.odgi = list(set(options.odgi)) if options.odgi else []
    for odgi in options.odgi:
        if odgi not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --odgi: {}. Must be one of {{clip, filter, full}}'.format(odgi))
        if odgi == 'clip' and not options.clip:
            raise RuntimError('--odgi cannot be set to clip since clipping is disabled')
        if odgi == 'filter' and not options.filter:
            raise RuntimeError('--odgi cannot be set to filter since filtering is disabled')
        
    if options.viz == []:
        options.viz = ['full']
    options.viz = list(set(options.viz)) if options.viz else []
    for viz in options.viz:
        if viz not in ['clip', 'full']:
            raise RuntimeError('Unrecognized value for --viz: {}. Must be one of {{clip, full}}'.format(viz))
        if viz == 'clip' and not options.clip:
            raise RuntimError('--viz cannot be set to clip since clipping is disabled')
        if viz == 'filter' and not options.filter:
            raise RuntimeError('--viz cannot be set to filter since filtering is disabled')

    if options.draw == []:
        options.draw = ['full']
    options.draw = list(set(options.draw)) if options.draw else []
    for draw in options.draw:
        if draw not in ['clip', 'full']:
            raise RuntimeError('Unrecognized value for --draw: {}. Must be one of {{clip, full}}'.format(draw))
        if draw == 'clip' and not options.clip:
            raise RuntimError('--draw cannot be set to clip since clipping is disabled')
        if draw == 'filter' and not options.filter:
            raise RuntimeError('--draw cannot be set to filter since filtering is disabled')                

    if options.chrom_vg == []:
        options.chrom_vg = ['clip'] if options.clip else ['full']
    options.chrom_vg = list(set(options.chrom_vg)) if options.chrom_vg else []
    for chrom_vg in options.chrom_vg:
        if chrom_vg not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --chrom-vg: {}. Must be one of {{clip, filter, full}}'.format(chrom_vg))
        if chrom_vg == 'clip' and not options.clip:
            raise RuntimError('--chrom-vg cannot be set to clip since clipping is disabled')
        if chrom_vg == 'filter' and not options.filter:
            raise RuntimeError('--chrom-vg cannot be set to filter since filtering is disabled')

    if options.chrom_og == []:
        options.chrom_og = ['full']
    options.chrom_og = list(set(options.chrom_og)) if options.chrom_og else []
    for chrom_og in options.chrom_og:
        if chrom_og not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --chrom-og: {}. Must be one of {{clip, filter, full}}'.format(chrom_og))
        if chrom_og == 'clip' and not options.clip:
            raise RuntimError('--chrom-og cannot be set to clip since clipping is disabled')
        if chrom_og == 'filter' and not options.filter:
            raise RuntimeError('--chrom-og cannot be set to filter since filtering is disabled')                
    
    if options.vcf == []:
        options.vcf = ['clip'] if options.clip else ['full']
    options.vcf = list(set(options.vcf)) if options.vcf else []    
    for vcf in options.vcf:
        if vcf not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --vcf: {}. Must be one of {{clip, filter, full}}'.format(vcf))
        if vcf == 'clip' and not options.clip:
            raise RuntimError('--vcf cannot be set to clip since clipping is disabled')
        if vcf == 'filter' and not options.filter:
            raise RuntimeError('--vcf cannot be set to filter since filtering is disabled')        

    # if someone used --vcfReference they probably want --vcf too
    if options.vcfReference and not options.vcf:
        raise RuntimeError('--vcfReference cannot be used without --vcf')

    if options.vcfwave:
        if not options.vcf:
            raise RuntimeError('--vcfwave cannot be used without --vcf')

        if options.batchSystem == 'single_machine':
            options.vcfwaveCores = min(options.vcfwaveCores, cactus_cpu_count())
        
    if not options.vcfReference:
        options.vcfReference = [options.reference[0]]
    else:
        for vcfref in options.vcfReference:
            if vcfref not in options.reference:
                raise RuntimeError('--vcfReference {} invalid because {} was not specified as a --reference'.format(vcfref, vcfref))

    if options.snarlStats == []:
        options.snarlStats = ['clip'] if options.clip else ['full']
    options.snarlStats = list(set(options.snarlStats)) if options.snarlStats else []    
    for snarl_stats in options.snarlStats:
        if snarl_stats not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --snarlStats: {}. Must be one of {{clip, filter, full}}'.format(snarl_stats))
        if snarl_stats == 'clip' and not options.clip:
            raise RuntimError('--snarlStats cannot be set to clip since clipping is disabled')
        if snarl_stats == 'filter' and not options.filter:
            raise RuntimeError('--snarlStats cannot be set to filter since filtering is disabled')

    if options.giraffe == []:
        options.giraffe = ['filter'] if options.filter else ['clip'] if options.clip else ['full']
    options.giraffe = list(set(options.giraffe)) if options.giraffe else []        
    for giraffe in options.giraffe:
        if giraffe not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --giraffe: {}. Must be one of {{clip, filter, full}}'.format(giraffe))
        if giraffe == 'clip' and not options.clip:
            raise RuntimError('--giraffe cannot be set to clip since clipping is disabled')
        if giraffe == 'filter' and not options.filter:
            raise RuntimeError('--giraffe cannot be set to filter since filtering is disabled')

    if options.haplo == []:
        options.haplo = ['clip'] if options.clip else ['full']
    options.haplo = list(set(options.haplo)) if options.haplo else []        
    for haplo in options.haplo:
        if haplo not in ['clip', 'full']:
            raise RuntimeError('Unrecognized value for --haplo: {}. Must be one of {{clip, full}}'.format(haplo))
        if haplo == 'clip' and not options.clip:
            raise RuntimError('--haplo cannot be set to clip since clipping is disabled')
        # you always want a gbz with haplo
        if haplo not in options.gbz:
            logger.warning("Activating --gbz {} since --haplo {} was specified".format(haplo, haplo))
            options.gbz.append(haplo)
        
    # Prevent some useless compute due to default param combos
    if options.clip and 'clip' not in options.gfa + options.gbz + options.odgi + options.chrom_vg + options.chrom_og + options.vcf + options.giraffe + options.viz + options.draw\
       and 'filter' not in options.gfa + options.gbz + options.odgi + options.chrom_vg + options.chrom_og + options.vcf + options.giraffe + options.viz + options.draw:
        options.clip = None
    if options.filter and 'filter' not in options.gfa + options.gbz + options.odgi + options.chrom_vg + options.chrom_og + options.vcf + options.giraffe + options.viz + options.draw:
        options.filter = None

    return options
    
def graphmap_join(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            wf_output = toil.restart()
        else:
            options.cactusDir = getTempDirectory()

            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals(options)

            if options.collapse:
                findRequiredNode(configNode, "graphmap").attrib["collapse"] = 'all'
                
            # load up the vgs
            vg_ids = []
            for vg_path in options.vg:
                logger.info("Importing {}".format(vg_path))
                vg_ids.append(toil.importFile(makeURL(vg_path)))
                
            # load up the hals
            hal_ids = []
            for hal_path in options.hal:
                logger.info("Importing {}".format(hal_path))
                hal_ids.append(toil.importFile(makeURL(hal_path)))

            # run the workflow
            wf_output = toil.start(Job.wrapJobFn(graphmap_join_workflow, options, config, vg_ids, hal_ids))
                
        #export the split data
        export_join_data(toil, options, wf_output[0], wf_output[1], wf_output[2], wf_output[3], wf_output[4], wf_output[5])

def graphmap_join_workflow(job, options, config, vg_ids, hal_ids):

    root_job = Job()
    job.addChild(root_job)

    # vcfwave isn't included in the static binary release, so we start by checking it's available
    if options.vcfwave and options.vcf:
        vcfwave_check_job = root_job.addFollowOnJobFn(check_vcfwave)
        root_job = vcfwave_check_job
        
    # make sure reference doesn't have a haplotype suffix, as it will have been changed upstream
    ref_base, ref_ext = os.path.splitext(options.reference[0])
    assert len(ref_base) > 0
    if ref_ext:
        assert len(ref_ext) > 1 and ref_ext[1:].isnumeric()
        options.reference[0] = ref_base
        if options.vcfReference:
            for i in range(len(options.vcfReference)):
                if options.vcfReference[i] == ref_base + ref_ext:
                    options.vcfReference[i] = ref_base

    # get the single machine memory
    config.setSystemMemory(options)    

    # use --indexMemory to cap memory on some jobs
    max_mem = options.indexMemory if options.indexMemory else sys.maxsize
    if not options.indexMemory:
        max_mem = cactus_clamp_memory(max_mem)
    
    # keep og ids in separate structure indexed on chromosome
    og_chrom_ids = {'full' : defaultdict(list), 'clip' : defaultdict(list), 'filter' : defaultdict(list)}
    # have a lot of trouble getting something working for small contigs, hack here:
    og_min_size = max(2**31, int(sum(vg_id.size for vg_id in vg_ids) / len(vg_ids)) * 32)
    # run the "full" phase to do pre-clipping stuff
    full_vg_ids = []
    assert len(options.vg) == len(vg_ids)
    for vg_path, vg_id in zip(options.vg, vg_ids):
        full_job = Job.wrapJobFn(clip_vg, options, config, vg_path, vg_id, 'full',
                                 disk=vg_id.size * 20, memory=max(2**31, min(vg_id.size * 20, max_mem)))
        root_job.addChild(full_job)
        full_vg_ids.append(full_job.rv(0))
        if 'full' in options.odgi + options.chrom_og + options.viz + options.draw:
            full_og_job = full_job.addFollowOnJobFn(vg_to_og, options, config, vg_path, full_job.rv(0),
                                                    disk=vg_id.size * 16, memory=min(max(og_min_size, vg_id.size * 32), max_mem))
            og_chrom_ids['full']['og'].append(full_og_job.rv())

    prev_job = root_job
    
    # join the ids
    join_job = prev_job.addFollowOnJobFn(join_vg, options, config, full_vg_ids,
                                         disk=sum([f.size for f in vg_ids]),
                                         memory=min(max([f.size for f in vg_ids]) * 4, max_mem))
    full_vg_ids = [join_job.rv(i) for i in range(len(vg_ids))]
    prev_job = join_job

    # take out the _MINIGRAPH_ paths
    if 'full' in options.chrom_vg:
        output_full_vg_ids = []
        for vg_path, vg_id, full_vg_id in zip(options.vg, vg_ids, full_vg_ids):
            drop_graph_event_job = join_job.addFollowOnJobFn(drop_graph_event, config, vg_path, full_vg_id,
                                                             disk=vg_id.size * 3, memory=min(vg_id.size * 6, max_mem))
            output_full_vg_ids.append(drop_graph_event_job.rv())
    else:
        output_full_vg_ids = full_vg_ids

    # run the "clip" phase to do the clip-vg clipping
    clip_vg_ids = []
    clipped_stats = None 
    if options.clip or options.filter:
        clip_root_job = Job()
        prev_job.addFollowOn(clip_root_job)
        clip_vg_stats = []
        assert len(options.vg) == len(full_vg_ids) == len(vg_ids)
        for vg_path, vg_id, input_vg_id in zip(options.vg, full_vg_ids, vg_ids):
            clip_job = Job.wrapJobFn(clip_vg, options, config, vg_path, vg_id, 'clip',
                                     disk=input_vg_id.size * 20, memory=max(2**31, min(input_vg_id.size * 20, max_mem)))
            clip_root_job.addChild(clip_job)
            clip_vg_ids.append(clip_job.rv(0))
            clip_vg_stats.append(clip_job.rv(1))
            if 'clip' in options.odgi + options.chrom_og + options.viz + options.draw:
                clip_og_job = clip_job.addFollowOnJobFn(vg_to_og, options, config, vg_path, clip_job.rv(0),
                                                        disk=input_vg_id.size * 16,
                                                        memory=min(max(og_min_size, input_vg_id.size * 32), max_mem))
                og_chrom_ids['clip']['og'].append(clip_og_job.rv())

        # join the stats
        clipped_stats = clip_root_job.addFollowOnJobFn(cat_stats, clip_vg_stats).rv()
        prev_job = clip_root_job

    # run the "filter" phase to do the vg clip clipping
    filter_vg_ids = []
    if options.filter:
        filter_root_job = Job()
        prev_job.addFollowOn(filter_root_job)
        assert len(options.vg) == len(clip_vg_ids) == len(vg_ids)
        for vg_path, vg_id, input_vg_id in zip(options.vg, clip_vg_ids, vg_ids):
            filter_job = filter_root_job.addChildJobFn(vg_clip_vg, options, config, vg_path, vg_id, 
                                                       disk=input_vg_id.size * 20,
                                                       memory=max(2**31, min(input_vg_id.size * 22, max_mem)))
            filter_vg_ids.append(filter_job.rv())
            if 'filter' in options.odgi + options.chrom_og + options.viz + options.draw:
                filter_og_job = filter_job.addFollowOnJobFn(vg_to_og, options, config, vg_path, filter_job.rv(),
                                                            disk=input_vg_id.size * 16,
                                                            memory=min(max(og_min_size, input_vg_id.size * 64), max_mem))
                og_chrom_ids['filter']['og'].append(filter_og_job.rv())                


        prev_job = filter_root_job

    # set up our whole-genome output
    out_dicts = []

    # optional hal merge
    if hal_ids:
        hal_merge_job = job.addChildJobFn(merge_hal, options, hal_ids,
                                          cores = 1,
                                          disk=sum(f.size for f in hal_ids) * 2,
                                          memory=min(max(f.size for f in hal_ids) * 2, max_mem))
        hal_id_dict = hal_merge_job.rv()
        out_dicts.append(hal_id_dict)
        # delete the chromosome hals
        hal_merge_job.addFollowOnJobFn(clean_jobstore_files, file_ids=hal_ids)

    if options.indexMemory:
        index_mem = options.indexMemory
    else:
        # hacky heuristic to get memory for vg jobs, where vg minimizer is tricky as
        # it's not really based on graph size so much (ie simple graphs can take lots of mem)
        index_mem = sum(f.size for f in vg_ids)
        for tot_size, fac in zip([1e6, 1e9, 10e9, 50e9, 100e9, 500e9], [1000, 32, 10, 6, 5, 2]):
            if index_mem < tot_size:
                index_mem *= fac
                break
    max_system_memory = config.getSystemMemory()
    if max_system_memory and index_mem > max_system_memory:
        RealtimeLogger.info("Clamping index memory {} to system limit of {}".format(index_mem, max_system_memory))
        index_mem = max_system_memory

    # keep track of unclipped reference fastas for bcftools norm
    need_ref_fasta = (options.vcfbub and getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "bcftoolsNorm", typeFn=bool, default=False)) or\
        (options.vcfwave and getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "vcfwaveNorm", typeFn=bool, default=True))        
    ref_fasta_job = None
        
    workflow_phases = [('full', full_vg_ids, join_job)]
    if options.clip:
        workflow_phases.append(('clip', clip_vg_ids, clip_root_job))
    if options.filter:
        workflow_phases.append(('filter', filter_vg_ids, filter_root_job))        
    for workflow_phase, phase_vg_ids, phase_root_job in workflow_phases:
            
        # make a gfa for each
        gfa_root_job = Job()        
        phase_root_job.addFollowOn(gfa_root_job)
        gfa_ids = []
        current_out_dict = None
        do_gbz = workflow_phase in options.gbz + options.vcf + options.giraffe + options.xg
        if workflow_phase == 'full' and need_ref_fasta:
            do_gbz = True
        if do_gbz or workflow_phase in options.gfa:
            assert len(options.vg) == len(phase_vg_ids) == len(vg_ids)
            for vg_path, vg_id, input_vg_id in zip(options.vg, phase_vg_ids, vg_ids):
                gfa_job = gfa_root_job.addChildJobFn(vg_to_gfa, options, config, vg_path, vg_id,
                                                     disk=input_vg_id.size * 10,
                                                     memory=min(max(2**31, input_vg_id.size * 16), max_mem))
                gfa_ids.append(gfa_job.rv())

            gfa_merge_job = gfa_root_job.addFollowOnJobFn(make_vg_indexes, options, config, gfa_ids,
                                                          tag=workflow_phase + '.',
                                                          do_gbz=do_gbz,
                                                          cores=options.indexCores,
                                                          disk=sum(f.size for f in vg_ids) * 6,
                                                          memory=index_mem)
            out_dicts.append(gfa_merge_job.rv())
            prev_job = gfa_merge_job
            current_out_dict = gfa_merge_job.rv()
            if workflow_phase == 'full' and need_ref_fasta:
                # we need fasta references from the gbz for vcf normalization
                ref_fasta_job = gfa_merge_job.addFollowOnJobFn(extract_gbz_fasta, options, current_out_dict,
                                                               tag=workflow_phase + '.',
                                                               memory=index_mem,
                                                               disk=sum(f.size for f in vg_ids) * 2)
        # optional unchopped gfa
        if workflow_phase in options.unchopped_gfa:
            unchopped_gfa_ids = []
            for vg_path, vg_id, input_vg_id in zip(options.vg, phase_vg_ids, vg_ids):
                unchopped_gfa_job = gfa_root_job.addChildJobFn(vg_to_gfa, options, config, vg_path, vg_id, unchopped=True,
                                                               disk=input_vg_id.size * 10,
                                                               memory=min(max(2**31, input_vg_id.size * 16), max_mem))
                unchopped_gfa_ids.append(unchopped_gfa_job.rv())

            unchopped_gfa_merge_job = gfa_root_job.addFollowOnJobFn(make_vg_indexes, options, config, unchopped_gfa_ids,
                                                                    tag=workflow_phase + '.unchopped.',
                                                                    do_gbz=False,
                                                                    cores=1,
                                                                    disk=sum(f.size for f in vg_ids) * 3)
            out_dicts.append(unchopped_gfa_merge_job.rv())                
            

        # optional xg
        if workflow_phase in options.xg:
            xg_job = prev_job.addFollowOnJobFn(make_xg, config, options.outName, current_out_dict,
                                               tag=workflow_phase + '.',
                                               cores=max(options.indexCores, 4),
                                               disk = sum(f.size for f in vg_ids) * 10,
                                               memory=index_mem)
            out_dicts.append(xg_job.rv())

        # optional vcf
        if workflow_phase in options.vcf:
            for vcf_ref in options.vcfReference:
                vcf_job = gfa_root_job.addFollowOnJobFn(make_vcf, config, options, workflow_phase,
                                                        index_mem, vcf_ref, phase_vg_ids,
                                                        ref_fasta_job.rv() if ref_fasta_job else None)
                if ref_fasta_job:
                    ref_fasta_job.addFollowOn(vcf_job)

                out_dicts.append(vcf_job.rv())
                    
        # optional giraffe
        giraffe_job = None
        if workflow_phase in options.giraffe:
            giraffe_job = gfa_merge_job.addFollowOnJobFn(make_giraffe_indexes, options, config, current_out_dict,
                                                         tag=workflow_phase + '.',
                                                         cores=options.indexCores,
                                                         disk = sum(f.size for f in vg_ids) * 16,
                                                         memory=index_mem)
            out_dicts.append(giraffe_job.rv())

        # optional haplo index
        if workflow_phase in options.haplo:
            prec_job = giraffe_job if giraffe_job else gfa_merge_job
            haplo_job = prec_job.addFollowOnJobFn(make_haplo_index, options, config, current_out_dict,
                                                  giraffe_job.rv() if giraffe_job else None,
                                                  tag=workflow_phase + '.',
                                                  cores=options.indexCores,
                                                  disk = sum(f.size for f in vg_ids) * 16,
                                                  memory=index_mem)
            out_dicts.append(haplo_job.rv())

        # optional full-genome odgi
        if workflow_phase in options.odgi:
            odgi_job = gfa_root_job.addChildJobFn(odgi_squeeze, config, options.vg, og_chrom_ids[workflow_phase]['og'],
                                                  tag=workflow_phase + '.', disk=sum(f.size for f in vg_ids) *4,
                                                  memory=index_mem, cores=options.indexCores)
            out_dicts.append(odgi_job.rv())                                                  

        # optional viz
        other_contig = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "otherContigName", typeFn=str, default="chrOther")
        if workflow_phase in options.viz + options.draw:
            do_viz = workflow_phase in options.viz
            do_draw = workflow_phase in options.draw
            for vg_path, input_vg_id, og_id in zip(options.vg, vg_ids, og_chrom_ids[workflow_phase]['og']):
                # chrOther can contain tons of components, not sure if it ever makes sense to try to visualize
                if not os.path.basename(vg_path).startswith(other_contig + '.'):
                    viz_job = gfa_root_job.addChildJobFn(make_odgi_viz, config, options, vg_path, og_id, tag=workflow_phase,
                                                         viz=do_viz, draw=do_draw,                                                     
                                                         cores=options.indexCores, disk = input_vg_id.size * 10,
                                                         memory=min(max(og_min_size, input_vg_id.size * 32), max_mem))
                else:
                    viz_job = None
                if do_viz:
                    og_chrom_ids[workflow_phase]['viz'].append(viz_job.rv(0) if viz_job else None)
                if do_draw:
                    og_chrom_ids[workflow_phase]['draw'].append(viz_job.rv(1) if viz_job else None)

        # optional snarl stats table
        if workflow_phase in options.snarlStats:
            snarl_stats_ids = []
            for vg_path, vg_id, input_vg_id in zip(options.vg, phase_vg_ids, vg_ids):
                
                snarl_stats_job = gfa_root_job.addChildJobFn(snarl_stats, options, config, vg_path, vg_id,
                                                             disk=input_vg_id.size * 2,
                                                             memory=cactus_clamp_memory(input_vg_id.size * 5),
                                                             cores=options.indexCores)
                snarl_stats_ids.append(snarl_stats_job.rv())

            snarl_stats_merge_job = gfa_root_job.addFollowOnJobFn(merge_snarl_stats, options.vg, snarl_stats_ids,
                                                                  tag=workflow_phase + '.',
                                                                  disk=sum(f.size for f in vg_ids) * 2)
            out_dicts.append(snarl_stats_merge_job.rv())
    
    return output_full_vg_ids, clip_vg_ids, clipped_stats, filter_vg_ids, out_dicts, og_chrom_ids

def clip_vg(job, options, config, vg_path, vg_id, phase):
    """ run clip-vg 
    """
    assert phase in ['full', 'clip']
    
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)

    clipped_path = vg_path + '.clip'
    clipped_bed_path = vg_path + '.clip.bed'

    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")

    # optional gfaffixing -- only happens in full
    if phase == 'full' and getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "gfaffix", typeFn=bool, default=True):
        normalized_path = vg_path + '.gfaffixed'
        gfa_in_path = vg_path + '.gfa'
        gfa_out_path = normalized_path + '.gfa'
        # GFAffix's --dont_collapse option doesn't work on W-lines, so we strip them for now with -W
        cactus_call(parameters=['vg', 'convert', '-W', '-f', vg_path], outfile=gfa_in_path, job_memory=job.memory)
        fix_cmd = ['gfaffix', gfa_in_path, '--output_refined', gfa_out_path, '--check_transformation', '--threads', str(job.cores)]
        if options.reference:
            fix_cmd += ['--dont_collapse', options.reference[0] + '#[.]*']
        cactus_call(parameters=fix_cmd, job_memory=job.memory)
        # GFAffix strips the header, until this is fixed we need to add it back (for the RS tag)
        gfa_header = cactus_call(parameters=['head', '-1', gfa_in_path], check_output=True).strip()
        cactus_call(parameters=['sed', '-i', gfa_out_path, '-e', '1s/.*/{}/'.format(gfa_header)])
        # Come back from gfa to vg
        conv_cmd = [['vg', 'convert', '-g', '-p', gfa_out_path]]
        # GFAFfix doesn't unchop, so we do that in vg after
        conv_cmd.append(['vg', 'mod', '-u', '-'])
        cactus_call(parameters=conv_cmd, outfile=normalized_path)        
        vg_path = normalized_path

    # run clip-vg no matter what, but we don't actually remove anything in full
    cmd = []
    clip_vg_cmd = ['clip-vg', vg_path, '-f']

    # don't clip the (first) reference, also enforces it's in forward orientation (for rGFA)
    clip_vg_cmd += ['-e', options.reference[0]]

    # our vg file has minigraph sequences -- we'll filter them out, along with any nodes
    # that don't appear in a non-minigraph path
    clip_vg_cmd += ['-d', graph_event]
    if phase == 'full':
        # but... in the full graph we'll leave the minigraph path fragments that are aligned to anything else in the vg's
        # since we'll need them to run the actual clipping
        clip_vg_cmd += ['-L']

    if phase == 'clip':
        if options.clip:
            clip_vg_cmd += ['-u', str(options.clip)]
            if getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "clipNonMinigraph", typeFn=bool, default=True):
                clip_vg_cmd += ['-a', graph_event]
            clip_vg_cmd += ['-o', clipped_bed_path]

    # disable reference cycle check if desiired
    if getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "collapse", typeFn=str, default="none") in ["all", "reference"]:
        clip_vg_cmd += ['-c']

    cmd.append(clip_vg_cmd)

    # enforce chopping
    max_node_len = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "maxNodeLength", typeFn=int, default=-1)
    if phase == 'full':
        if max_node_len > 0:
            cmd.append(['vg', 'mod', '-X', str(max_node_len), '-'])
    
    if options.reference:
        if phase == 'full' and getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "pathNormalize", typeFn=bool, default=True):
            # run path normalization (important to do before vg clip -d1)
            # if node-chopping isn't enabled, we take the time to mod before/after since path normalization requires chopping
            if max_node_len <= 0 or max_node_len > 1024:
                cmd.append(['vg', 'mod', '-X', '1024', '-'])
            norm_command = ['vg', 'paths', '-x', '-', '-n', '-Q', options.reference[0], '-t', str(job.cores)]
            cmd.append(norm_command)
            if max_node_len <= 0 or max_node_len > 1024:
                cmd.append(['vg', 'mod', '-u', '-'])
                if max_node_len > 1024:
                    cmd.append(['vg', 'mod', '-X', str(max_node_len), '-'])
            
        # GFAFfix can leave uncovered nodes with --dont_collapse.  We filter out here so they dont cause trouble later
        # Also: any kind of clipping or path dropping can leave uncovered edges, so we remove them with vg clip        
        clip_cmd = ['vg', 'clip', '-d', '1', '-', '-P', options.reference[0]]
        cmd.append(clip_cmd)
        if phase == 'clip' and getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "removeStubs", typeFn=bool, default=True):
            # todo: could save a little time by making vg clip smart enough to do two things at once
            stub_cmd = ['vg', 'clip', '-sS', '-', '-P', options.reference[0]]

            # todo: do we want to add the minigraph prefix to keep stubs from minigraph? but I don't think it makes stubs....
            cmd.append(stub_cmd)

    # and we sort by id on the first go-around
    if phase == 'full':
        cmd.append(['vg', 'ids', '-s', '-'])
        
    cactus_call(parameters=cmd, outfile=clipped_path, job_memory=job.memory)

    # worth it
    cactus_call(parameters=['vg', 'validate', clipped_path])

    if phase == 'full':
        job.fileStore.deleteGlobalFile(vg_id)
        
    # keep some stats of the chromosomal vg file:
    # Path Lengths
    if phase == 'clip':
        chr_name = os.path.splitext(os.path.basename(vg_path))[0]
        path_stats_path = vg_path + '.path-stats.tsv'
        cactus_call(parameters=[['vg', 'paths', '-E', '-v', clipped_path],
                                ['awk', '{{print "{}\t" $0}}'.format(chr_name)]],
                    outfile=path_stats_path)
        # Nodes, edges and total length
        graph_stats_path = vg_path + '.graph-stats.tsv'
        cactus_call(parameters=[['vg', 'stats', '-l', '-z', clipped_path],
                                ['awk', '{{print "{}\t" $0}}'.format(chr_name)]],
                    outfile=graph_stats_path)
        # Stick the contig identifier onto the clipped regions bed
        clipped_bed_chr_path = clipped_bed_path + '.named'
        cactus_call(parameters=['awk',  '{{print $0 "\t{}"}}'.format(chr_name), clipped_bed_path],
                    outfile=clipped_bed_chr_path)
        sample_stats_path = vg_path + '.sample-stats.tsv'
        sample_stats = {}
        with open(path_stats_path, 'r') as path_stats_file:
            for line in path_stats_file:
                toks = line.split()
                contig_name = toks[1]
                contig_length = int(toks[2])
                sample_name = contig_name if '#' not in contig_name else '#'.join(contig_name.split('#')[:2])
                if sample_name not in sample_stats:
                    sample_stats[sample_name] = contig_length
                else:
                    sample_stats[sample_name] += contig_length
        with open(sample_stats_path, 'w') as sample_stats_file:
            for sample in sorted(sample_stats.keys()):
                sample_stats_file.write("{}\t{}\t{}\n".format(chr_name, sample, sample_stats[sample]))

        out_stats = { 'clip-stats.bed' : job.fileStore.writeGlobalFile(clipped_bed_chr_path),
                      'path-stats.tsv' : job.fileStore.writeGlobalFile(path_stats_path),
                      'graph-stats.tsv' : job.fileStore.writeGlobalFile(graph_stats_path),
                      'sample-stats.tsv' : job.fileStore.writeGlobalFile(sample_stats_path) }
    else:
        out_stats = None
    return job.fileStore.writeGlobalFile(clipped_path), out_stats

def vg_clip_vg(job , options, config, vg_path, vg_id):
    """ run vg clip, chaining multiple invocations if desired.
    """
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    clipped_path = vg_path + '.clip'

    clip_cmd = ['vg', 'clip', vg_path, '-d', str(options.filter)]
    for ref in options.reference:
        clip_cmd += ['-P', ref]
    min_fragment = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "minFilterFragment", typeFn=int, default=None)
    if min_fragment:
        clip_cmd += ['-m', str(min_fragment)]
        # the min-fragment option is kind of ill thought out, and works as post-processing.  therefore it can
        # actually leave uncovered nodes, we clean here but note that it can leave fragments smaller than min-fragment
        clean_cmd = ['vg', 'clip', '-d', '1', '-']
        for ref in options.reference:
            clean_cmd += ['-P', ref]
        clip_cmd = [clip_cmd, clean_cmd]

    if getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "removeStubs", typeFn=bool, default=True):
        # this command can also leave fragments smaller than min-fragment
        stub_cmd = ['vg', 'clip', '-sS', '-']
        if options.reference:
            stub_cmd += ['-P', options.reference[0]]
        if min_fragment:
            clip_cmd.append(stub_cmd)
        else:
            clip_cmd = [clip_cmd, stub_cmd]

    cactus_call(parameters=clip_cmd, outfile=clipped_path, job_memory=job.memory)

    # worth it
    cactus_call(parameters=['vg', 'validate', clipped_path])

    return job.fileStore.writeGlobalFile(clipped_path)
    
def join_vg(job, options, config, clipped_vg_ids):
    """ run vg ids -j
    """
    work_dir = job.fileStore.getLocalTempDir()
    vg_paths = []

    assert len(options.vg) == len(clipped_vg_ids)
    for vg_path, vg_id in zip(options.vg, clipped_vg_ids):
        vg_path = os.path.join(work_dir, os.path.basename(vg_path))
        job.fileStore.readGlobalFile(vg_id, vg_path, mutable=True)
        vg_paths.append(vg_path)
        
    cactus_call(parameters=['vg', 'ids', '-j'] + vg_paths, job_memory=job.memory)

    return [job.fileStore.writeGlobalFile(f) for f in vg_paths]

def drop_graph_event(job, config, vg_path, full_vg_id):
    """ take the _MINIGRAPH_ paths out of a chrom-vg full output graph """
    work_dir = job.fileStore.getLocalTempDir()
    full_vg_path = os.path.join(work_dir, os.path.splitext(os.path.basename(vg_path))[0]) + '.full.vg'
    job.fileStore.readGlobalFile(full_vg_id, full_vg_path)
    out_path = full_vg_path + '.drop'
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    
    cactus_call(parameters=['vg', 'paths', '-d', '-S', graph_event, '-x', full_vg_path], outfile=out_path)
    return job.fileStore.writeGlobalFile(out_path)
    
def vg_to_gfa(job, options, config, vg_path, vg_id, unchopped=False):
    """ run gfa conversion """
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    out_path = vg_path + '.unchopped.gfa' if unchopped else vg_path + '.gfa'

    input_path = '-' if unchopped else os.path.basename(vg_path)
    cmd = ['vg', 'convert', '-f', input_path]
    if getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "collapse", typeFn=str, default="none") not in ["all", "reference"]:
        cmd += ['-Q', options.reference[0], '-B']
    if unchopped:
        cmd = [['vg', 'mod', '-u', os.path.basename(vg_path)], cmd]
    
    cactus_call(parameters=cmd, outfile=out_path, work_dir=work_dir, job_memory=job.memory)

    return job.fileStore.writeGlobalFile(out_path)

def vg_to_og(job, options, config, vg_path, vg_id):
    """ run odgi conversion """
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    gfa_path = vg_path + '.gfa'
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    cactus_call(parameters=[['vg', 'convert', '-f', '-W', os.path.basename(vg_path)],
                            ['grep', '-v', '^P	{}'.format(graph_event)]],
                outfile=gfa_path, work_dir=work_dir, job_memory=job.memory)
    og_path = vg_path + '.og'
    cactus_call(parameters=['odgi', 'build', '-g', os.path.basename(gfa_path), '-o',
                            os.path.basename(og_path), '-t', str(job.cores)], work_dir=work_dir, job_memory=job.memory)
    return job.fileStore.writeGlobalFile(og_path)

def make_vg_indexes(job, options, config, gfa_ids, tag="", do_gbz=False):
    """ merge of the gfas, then make gbz / snarls / trans
    """
    work_dir = job.fileStore.getLocalTempDir()
    vg_paths = []
    merge_gfa_path = os.path.join(work_dir, '{}merged.gfa'.format(tag))

    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    
    # merge the gfas
    assert len(options.vg) == len(gfa_ids)
    for i, (vg_path, gfa_id) in enumerate(zip(options.vg, gfa_ids)):
        gfa_path = os.path.join(work_dir, os.path.basename(vg_path) +  '.gfa')
        job.fileStore.readGlobalFile(gfa_id, gfa_path, mutable=True)
        if i == 0:
            # make sure every --reference sample is in the GFA header RS tag (which won't be the case if one sample
            # is completely missing from the first graph, as vg may filter it out apparently)
            cmd = [['head', '-1', gfa_path], ['sed', '-e', '1s/{}//'.format(graph_event)]]
            gfa_header = cactus_call(parameters=cmd, check_output=True).strip().split('\t')
            for i in range(len(gfa_header)):
                if gfa_header[i].startswith('RS:Z:'):
                    header_refs = set(gfa_header[i][len('RS:Z:'):].split(' '))
                    for ref_sample in options.reference:
                        if ref_sample not in header_refs:
                            gfa_header[i] += ' ' + ref_sample
            with open(merge_gfa_path, 'w') as merge_gfa_file:
                merge_gfa_file.write('\t'.join(gfa_header) + '\n')
        # strip out header and minigraph paths
        cmd = ['grep', '-v', '^H\|^W	{}'.format(graph_event), gfa_path]
        cactus_call(parameters=cmd, outfile=merge_gfa_path, outappend=True, job_memory=job.memory)
        job.fileStore.deleteGlobalFile(gfa_id)

    # make the gbz
    if do_gbz:
        gbz_path = os.path.join(work_dir, '{}merged.gbz'.format(tag))
        gbz_cmd = ['vg', 'gbwt', '-G', merge_gfa_path, '--gbz-format', '-g', gbz_path]
        # sanity check to make sure the gbz didn't chop anything
        check_trans = 0 < getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "maxNodeLength", typeFn=int, default=-1) <= 1024
        if check_trans:
            trans_path = os.path.join(work_dir, '{}merged.gbz.trans'.format(tag))
            gbz_cmd += ['--translation', trans_path]        
        cactus_call(parameters=gbz_cmd, job_memory=job.memory)
        if check_trans:
            with open(trans_path, 'r') as trans_file:
                for line in trans_file:
                    node_id, trans_id = line.strip().split()[1:3]
                    if node_id != trans_id:
                        raise RuntimeError('GBZ translated {} to {}, but since chopping is enabled this should not have happened. Please report this bug!!'.format(node_id, trans_id))

    # zip the gfa
    cactus_call(parameters=['bgzip', merge_gfa_path, '--threads', str(job.cores)])
    gfa_path = merge_gfa_path + '.gz'

    # make the snarls
    if do_gbz:
        snarls_path = os.path.join(work_dir, '{}merged.snarls'.format(tag))
        cactus_call(parameters=['vg', 'snarls', gbz_path, '-T', '-t', str(job.cores)], outfile=snarls_path, job_memory=job.memory)

    out_dict = { '{}gfa.gz'.format(tag) : job.fileStore.writeGlobalFile(gfa_path) }
    if do_gbz:
        out_dict['{}gbz'.format(tag)] = job.fileStore.writeGlobalFile(gbz_path)
        out_dict['{}snarls'.format(tag)] =  job.fileStore.writeGlobalFile(snarls_path)
    return out_dict

def extract_gbz_fasta(job, options, index_dict, tag):
    """ get the reference fasta files from the unclipped gbz so they can be used with bcftools norm
    """
    work_dir = job.fileStore.getLocalTempDir()
    gbz_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.gbz')
    job.fileStore.readGlobalFile(index_dict['{}gbz'.format(tag)], gbz_path)

    out_dict = {}
    for vcf_ref in options.vcfReference:
        fa_ref_path = os.path.join(work_dir, vcf_ref + '.fa.gz')            
        cactus_call(parameters=[['vg', 'paths', '-x', gbz_path, '-S', vcf_ref, '-F'],
                                ['sed', '-e', 's/{}#.\+#//g'.format(vcf_ref)],
                                ['bgzip', '--threads', str(job.cores)]],
                    outfile=fa_ref_path)
        out_dict[vcf_ref] = job.fileStore.writeGlobalFile(fa_ref_path)

    return out_dict

def make_xg(job, config, out_name, index_dict, tag='', drop_haplotypes=False):
    """ make the xg from the gbz
    """
    work_dir = job.fileStore.getLocalTempDir()
    gbz_path = os.path.join(work_dir, tag + os.path.basename(out_name) + '.gbz')
    xg_path = os.path.join(work_dir, tag + os.path.basename(out_name) + '.xg')
    job.fileStore.readGlobalFile(index_dict['{}gbz'.format(tag)], gbz_path)

    # make the xg
    xg_cmd = ['vg', 'convert', '-x', gbz_path, '-t', str(job.cores)]
    if drop_haplotypes:
        xg_cmd += ['-H']
    cactus_call(parameters=xg_cmd, outfile=xg_path)

    return { '{}xg'.format(tag) : job.fileStore.writeGlobalFile(xg_path) }

def make_vcf(job, config, options, workflow_phase, index_mem, vcf_ref, vg_ids, ref_fasta_dict):
    """ make the raw vcf with deconstruct. optionally add in the bub and wave vcfs too
    this is done in parallel on each chrom .vg graph
    """
    root_job = Job()
    job.addChild(root_job)
    vcftag = vcf_ref + '.' + workflow_phase if vcf_ref != options.reference[0] else workflow_phase
    raw_vcf_tbi_ids, bub_vcf_tbi_ids, wave_vcf_tbi_ids = [], [], []
    for vg_path, vg_id, in zip(options.vg, vg_ids):
        deconstruct_job = root_job.addChildJobFn(deconstruct, config, options.outName,
                                                 vcf_ref, vg_id,
                                                 tag=os.path.splitext(os.path.basename(vg_path))[0] + '.' + vcftag + '.',
                                                 cores=options.indexCores,
                                                 disk = vg_id.size * 6,
                                                 memory=index_mem)

        raw_vcf_id, raw_tbi_id = deconstruct_job.rv(0), deconstruct_job.rv(1)
        raw_vcf_tbi_ids.append((raw_vcf_id, raw_tbi_id))

        if options.vcfbub:
            vcfbub_job = deconstruct_job.addFollowOnJobFn(vcfbub, config, options.outName, vcf_ref,
                                                                  raw_vcf_id, raw_tbi_id,
                                                                  options.vcfbub,
                                                                  ref_fasta_dict,
                                                                  tag=os.path.splitext(os.path.basename(vg_path))[0] + '.' + vcftag + '.',
                                                                  disk = vg_id.size * 6,
                                                                  memory=cactus_clamp_memory(vg_id.size * 2))
            bub_vcf_id, bub_tbi_id = vcfbub_job.rv(0), vcfbub_job.rv(1)
            bub_vcf_tbi_ids.append((bub_vcf_id, bub_tbi_id))

        if options.vcfwave:
            vcfwave_job = deconstruct_job.addFollowOnJobFn(vcfwave, config, options.outName, vcf_ref,
                                                           raw_vcf_id, raw_tbi_id,
                                                           options.vcfbub,
                                                           ref_fasta_dict,
                                                           tag=os.path.splitext(os.path.basename(vg_path))[0] + '.' + vcftag + '.',
                                                           cores=options.vcfwaveCores,
                                                           disk=vg_id.size * 6,
                                                           memory=cactus_clamp_memory(options.vcfwaveMemory))
            wave_vcf_id, wave_tbi_id = vcfwave_job.rv(0), vcfwave_job.rv(1)
            wave_vcf_tbi_ids.append((wave_vcf_id, wave_tbi_id))

    merge_vcf_job = root_job.addFollowOnJobFn(vcf_cat, raw_vcf_tbi_ids, vcftag + '.raw.',
                                                          disk = sum(f.size for f in vg_ids) * 16)
    out_dict = {'{}.raw.vcf.gz'.format(vcftag) : merge_vcf_job.rv(0),
                '{}.raw.vcf.gz.tbi'.format(vcftag) : merge_vcf_job.rv(1) }
    if bub_vcf_tbi_ids:
        merge_bub_job = root_job.addFollowOnJobFn(vcf_cat, bub_vcf_tbi_ids, vcftag + '.bub.',
                                                  disk = sum(f.size for f in vg_ids) * 16)
        out_dict['{}.vcf.gz'.format(vcftag)] = merge_bub_job.rv(0)
        out_dict['{}.vcf.gz.tbi'.format(vcftag)] = merge_bub_job.rv(1)
    if wave_vcf_tbi_ids:
        merge_wave_job = root_job.addFollowOnJobFn(vcf_cat, wave_vcf_tbi_ids, vcftag + '.wave.',
                                                   disk = sum(f.size for f in vg_ids) * 16)            
        out_dict['{}.wave.vcf.gz'.format(vcftag)] = merge_wave_job.rv(0)
        out_dict['{}.wave.vcf.gz.tbi'.format(vcftag)] = merge_wave_job.rv(1)
        
    return out_dict

def deconstruct(job, config, out_name, vcf_ref, vg_id, tag):
    """ make the raw vcf
    """ 
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(out_name) + '.' + tag + 'vg')    
    job.fileStore.readGlobalFile(vg_id, vg_path)

    # make the vcf
    vcf_path = os.path.join(work_dir, os.path.basename(out_name) + '.' + tag + 'raw.vcf.gz')
    decon_cmd = ['vg', 'deconstruct', vg_path, '-P', vcf_ref, '-C', '-a', '-t', str(job.cores)]
    cactus_call(parameters=[decon_cmd, ['bgzip', '--threads', str(job.cores)]], outfile=vcf_path, job_memory=job.memory)
    try:
        cactus_call(parameters=['tabix', '-p', 'vcf', vcf_path])
        tbi_path = vcf_path + '.tbi'
    except Exception as e:
        cactus_call(parameters=['bcftools', 'index', '-c', vcf_path])
        tbi_path = vcf_path + '.csi'        

    return job.fileStore.writeGlobalFile(vcf_path), job.fileStore.writeGlobalFile(tbi_path)

def vcfbub(job, config, out_name, vcf_ref, vcf_id, tbi_id, max_ref_allele, fasta_ref_dict, tag):
    """ make the vcfbub vcf
    """
    work_dir = job.fileStore.getLocalTempDir()
    vcf_path = os.path.join(work_dir, os.path.basename(out_name) + '.' + tag + 'raw.vcf.gz')
    job.fileStore.readGlobalFile(vcf_id, vcf_path)
    job.fileStore.readGlobalFile(tbi_id, vcf_path + '.tbi')

    # short circuit on empty file (note zcat -> head exits 141, so we can't use cactus_call)
    if int(subprocess.check_output('gzip -dc {} | grep -v ^# | head | wc -l'.format(vcf_path),
                                   shell=True).decode('utf-8').strip()) == 0:    
        return vcf_id, tbi_id 

    # run vcfbub
    vcfbub_path = os.path.join(work_dir, os.path.basename(out_name) + '.' + tag + 'bub.vcf.gz')
    assert max_ref_allele
    bub_cmd = [['vcfbub', '--input', vcf_path, '--max-ref-length', str(max_ref_allele), '--max-level', '0']]
    if getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "bcftoolsNorm", typeFn=bool, default=False):
        fa_ref_path = os.path.join(work_dir, tag + 'fa.gz')
        job.fileStore.readGlobalFile(fasta_ref_dict[vcf_ref], fa_ref_path)
        bub_cmd.append(['bcftools', 'norm', '-f', fa_ref_path])
        bub_cmd.append(['bcftools', 'norm', '-m', '+any'])
        bub_cmd.append(['bcftools', 'sort', '-T', os.path.join(work_dir, 'bcftools.XXXXXX')])
    bub_cmd.append(['bgzip', '--threads', str(job.cores)])
    cactus_call(parameters=bub_cmd, outfile=vcfbub_path)
    try:
        cactus_call(parameters=['tabix', '-p', 'vcf', vcfbub_path])
        tbi_path = vcfbub_path + '.tbi'
    except Exception as e:
        cactus_call(parameters=['bcftools', 'index', '-c', vcfbub_path])
        tbi_path = vcfbub_path + '.csi'        

    return job.fileStore.writeGlobalFile(vcfbub_path), job.fileStore.writeGlobalFile(tbi_path)

def check_vcfwave(job):
    """ check to make sure vcfwave is installed """
    try:
        cactus_call(parameters=['vcfwave', '-h'])
    except:
        raise RuntimeError('--vcfwave option used, but vcfwave tool not found in PATH. vcfwave is *not* included in the cactus binary release, but it is in the cactus Docker image. If you have Docker installed, you can try running again with --binariesMode docker. Or running your whole command with docker run. If you cannot use Docker, then you will need to build vcflib yourself before retrying: source code and details here: https://github.com/vcflib/vcflib')

def vcfwave(job, config, out_name, vcf_ref, vcf_id, tbi_id, max_ref_allele, fasta_ref_dict, tag):
    """ run vcfwave """
    work_dir = job.fileStore.getLocalTempDir()
    vcf_path = os.path.join(work_dir, os.path.basename(out_name) + '.' + vcf_ref + tag + 'raw.vcf.gz')
    job.fileStore.readGlobalFile(vcf_id, vcf_path)
    job.fileStore.readGlobalFile(tbi_id, vcf_path + '.tbi')

    # short circuit on empty file (note zcat -> head exits 141, so we can't use cactus_call)
    if int(subprocess.check_output('gzip -dc {} | grep -v ^# | head | wc -l'.format(vcf_path),
                                   shell=True).decode('utf-8').strip()) == 0:    
        return vcf_id, tbi_id 

    # run vcfbub and vcfwave, using original HPRC recipe
    vcfwave_path = os.path.join(work_dir, os.path.basename(out_name) + '.' + vcf_ref + tag + 'wave.vcf.gz')
    wave_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "vcfwaveOptions", typeFn=str, default=None)
    assert wave_opts
    bubwave_cmd = [['vcfbub', '--input', vcf_path, '-l', '0', '-a', str(max_ref_allele)],
                   ['bcftools', 'annotate', '-x', 'INFO/AT'],                   
                   ['vcfwave'] + wave_opts.split(' ') + ['-t', str(job.cores)]]

    if getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "vcfwaveNorm", typeFn=bool, default=True):
        fa_ref_path = os.path.join(work_dir, tag + 'fa.gz')
        job.fileStore.readGlobalFile(fasta_ref_dict[vcf_ref], fa_ref_path)
        bubwave_cmd.append(['bcftools', 'norm', '-f', fa_ref_path])
        bubwave_cmd.append(['bcftools', 'sort', '-T', os.path.join(work_dir, 'bcftools.XXXXXX')])
    bubwave_cmd.append(['bgzip'])

    cactus_call(parameters=bubwave_cmd, outfile=vcfwave_path)
    try:
        cactus_call(parameters=['tabix', '-p', 'vcf', vcfwave_path])
        tbi_path = vcfwave_path + '.tbi'
    except Exception as e:
        cactus_call(parameters=['bcftools', 'index', '-c', vcfwave_path])
        tbi_path = vcfwave_path + '.csi'

    return job.fileStore.writeGlobalFile(vcfwave_path), job.fileStore.writeGlobalFile(tbi_path)

def vcf_cat(job, vcf_tbi_ids, tag):
    """ concat some vcf files """
    work_dir = job.fileStore.getLocalTempDir()
    vcf_paths = []
    sample_sets = []
    all_sample_set = set()
    for i, (vcf_id, tbi_id) in enumerate(vcf_tbi_ids):
        vcf_path = os.path.join(work_dir, '{}.{}vcf.gz'.format(i, tag))
        job.fileStore.readGlobalFile(vcf_id, vcf_path)
        job.fileStore.readGlobalFile(tbi_id, vcf_path + '.tbi')
        vcf_paths.append(vcf_path)
        samples = cactus_call(parameters=['bcftools', 'query', '-l', vcf_path], check_output=True).strip().split('\n')
        all_sample_set.update(samples)
        sample_sets.append(set(samples))

    all_sample_list = sorted(list(all_sample_set))
    all_sample_list_path = os.path.join(work_dir, 'samples.txt')
    with open(all_sample_list_path, 'w') as sample_list_file:
        for sample in all_sample_list:
            sample_list_file.write(sample + '\n')

    # bcftools concat doesn't work unless every input file has the same samples
    # but his may not happen since deconstruct is run individually on the chroms
    # so we use bcftools merge to add any missing samples to each chromsome
    # then view -S to sort them
    # Note: This only works if the samples are alphabetically sorted coming
    # out of deconstruct (since files with missing samples aren't touched)
    for i, (vcf_path, sample_set) in enumerate(zip(vcf_paths, sample_sets)):
        missing_samples = all_sample_set - sample_set
        if missing_samples:            
            header = cactus_call(parameters=['bcftools', 'view', '-h', vcf_path, '-s', '', '--force-samples'],
                                 check_output=True).strip() + '\tFORMAT'
            for s in missing_samples:
                header += '\t{}'.format(s)
            header += '\n\n'
            header_path = vcf_path + '.missing-sample-header'
            with open(header_path, 'w') as header_file:
                header_file.write(header)
            cactus_call(parameters=['bgzip', header_path])
            cactus_call(parameters=['tabix', '-fp', 'vcf', header_path + '.gz'])
            updated_vcf_path = vcf_path + '.fix'
            cactus_call(parameters=[['bcftools', 'merge', vcf_path, header_path + '.gz'],
                                    ['bcftools', 'view', '-S', all_sample_list_path, '-O', 'z']],
                        outfile=updated_vcf_path)
            vcf_paths[i] = updated_vcf_path
        
    cat_vcf_path = os.path.join(work_dir, '{}vcf.gz'.format(tag))
    cactus_call(parameters=['bcftools', 'concat', '-O', 'z', '--threads', str(job.cores)] + \
                [os.path.basename(vcf_path) for vcf_path in vcf_paths],
                work_dir=work_dir, outfile=cat_vcf_path)

    try:
        cactus_call(parameters=['tabix', '-p', 'vcf', cat_vcf_path])
        tbi_path = cat_vcf_path + '.tbi'
    except Exception as e:
        cactus_call(parameters=['bcftools', 'index', '-c', cat_vcf_path])
        tbi_path = cat_vcf_path + '.csi'

    for vcf_id, tbi_id in vcf_tbi_ids:
        job.fileStore.deleteGlobalFile(vcf_id)
        job.fileStore.deleteGlobalFile(tbi_id)

    return job.fileStore.writeGlobalFile(cat_vcf_path), job.fileStore.writeGlobalFile(tbi_path)

def make_giraffe_indexes(job, options, config, index_dict, tag=''):
    """ make giraffe-specific indexes: distance and minimaer """
    work_dir = job.fileStore.getLocalTempDir()
    gbz_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.gbz')
    job.fileStore.readGlobalFile(index_dict['{}gbz'.format(tag)], gbz_path)

    # make the distance index
    dist_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.dist')
    cactus_call(parameters=['vg', 'index', '-t', str(job.cores), '-j', dist_path, gbz_path], job_memory=job.memory)

    # make the minimizer index
    min_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.min')
    min_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "minimizerOptions", default='').split()
    cactus_call(parameters=['vg', 'minimizer'] + min_opts + ['-t', str(job.cores), '-d', dist_path, '-o', min_path, gbz_path], job_memory=job.memory)

    return { '{}min'.format(tag) : job.fileStore.writeGlobalFile(min_path),
             '{}dist'.format(tag) : job.fileStore.writeGlobalFile(dist_path) }

def make_haplo_index(job, options, config, index_dict, giraffe_dict, tag=''):
    """ make new haplotype subsampling (.hapl) index.  this generally replaces the giraffe (.min .dist) indexes """
    work_dir = job.fileStore.getLocalTempDir()
    gbz_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.gbz')
    job.fileStore.readGlobalFile(index_dict['{}gbz'.format(tag)], gbz_path)

    # make or download the distance index
    dist_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.dist')
    if giraffe_dict:
        job.fileStore.readGlobalFile(giraffe_dict['{}dist'.format(tag)], dist_path)
    else:
        # note we are using --no-nested-distance here to reduce memory, and its all haplo sampling needs
        dist_path += '1'
        cactus_call(parameters=['vg', 'index', '-t', str(job.cores), '-j', dist_path, gbz_path, '--no-nested-distance'], job_memory=job.memory)
        
    # make the r-index
    # note: we don't bother adding it to the output_dict since it's only used for making the .hapl
    ri_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.ri')
    ri_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "rindexOptions", default='').split()
    cactus_call(parameters=['vg', 'gbwt'] + ri_opts + ['--num-threads', str(job.cores), '-r', ri_path, '-Z', gbz_path], job_memory=job.memory)

    # make the haplotype index
    hapl_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.hapl')
    hapl_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "haplOptions", default='').split()
    cactus_call(parameters=['vg', 'haplotypes'] + hapl_opts + ['-t', str(job.cores), '-H', hapl_path, '-d', dist_path, '-r', ri_path, gbz_path])

    return { '{}hapl'.format(tag) : job.fileStore.writeGlobalFile(hapl_path) }

def odgi_squeeze(job, config, vg_paths, og_ids, tag=''):
    """ combine chrom odgis into one big odgi """
    work_dir = job.fileStore.getLocalTempDir()
    merged_path = os.path.join(work_dir, tag + 'og')
    og_paths = []
    for vg_path, og_id in zip(vg_paths, og_ids):
        og_path = os.path.join(work_dir, os.path.basename(os.path.splitext(vg_path)[0]) + '.{}og'.format(tag))
        job.fileStore.readGlobalFile(og_id, og_path)
        og_paths.append(og_path)
    list_path = os.path.join(work_dir, '{}squeeze.input'.format(tag))
    with open(list_path, 'w') as list_file:
        for og_path in og_paths:
            # important to use relative path for docker
            list_file.write(os.path.basename(og_path) +'\n')
    cactus_call(parameters=['odgi', 'squeeze', '-f', list_path, '-o', merged_path, '-t', str(job.cores)],
                work_dir=work_dir, job_memory=job.memory)
    return { '{}og'.format(tag) : job.fileStore.writeGlobalFile(merged_path) }    

        
def make_odgi_viz(job, config, options, vg_path, og_id, tag='', viz=True, draw=True):
    """ use odgi viz and draw to make some images """
    work_dir = job.fileStore.getLocalTempDir()
    og_path = os.path.join(work_dir, os.path.basename(os.path.splitext(vg_path)[0]) + '.{}.og'.format(tag))
    job.fileStore.readGlobalFile(og_id, og_path)

    og_sort_path = og_path + '.sort'
    cactus_call(parameters=['odgi', 'sort', '-i', og_path, '-o', og_sort_path, '-t', str(job.cores)], job_memory=job.memory)

    if viz:     
        # determine prefixes to merge together (ie sample name + haplotype, splitting on #)
        prefix_path = os.path.join(work_dir, 'path_sample_names')
        odgi_paths_output = cactus_call(parameters=['odgi', 'paths', '-i', og_sort_path, '-L'], check_output=True).strip()
        prefixes = set()
        for line in odgi_paths_output.split('\n'):
            m = re.match('.+#[0-9]+#', line)
            if m and m.span()[0] == 0:
                cut_pos = m.span()[1]
            else:
                cut_pos = line.find('#')
            if cut_pos > 0:
                prefixes.add(line[0:cut_pos])

        with open(prefix_path, 'w') as prefix_file:
            for prefix in sorted(list(prefixes)):
                prefix_file.write(prefix + '\n')

        viz_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "odgiVizOptions", default='').split()
        viz_png_path = og_path + '.viz.png'
        cactus_call(parameters=['odgi', 'viz', '-i', og_sort_path, '-o', viz_png_path, '-M', prefix_path, '-t', str(job.cores)] + viz_opts, job_memory=job.memory)
        
    if draw:    
        lay_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "odgiLayoutOptions", default='').split()
        lay_path = og_path + '.lay'
        cactus_call(parameters=['odgi', 'layout', '-i', og_sort_path, '-o', lay_path, '-t', str(job.cores)] + lay_opts, job_memory=job.memory)
    
        draw_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "odgiDrawOptions", default='').split()
        draw_png_path = og_path + '.draw.png'
        cactus_call(parameters=['odgi', 'draw', '-i', og_sort_path, '-c', lay_path, '-p', draw_png_path, '-t', str(job.cores)] + draw_opts, job_memory=job.memory)
    
    return (job.fileStore.writeGlobalFile(viz_png_path) if viz else None,
            job.fileStore.writeGlobalFile(draw_png_path) if draw else None)
    

def snarl_stats(job, options, config, vg_path, vg_id):
    """ use vg stats to make a table of snarls, including reference path intervals where possible"""
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)

    snarl_stats_path = vg_path + '.snarl-stats.tsv'
    cactus_call(parameters=['vg', 'stats', '-R', vg_path, '--snarl-sample', options.reference[0], '--threads', str(job.cores)],
                outfile=snarl_stats_path)

    return job.fileStore.writeGlobalFile(snarl_stats_path)

def merge_snarl_stats(job, vg_paths, snarl_stats_ids, tag):
    """ sort (decreasing reference interval sizer) and merge the different stats files """
    work_dir = job.fileStore.getLocalTempDir()
    local_paths = []
    for in_path, in_id in zip(vg_paths, snarl_stats_ids):
        name = os.path.splitext(os.path.basename(in_path))[0] + '.' + tag + 'ss.tsv'
        vg_path = os.path.join(work_dir, name)
        job.fileStore.readGlobalFile(in_id, vg_path)
        local_paths.append(os.path.basename(vg_path))

    merged_stats_path = 'merge{}.snarl-stats.tsv.gz'.format(tag)
    # get the header
    cactus_call(parameters=['head', '-1', local_paths[0]], outfile=merged_stats_path, work_dir=work_dir)

    # sort by interval size
    os.makedirs(os.path.join(work_dir, 'sorttmp'))
    cactus_call(parameters=[['cat'] + local_paths,
                            ['grep', '-v', '^#'],
                            ['awk', '{print $3-$2 \"\\t\" $0}'],
                            ['sort', '-k1', '-r', '-n', '-T', os.path.join(work_dir, 'sorttmp')],
                            ['cut', '-f2-'],
                            ['bgzip']],
                outfile=merged_stats_path,
                outappend=True,
                work_dir=work_dir)

    return { '{}snarl-stats.tsv.gz'.format(tag) : job.fileStore.writeGlobalFile(merged_stats_path) }

def merge_hal(job, options, hal_ids):
    """ call halMergeChroms to make one big hal file out of the chromosome hal files """
    work_dir = job.fileStore.getLocalTempDir()
    hal_paths = []
    assert len(options.hal) == len(hal_ids)
    for in_path, hal_id in zip(options.hal, hal_ids):
        hal_path = os.path.join(work_dir, os.path.basename(in_path))
        job.fileStore.readGlobalFile(hal_id, hal_path)
        hal_paths.append(hal_path)

    merged_path = os.path.join(work_dir, '__merged__.hal')
    assert merged_path not in hal_paths

    # note: cactus_call tries to sort out relative paths by itself for docker.  but the comma-separated list
    # will likely throw it off, so we take care to specify it relative manually.
    # also note: most hal commands need --inMemory to run at scale, but the access patterns for chrom
    # merging are linear enough that it shouldn't be needed
    cmd = ['halMergeChroms',
           ','.join([os.path.basename(p) for p in hal_paths]),
           os.path.basename(merged_path),
           '--progress']
    cactus_call(parameters=cmd, work_dir = work_dir, job_memory=job.memory)

    return { 'full.hal' : job.fileStore.writeGlobalFile(merged_path) }

def unclip_hal(job, hal_id_dict, seq_id_map):
    """ run halUnclip """
    work_dir = job.fileStore.getLocalTempDir()
    in_hal_path = os.path.join(work_dir, "in.hal")
    job.fileStore.readGlobalFile(hal_id_dict['hal'], in_hal_path)
    out_hal_path = os.path.join(work_dir, "out.hal")
    seqfile_path = os.path.join(work_dir, "seqfile.txt")
    with open(seqfile_path, "w") as seqfile:
        for event, name_id in seq_id_map.items():
            seqfile.write("{}\t{}\n".format(event, job.fileStore.readGlobalFile(name_id[1])))

    cactus_call(parameters=["halUnclip", in_hal_path, seqfile_path, out_hal_path, "--progress", "--validate"])

    return { 'hal' : job.fileStore.writeGlobalFile(out_hal_path) }

def cat_bed_files(job, bed_ids):
    bed_paths = [job.fileStore.readGlobalFile(bed_id) for bed_id in bed_ids]
    cat_bed_path = job.fileStore.getLocalTempFile()
    catFiles(bed_paths, cat_bed_path)
    renamed_bed_path = cat_bed_path + '.renamed'
    # cactus-preprocess will make a bed like "ID=SAMPLE.1|CONTIG"
    # but by the name its been through cactus/hal it'll be SAMPLE.1.CONTIG
    # so we need to make the same transformation in the bed
    # if the id=| prefix is not present, it's up to the user to make sure it's consistent on input
    with open(cat_bed_path, 'r') as ifile, open(renamed_bed_path, 'w') as ofile:
        for line in ifile:
            if line.startswith('ID='):
                p = line.find('|')
                if p > 3:
                    sample = line[3:p]
                    rest = line[p+1:]
                    line = '{}.{}'.format(sample, rest)
            ofile.write(line)    
    return job.fileStore.writeGlobalFile(renamed_bed_path)

def cat_stats(job, stats_dict_list, zip_stats=True):
    work_dir = job.fileStore.getLocalTempDir()
    merged_dict = {}
    for key in stats_dict_list[0].keys():
        merged_dict[key] = os.path.join(work_dir, key)
        catFiles([job.fileStore.readGlobalFile(sd[key]) for sd in stats_dict_list], merged_dict[key])
    if zip_stats:
        cactus_call(parameters=['tar', 'czf', os.path.join(work_dir, 'stats.tgz')] + list(merged_dict.values()))
        return { 'stats.tgz': job.fileStore.writeGlobalFile(os.path.join(work_dir, 'stats.tgz')) }
    else:
        for key in stats_dict_list[0].keys():
            merged_dict[key] = job.fileStore.writeGlobalFile(merged_dict[key])
        return merged_dict

def export_join_data(toil, options, full_ids, clip_ids, clip_stats, filter_ids, idx_maps, og_chrom_ids):
    """ download all the output data
    """

    # make a directory for the chromosomal vgs
    if options.chrom_vg:
        clip_base = os.path.join(options.outDir, '{}.chroms'.format(options.outName))
        if not clip_base.startswith('s3://') and not os.path.isdir(clip_base):
            os.makedirs(clip_base)

        if 'full' in options.chrom_vg:
            # download the "full" vgs
            assert len(options.vg) == len(full_ids)
            for vg_path, full_id,  in zip(options.vg, full_ids):
                name = os.path.splitext(vg_path)[0] + '.full.vg'
                toil.exportFile(full_id, makeURL(os.path.join(clip_base, os.path.basename(name))))

        if 'clip' in options.chrom_vg:
            # download the "clip" vgs
            assert len(options.vg) == len(clip_ids)
            for vg_path, clip_id,  in zip(options.vg, clip_ids):
                name = os.path.splitext(vg_path)[0] + '.vg'
                toil.exportFile(clip_id, makeURL(os.path.join(clip_base, os.path.basename(name))))

        if 'filter' in options.chrom_vg:
            # download the "filter" vgs
            assert len(options.vg) == len(filter_ids)
            for vg_path, filter_id,  in zip(options.vg, filter_ids):
                name = os.path.splitext(vg_path)[0] + '.d{}.vg'.format(options.filter)
                toil.exportFile(filter_id, makeURL(os.path.join(clip_base, os.path.basename(name))))

    # make a directory for the chromosomal ogs
    if options.chrom_og:
        clip_base = os.path.join(options.outDir, '{}.chroms'.format(options.outName))
        if not clip_base.startswith('s3://') and not os.path.isdir(clip_base):
            os.makedirs(clip_base)

        for gtype in options.chrom_og:
            # download all the chromosomal ogs
            tag = 'd{}'.format(options.filter) if gtype == 'filter' else gtype
            og_ids = og_chrom_ids[gtype]['og']
            assert len(options.vg) == len(og_ids)
            for vg_path, og_id in zip(options.vg, og_ids):
                name = os.path.splitext(vg_path)[0] + '{}.og'.format( '.' + tag if tag != 'clip' else '')
                toil.exportFile(og_id, makeURL(os.path.join(clip_base, os.path.basename(name))))

    # make a directory for the viz
    if options.viz + options.draw:
        viz_base = os.path.join(options.outDir, '{}.viz'.format(options.outName))
        if not viz_base.startswith('s3://') and not os.path.isdir(viz_base):
            os.makedirs(viz_base)

        for gtype in options.viz:
            # download all the chromosomal 1D visualizations
            assert len(options.vg) == len(og_chrom_ids[gtype]['viz'])
            tag = 'd{}'.format(options.filter) if gtype == 'filter' else gtype            
            for vg_path, viz_id in zip(options.vg, og_chrom_ids[gtype]['viz']):
                if viz_id:
                    viz_name = os.path.splitext(vg_path)[0] + '{}.viz.png'.format('.' + tag if tag != 'clip' else '')
                    toil.exportFile(viz_id, makeURL(os.path.join(viz_base, os.path.basename(viz_name))))

        for gtype in options.draw:
            # download all the chromosomal 2D visualizations
            assert len(options.vg) == len(og_chrom_ids[gtype]['draw'])
            tag = 'd{}'.format(options.filter) if gtype == 'filter' else gtype            
            for vg_path, draw_id in zip(options.vg, og_chrom_ids[gtype]['draw']):
                if draw_id:
                    draw_name = os.path.splitext(vg_path)[0] + '{}.draw.png'.format('.' + tag if tag != 'clip' else '')
                    toil.exportFile(draw_id, makeURL(os.path.join(viz_base, os.path.basename(draw_name))))
                
    # download the stats files
    if clip_stats:
        for stats_file in clip_stats.keys():
            toil.exportFile(clip_stats[stats_file], makeURL(os.path.join(options.outDir, '{}.{}'.format(options.outName, stats_file))))
        
    # download everything else
    for idx_map in idx_maps:
        for ext, idx_id in idx_map.items():
            # hacky filtering of intermediate indexes that the user doesn't want
            # ex if someone did --vcf clip --gfa full this would filter out clip.gfa etc.
            is_intermediate = False
            for out_ext, out_sel in [('gfa', options.gfa), ('gbz', options.gbz + options.giraffe), ('snarls', options.gbz)]:
                for phase in ['full', 'clip', 'filter']:
                    if '{}.{}'.format(phase, out_ext) in ext and phase not in out_sel:
                        is_intermediate = True
            if not is_intermediate:
                out_ext = ext
                if 'clip.' in out_ext:
                    out_ext = out_ext.replace('clip.', '')
                if 'filter.' in out_ext:
                    out_ext = out_ext.replace('filter.', 'd{}.'.format(options.filter))
                toil.exportFile(idx_id, makeURL(os.path.join(options.outDir, '{}.{}'.format(options.outName, out_ext))))
        
if __name__ == "__main__":
    main()
