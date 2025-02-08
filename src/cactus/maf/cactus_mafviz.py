#!/usr/bin/env python3

"""
This is a script to view a MAF (from cactus-hal2maf) using taffy
"""

import os, sys
from argparse import ArgumentParser
import copy
import timeit, time
import math
import xml.etree.ElementTree as ET

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.version import cactus_commit
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from cactus.maf.cactus_hal2maf import get_sed_rename_scripts
from toil.lib.humanize import bytes2human
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("mafFile", help = "MAF file to convert to BigMaf (can be gzipped)")
    parser.add_argument("outDir", help = "Output bigMaf file (.bb)")

    parser.add_argument("--refGenome", required=True, help = "Name of reference genome (1st row in MAF)")
    parser.add_argument("--chroms", nargs='*',
                        help = "Subset to input chromosome(s). If used without arguments, split on all chromosomes")

    parser.add_argument("--genomes", nargs='+',
                        help = "Split into separate output for each given genome.  Ex. --genomes $(halStats --genomes <hal>)")

    parser.add_argument("--haplotypes", nargs='+',
                        help = "Split into separate output for each given sample.  Same as --genomes except "
                        " samples will be treaded together. ex hap.1,hap.2,hap.3 will be in the same output file "
                        "and renamed to hap_1,hap_2,hap_3")
    
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

    if not (options.mafFile.endswith('.maf') or options.mafFile.endswith('.maf.gz') or
            options.mafFile.endswith('.taf') or options.mafFile.endswith('.taf.gz')):
        raise RuntimeError('Input MAF file must have .maf/.maf.gz/.taf/.taf.gz extension')
    
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            viz_id_dict = toil.restart()
        else:
            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals(options)
            
            logger.info("Importing {}".format(options.mafFile))
            maf_id = toil.importFile(makeURL(options.mafFile))
            try:
                logger.info("Importing {}".format(options.mafFile + '.tai'))
                tai_id = toil.importFile(makeURL(options.mafFile + '.tai'))
            except:
                raise RuntimeError('Maf Index {} not found: Please create with taffy index'.format(
                    options.mafFile + '.tai'))

            viz_id_dict = toil.start(Job.wrapJobFn(mafviz_workflow, config, options, maf_id, tai_id))

        #export the results
        export_mafviz(toil, options, viz_id_dict)
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-mafviz has finished after {} seconds".format(run_time))


def mafviz_workflow(job, config, options, maf_id, tai_id):
    root_job = Job()
    job.addChild(root_job)
    # get the reference chromosomes
    if options.chroms == []:
        chrom_job = root_job.addFollowOnJobFn(tai_ref_croms, options, maf_id, tai_id)
        root_job = chrom_job
        options.chroms = chrom_job.rv()

    viz_all_job = root_job.addFollowOnJobFn(maf_viz_all, config, options, maf_id, tai_id)

    return viz_all_job.rv()

def tai_ref_chroms(job, options, maf_id, tai_id):
    """ use taffy to get the reference chromosomes from the tai """
    work_dir = job.fileStore.getLocalTempDir()
    maf_path = os.path.join(work_dir, os.path.basename(options.mafFile))
    tai_path = maf_path + '.tai'
    job.fileStore.readGlobalFile(maf_id, maf_path)
    job.fileStore.readGlobalFile(tai_id, tai_path)
    seq_stats = cactus_call(parameters=['taffy' 'stats', '-s', '-i', maf_path], check_output=True)
    seq_list = []
    for line in seq_stats.split('\n'):
        seq_list.append(line.split()[0])
    return seq_list

def maf_viz_all(job, config, options, maf_id, tai_id):
    """ make visualization, for each chrom / genome as specified, in parallel """
    out_dict = {}

    # "None" chromsome means whole-genome
    if not options.chroms:
        options.chroms = [None]
    if options.chroms:
        for chrom in options.chroms:
            chrom_viz_job = job.addChildJobFn(maf_viz_chrom, config, options, maf_id, tai_id, chrom,
                                              disk=maf_id.size*2)
            out_dict[chrom] = chrom_viz_job.rv()

    return out_dict

def maf_viz_chrom(job, config, options, maf_id, tai_id, chrom):
    """ make visualizations for given chromosome """
    work_dir = job.fileStore.getLocalTempDir()
    maf_path = os.path.join(work_dir, os.path.basename(options.mafFile))
    tai_path = maf_path + '.tai'
    job.fileStore.readGlobalFile(maf_id, maf_path)
    job.fileStore.readGlobalFile(tai_id, tai_path)

    # extract the chromosome if specified
    if chrom:
        chrom_maf_path = os.path.join(work_dir, chrom + '.' + os.path.basename(options.mafFile))
        maf_view_command = ['taffy', 'view', '-m', '-i', maf_path, '-r', options.refGenome + '.' + chrom]
        if maf_path.endswith('.gz'):
            maf_view_command += ['-c']
        cactus_call(parameters=maf_view_command, outfile=chrom_maf_path)
        maf_path = chrom_maf_path
        maf_id = job.fileStore.writeGlobalFile(maf_path)

    if not options.genomes and not options.haplotypes:
        options.genomes = [None]

    genome_groups = []
    if options.genomes:
        for genome in options.genomes:
            genome_groups.append([genome])

    # group haplotypes by sample, using .1 .2 etc naming convention from mc
    if options.haplotypes:
        sample_map = {}
        for h in options.haplotypes:
            assert '.' in h
            sample = h[0:h.rfind('.')]
            if sample in sample_map:
                sample_map[sample].append(h)
            else:
                sample_map[sample] = [h]
        for hap_list in sample_map.values():
            genome_groups.append(hap_list)

    print(genome_groups)
            
    genome_dict = {}
    for genome_group in genome_groups:
        viz_job = job.addChildJobFn(maf_viz, config, options, maf_id, chrom, genome_group,
                                    disk=maf_id.size*3)
        genome_dict[sorted(genome_group)[0]] = viz_job.rv()

    print(genome_dict)
    return genome_dict

def maf_viz(job, config, options, maf_id, chrom, genomes):
    """ make the visualizations for the given maf and genomes """

    work_dir = job.fileStore.getLocalTempDir()
    prefix = chrom + '.' if chrom else ''
    maf_path = os.path.join(work_dir, prefix + os.path.basename(options.mafFile))
    job.fileStore.readGlobalFile(maf_id, maf_path)

    if genomes != [None]:
        sed_cmd = ['sed']
        extract_cmd = [['gzip', '-dc', maf_path]] if options.mafFile.endswith('.gz') else [['cat', maf_path]]
        new_genomes = []
        for genome in genomes:
            if '.' in genome:
                new_genome = '_'.join(genome.rsplit('.', 1))
                sed_cmd += ['-e', 's/{}/{}/g'.format(genome, new_genome)]
            else:
                new_genome = genome
            new_genomes.append(new_genome)
        if len(sed_cmd) > 1:
            extract_cmd.append(sed_cmd)
        extract_cmd.append(['mafFilter', '-m', '-', '-i', ','.join(new_genomes + [options.refGenome])])
        extract_cmd.append(['mafFilter', '-m', '-', '-l', '2'])
        
        extract_cmd.append(['bgzip'])

        extract_maf_path = os.path.join(work_dir, prefix + genomes[0] + '.' + os.path.basename(options.mafFile))
        if not extract_maf_path.endswith('.gz'):
            extract_maf_path += '.gz'
        cactus_call(parameters=extract_cmd, outfile=extract_maf_path)
    else:
        extract_maf_path = maf_path

    # now we can do this viz
    dotfile = os.path.join(work_dir, 'dotplot.pdf')
    cactus_call(parameters=['python3', '/home/hickey/dev/taffy/scripts/taffy_alignment_plot.py',
                            extract_maf_path, '--show_dot_plot', '--output_format', 'png',
                            '--out_file', dotfile])

    syntfile = os.path.join(work_dir, 'synteny.pdf')
    cactus_call(parameters=['python3', '/home/hickey/dev/taffy/scripts/taffy_alignment_plot.py',
                            extract_maf_path, '--show_synteny', '--output_format', 'png',
                            '--bin_number', '500', '--show_sequence_boundaries',
                            '--out_file', syntfile])

    return { 'dot' : job.fileStore.writeGlobalFile(dotfile),
             'synteny' : job.fileStore.writeGlobalFile(syntfile) }
                            
def export_mafviz(toil, options, output_dict):
    """ export the files into a directory """

    name = os.path.splitext(os.path.basename(options.mafFile))[0]
    if options.mafFile.endswith('.gz'):
        name = os.path.splitext(name)[0]

    for chrom, chrom_dict in output_dict.items():
        chrom_base = os.path.join(options.outDir, chrom) if chrom else options.outDir
        if not chrom_base.startswith('s3://') and not os.path.isdir(chrom_base):
            os.makedirs(chrom_base)

        for genome, genome_dict in chrom_dict.items():
            for viz_type, viz_id in genome_dict.items():
                file_name = name
                if chrom:
                    file_name += '.{}'.format(chrom)
                if genome:
                    file_name += '.{}'.format(genome)
                file_name += '.{}.png'.format(viz_type)
                
                toil.exportFile(viz_id, makeURL(os.path.join(chrom_base, file_name)))

