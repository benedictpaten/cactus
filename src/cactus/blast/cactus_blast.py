#!/usr/bin/env python3

#Released under the MIT license, see LICENSE.txt

"""Run the pairwise blast stage only

"""
import os
import sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import timeit

from cactus.progressive.progressive_decomposition import compute_outgroups, parse_seqfile, get_subtree, get_spanning_subtree, get_event_set
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.version import cactus_commit
from cactus.progressive.cactus_prepare import human2bytesN

from cactus.paf.local_alignment import sanitize_then_make_paf_alignments

from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory, getTempFile

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("seqFile", help="Seq file")
    parser.add_argument("outputFile", type=str, help="Output pairwise alignment file")
    parser.add_argument("--pathOverrides", nargs="*", help="paths (multiple allowed) to override from seqFile")
    parser.add_argument("--pathOverrideNames", nargs="*", help="names (must be same number as --pathOverrides) of path overrides")

    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--root", dest="root", help="Name of ancestral node (which"
                        " must appear in NEWICK tree in <seqfile>) to use as a "
                        "root for the alignment.  Any genomes not below this node "
                        "in the tree may be used as outgroups but will never appear"
                        " in the output.  If no root is specifed then the root"
                        " of the tree is used. ", default=None, required=True)
    parser.add_argument("--includeRoot", action="store_true", help="Include the root's sequence in the alignment"
                        " (used only when running alignment update recipes)")
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)
    parser.add_argument("--gpu", nargs='?', const='all', default=None, help="toggle on GPU-enabled lastz, and specify number of GPUs (all available if no value provided)")
    parser.add_argument("--lastzCores", type=int, default=None, help="Number of cores for each lastz/segalign job, only relevant when running with --gpu")
    parser.add_argument("--lastzMemory", type=human2bytesN,
                        help="Memory in bytes for each lastz/segalign job (defaults to an estimate based on the input data size). "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))", default=None)
    parser.add_argument("--chromInfo",
                        help="Two-column file mapping genome (col 1) to comma-separated list of sex chromosomes. This information "
                        "will be used to guide outgroup selection so that, where possible, all chromosomes are present in"
                        " at least one outgroup.")    
    
    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    if (options.pathOverrides or options.pathOverrideNames):
        if not options.pathOverrides or not options.pathOverrideNames or \
           len(options.pathOverrideNames) != len(options.pathOverrides):
            raise RuntimeError('same number of values must be passed to --pathOverrides and --pathOverrideNames')

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()
    runCactusBlastOnly(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-blast has finished after {} seconds".format(run_time))

def runCactusBlastOnly(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            paf_id = toil.restart()
        else:

            # load up the seqfile and figure out the outgroups and schedule
            config_node = ET.parse(options.configFile).getroot()
            config_wrapper = ConfigWrapper(config_node)
            config_wrapper.substituteAllPredefinedConstantsWithLiterals(options)
            config_wrapper.applySlurmChunkScaling(options)
            # apply gpu override
            config_wrapper.initLastz(options)
            mc_tree, input_seq_map, og_candidates = parse_seqfile(options.seqFile, config_wrapper)
            og_map = compute_outgroups(mc_tree, config_wrapper, set(og_candidates), chrom_info_file = options.chromInfo)
            event_set = get_event_set(mc_tree, config_wrapper, og_map, options.root)
            if options.includeRoot:
                event_set.add(options.root)
            
            # apply path overrides.  this was necessary for wdl which doesn't take kindly to
            # text files of local paths (ie seqfile).  one way to fix would be to add support
            # for s3 paths and force wdl to use it.  a better way would be a more fundamental
            # interface shift away from files of paths throughout all of cactus
            if options.pathOverrides:
                for name, override in zip(options.pathOverrideNames, options.pathOverrides):
                    input_seq_map[name] = override

            # get the spanning tree (which is what it paf aligner wants)
            spanning_tree = get_spanning_subtree(mc_tree, options.root, config_wrapper, og_map)
                    
            #import the sequences
            input_seq_id_map = {}
            for (genome, seq) in input_seq_map.items():
                if genome in event_set:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    logger.info("Importing {}".format(seq))
                    input_seq_id_map[genome] = toil.importFile(seq)

            paf_id = toil.start(Job.wrapJobFn(sanitize_then_make_paf_alignments, NXNewick().writeString(spanning_tree),
                                              input_seq_id_map, options.root, config_node))

        # export the alignments
        toil.exportFile(paf_id, makeURL(options.outputFile))


if __name__ == '__main__':
    main()
