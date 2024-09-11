#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

"""Wrapper to run the cactus_workflow progressively, using the input species tree as a guide

tree.
"""

import logging
import os
import sys
import xml.etree.ElementTree as ET
import timeit
from argparse import ArgumentParser
from base64 import b64encode

from toil.lib.bioio import getTempFile
from cactus.shared.common import cactus_cpu_count
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from toil.job import Job
from toil.common import Toil

from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import findRequiredNode
from cactus.shared.common import makeURL
from cactus.shared.common import catFiles
from cactus.shared.common import cactus_call
from cactus.shared.version import cactus_commit
from cactus.shared.common import cactusRootPath
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import write_s3
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import clean_jobstore_files
from cactus.pipeline.cactus_workflow import cactus_cons_with_resources
from cactus.progressive.progressive_decomposition import compute_outgroups, parse_seqfile, get_subtree, get_spanning_subtree, get_event_set
from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor
from cactus.preprocessor.dnabrnnMasking import loadDnaBrnnModel
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from cactus.paf.local_alignment import make_paf_alignments, trim_unaligned_sequences
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.progressive.cactus_prepare import human2bytesN

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory

def logAssemblyStats(job, message, name, sequenceID, preemptable=True):
    sequenceFile = job.fileStore.readGlobalFile(sequenceID)
    analysisString = cactus_call(parameters=["cactus_analyseAssembly", sequenceFile], check_output=True)
    job.fileStore.logToMaster("%s, got assembly stats for genome %s: %s" % (message, name, analysisString))

def preprocess_all(job, options, config_node, input_seq_id_map):
    ''' run prepreprocessor on every input sequence '''
    root_job = Job()
    job.addChild(root_job)
    events = list(input_seq_id_map.keys())
    seq_ids = list(input_seq_id_map.values())
    preprocessor_job = root_job.addChild(CactusPreprocessor(seq_ids, config_node, eventNames=events))
    pp_seq_ids = {}
    for i, event in enumerate(events):
        pp_seq_ids[event] = preprocessor_job.rv(i)

    # do the logging and checkpointing
    root_job.addFollowOnJobFn(save_preprocessed_files, options, config_node, pp_seq_ids)
    
    return pp_seq_ids

def save_preprocessed_files(job, options, config_node, seq_id_map):
    ''' checkpoint the the intermediate url; also log the stats '''
    # Save preprocessed sequences
    if options.intermediateResultsUrl is not None:
        for genome, seqID in list(seq_id_map.items()):
            job.fileStore.exportFile(seqID, options.intermediateResultsUrl + '-preprocessed-' + genome)

    # Log the stats for the preprocessed assemblies
    for name, sequence in list(seq_id_map.items()):
        job.addChildJobFn(logAssemblyStats, "After preprocessing", name, sequence)

def progressive_schedule(job, options, config_node, seq_id_map, tree, og_map, root_event):
    ''' create job for every internal node, use child dependencies to make tree (+ outgroups)'''

    root_job = Job()
    job.addChild(root_job)

    config_wrapper = ConfigWrapper(config_node)

    # event -> dependencies
    # dependencies are its children + outgroups (which is what get_subtree() returns)
    dep_table = {}
    subtree_event_set = set(tree.getSubtreeRootNames())
    
    iteration_0_events = []
    for node in tree.postOrderTraversal(tree.getNodeId(root_event)):
        event = tree.getName(node)
        if event in subtree_event_set:
            dep_table[event] = []
            subtree = get_subtree(tree, event, config_wrapper, og_map)
            dep_set = set()
            input_count = 0
            for leaf_id in subtree.getLeaves():
                leaf_name = subtree.getName(leaf_id)
                dep_table[event].append(leaf_name)
                if leaf_name in seq_id_map:
                    input_count += 1
            if input_count == len(dep_table[event]):
                iteration_0_events.append(event)
    assert len(dep_table) > 0

    # make the jobs.  jobs must be made after their dependencies.  we just brute force it
    # todo: can speed up with better indexing/updating
    job_table = {}
    fasta_results = seq_id_map
    hal_results = {}
    while len(job_table) != len(dep_table):
        num_jobs_at_iteration_start = len(job_table)
        events = [k for k in dep_table.keys()]
        for event in events:
            if event in job_table:
                continue
            unscheduled_deps = False
            for dep in dep_table[event]:
                if dep not in fasta_results:
                    unscheduled_deps = True
                    break
            if not unscheduled_deps:
                # found an event that a) hasn't been made into a job yet an b) has an existing job for every dep
                event_id_map = {}
                for dep in dep_table[event]:
                    event_id_map[dep] = fasta_results[dep]
                # to be consistent with pre-refactor (and work with updating tests), we include the root when its id is input
                if event in seq_id_map and seq_id_map[event]:
                    event_id_map[event] = seq_id_map[event]                    
                event_job = Job.wrapJobFn(progressive_step, options, config_node, event_id_map, tree, og_map, event)
                job_table[event] = event_job
                for dep in dep_table[event]:
                    if dep in job_table:
                        # add follow-on edge between this job and all jobs that read its results
                        job_table[dep].addFollowOn(event_job)
                hal_results[event] = (event_job.rv(1), event_job.rv(2))
                fasta_results[event] = event_job.rv(3)
        if len(job_table) == num_jobs_at_iteration_start and len(job_table) != len(dep_table):
            raise RuntimeError("Unable to schedule job dependencies.  Please file a bug report on github")

    # hook up or starting jobs (dependent only on input)
    for event in iteration_0_events:
        root_job.addChild(job_table[event])
            
    return hal_results
                
def progressive_step(job, options, config_node, seq_id_map, tree, og_map, event):
    ''' run the blast -> consolidated workflow on an event'''

    # get our subtree (just ingroups and outgroups)
    subtree = get_subtree(tree, event, ConfigWrapper(config_node), og_map)

    # get our events
    subtree_eventmap = {}
    for leaf in subtree.getLeaves():
        subtree_eventmap[subtree.getName(leaf)] = seq_id_map[subtree.getName(leaf)]

    # to be consistent with pre-refactor (and work with updating tests), we include the root when its id is input
    if event in seq_id_map and seq_id_map[event]:
        subtree_eventmap[event] = seq_id_map[event]

    # get the spanning tree (which is what consolidated wants)
    spanning_tree = get_spanning_subtree(tree, event, ConfigWrapper(config_node), og_map)

    # do the blast
    paf_job = job.addChildJobFn(make_paf_alignments, NXNewick().writeString(spanning_tree),
                                subtree_eventmap, event, config_node).encapsulate()

    outgroups = og_map[event] if event in og_map else []
    # trim the outgroups
    if outgroups and int(config_node.find("blast").attrib["trimOutgroups"]):  # Trim the outgroup sequences
        trim_sequences = paf_job.addChildJobFn(trim_unaligned_sequences,
                                               [subtree_eventmap[i] for i in outgroups], paf_job.rv(), config_node)
        cons_job = paf_job.addFollowOnJobFn(progressive_step_2, trim_sequences.rv(), options, config_node, subtree_eventmap,
                                            spanning_tree, og_map, event)
        
    else:  # Without outgroup trimming (or if there are no outgroups to trim)
        cons_job = paf_job.addChildJobFn(cactus_cons_with_resources, spanning_tree, event, config_node, subtree_eventmap,
                                         og_map, paf_job.rv(), cons_cores=options.consCores, cons_memory=options.consMemory,
                                         intermediate_results_url=options.intermediateResultsUrl)
    # erase the paf since its now longer needed
    cons_job.addFollowOnJobFn(clean_jobstore_files, file_ids=[paf_job.rv()])
    return cons_job.rv()


def progressive_step_2(job, trimmed_outgroups_and_alignments, options, config_node, subtree_eventmap,
                       spanning_tree, og_map, event):
    trimmed_outgroup_seqs, pafs = trimmed_outgroups_and_alignments  # unpack the pafs and outgroup sequences
    # set the outgroup seqs
    for outgroup, sequence in zip(og_map[event] if event in og_map else [], trimmed_outgroup_seqs):
        subtree_eventmap[outgroup] = sequence

    # now do consolidated
    return job.addChildJobFn(cactus_cons_with_resources, spanning_tree, event, config_node, subtree_eventmap, og_map,
                             pafs, cons_cores=options.consCores, cons_memory=options.consMemory,
                             intermediate_results_url=options.intermediateResultsUrl).rv()


def export_hal(job, mc_tree, config_node, seq_id_map, og_map, results, event=None, cacheBytes=None,
               cacheMDC=None, cacheRDC=None, cacheW0=None, chunk=None, deflate=None, inMemory=False,
               checkpointInfo=None, acyclicEvent=None, has_resources=False, memory_override=None):

    # todo: going through list nonsense because (i think) it helps with promises, should at least clean up
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, '{}.hal'.format(event if event else mc_tree.getRootName()))

    assert isinstance(mc_tree, MultiCactusTree)
    subtree_roots = set(mc_tree.getSubtreeRootNames())

    # find subtree if event specified
    root_node = None
    if event is not None:
        assert event in mc_tree.nameToId and not mc_tree.isLeaf(mc_tree.nameToId[event])
        root_node = mc_tree.nameToId[event]

    hal_path = os.path.join(work_dir, '{}.hal'.format(event if event else mc_tree.getRootName()))

    if not has_resources:
        fa_file_ids = []
        c2h_file_ids = []
    
    for node in mc_tree.breadthFirstTraversal(root_node):
        genome_name = mc_tree.getName(node)              
        if genome_name in subtree_roots:
            outgroups = og_map[genome_name] if genome_name in og_map else []
            subtree = get_subtree(mc_tree, genome_name, ConfigWrapper(config_node), og_map, include_outgroups=False)
            tree_path = os.path.join(work_dir, '{}.newick'.format(genome_name))
            NXNewick().writeFile(tree_path, nxTree=subtree)
            sub_hal_path = os.path.join(work_dir, '{}.hal.c2h'.format(genome_name))
            hal_fasta_path = os.path.join(work_dir, '{}.hal.fa'.format(genome_name))
            assert genome_name in results
            if not has_resources:
                c2h_file_ids.append(results[genome_name][0])
                fa_file_ids.append(results[genome_name][1])                
            else:
                job.fileStore.readGlobalFile(results[genome_name][0], sub_hal_path)
                job.fileStore.readGlobalFile(results[genome_name][1], hal_fasta_path)
                
                args = ['halAppendCactusSubtree', sub_hal_path, hal_fasta_path, tree_path, hal_path]

                if len(outgroups) > 0:
                    args += ["--outgroups", ",".join(outgroups)]
                if cacheBytes is not None:
                    args += ["--cacheBytes", cacheBytes]
                if cacheMDC is not None:
                    args += ["--cacheMDC", cacheMDC]
                if cacheRDC is not None:
                    args += ["--cacheRDC", cacheRDC]
                if cacheW0 is not None:
                    args += ["--cacheW0", cacheW0]
                if chunk is not None:
                    args += ["--chunk", chunk]
                if deflate is not None:
                    args += ["--deflate", deflate]
                if inMemory is True:
                    args += ["--inMemory"]

                cactus_call(parameters=args, job_memory=job.memory)

    if not has_resources:
        disk = 3 * sum([file_id.size for file_id in fa_file_ids + c2h_file_ids])
        mem = cactus_clamp_memory(5 * (max([file_id.size for file_id in fa_file_ids]) + max([file_id.size for file_id in c2h_file_ids])))
        # allows pass-through of memory override from --consMemory
        if memory_override:
            mem = memory_override
        return job.addChildJobFn(export_hal, mc_tree, config_node, seq_id_map, og_map, results, event=event,
                                 cacheBytes=cacheBytes, cacheMDC=cacheMDC, cacheRDC=cacheRDC, cacheW0=cacheW0,
                                 chunk=chunk, deflate=deflate, inMemory=inMemory, checkpointInfo=checkpointInfo,
                                 acyclicEvent=acyclicEvent, has_resources=True,
                                 disk=disk, memory=mem).rv()

    cactus_call(parameters=["halSetMetadata", hal_path, "CACTUS_COMMIT", cactus_commit])
    config_path = os.path.join(work_dir, 'config.xml')
    ConfigWrapper(config_node).writeXML(config_path)
    with open(config_path, 'rb') as configFile:
        cactus_call(parameters=["halSetMetadata", hal_path, "CACTUS_CONFIG", b64encode(configFile.read()).decode()])

    if acyclicEvent:
        cactus_call(parameters=["halRemoveDupes", hal_path, acyclicEvent], job_memory=job.memory)

    if checkpointInfo:
        write_s3(hal_path, checkpointInfo[1], region=checkpointInfo[0])

    return job.fileStore.writeGlobalFile(hal_path)
    
def progressive_workflow(job, options, config_node, mc_tree, og_map, input_seq_id_map):
    ''' run the entire progressive workflow '''

    # run the usual unzip / rename, even before preprocessing
    sanitize_job = job.addChildJobFn(sanitize_fasta_headers, input_seq_id_map)
    
    # start with the preprocessor
    if not options.skipPreprocessor:
        pp_job = sanitize_job.addFollowOnJobFn(preprocess_all, options, config_node, sanitize_job.rv())
        seq_id_map = pp_job.rv()
        sanitize_job = pp_job
    else:
        seq_id_map = sanitize_job.rv()

    # then do the progressive workflow
    root_event = options.root if options.root else mc_tree.getRootName()
    progressive_job = sanitize_job.addFollowOnJobFn(progressive_schedule, options, config_node, seq_id_map, mc_tree, og_map, root_event)

    # then do the hal export
    hal_export_job = progressive_job.addFollowOnJobFn(export_hal, mc_tree, config_node, seq_id_map, og_map,
                                                      progressive_job.rv(), event=root_event, memory_override=options.consMemory)

    return hal_export_job.rv()

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("outputHal", type=str, help = "Output HAL file")

    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--root", dest="root", help="Name of ancestral node (which"
                      " must appear in NEWICK tree in <seqfile>) to use as a "
                      "root for the alignment.  Any genomes not below this node "
                      "in the tree may be used as outgroups but will never appear"
                      " in the output.  If no root is specifed then the root"
                      " of the tree is used. ", default=None)
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
    parser.add_argument("--consCores", type=int, 
                        help="Number of cores for each cactus_consolidated job (defaults to all available / maxCores on single_machine)", default=None)
    parser.add_argument("--consMemory", type=human2bytesN,
                        help="Memory in bytes for each cactus_consolidated job (defaults to an estimate based on the input data size). "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))", default=None)
    parser.add_argument("--intermediateResultsUrl",
                        help="URL prefix to save intermediate results like DB dumps to (e.g. "
                        "prefix-dump-caf, prefix-dump-avg, etc.)", default=None)
    parser.add_argument("--skipPreprocessor",
                        help="Do not run any preprocessor jobs",
                        action="store_true")
    parser.add_argument("--maxOutgroups",
                        help="Number of outgroups to use for each alignment. Overrides `max_num_outgroups` from the config XML."
                        " Using many outgroups can significantly increase alignment time, especially around long branchecs.",
                        type=int, default=None)
    parser.add_argument("--chromInfo",
                        help="Two-column file mapping genome (col 1) to comma-separated list of sex chromosomes. This information "
                        "will be used to guide outgroup selection so that, where possible, all chromosomes are present in"
                        " at least one outgroup.")

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    # Todo: do we even need --consCores anymore?
    if options.maxCores is not None and options.maxCores < 2**20 and options.consCores is None:
        options.consCores = int(options.maxCores)

    # Try to juggle --maxCores and --consCores to give some reasonable defaults where possible
    if options.batchSystem.lower() in ['single_machine', 'singlemachine']:
        if options.maxCores is not None and options.maxCores < 2**20:
            if int(options.maxCores) <= 0:
                raise RuntimeError('Cactus requires --maxCores >= 1')
        if options.consCores is None:
            if options.maxCores is not None and options.maxCores < 2**20:
                options.consCores = int(options.maxCores)
            else:
                options.consCores = cactus_cpu_count()
        elif options.maxCores is not None and options.maxCores < 2**20 and options.consCores > int(options.maxCores):
            raise RuntimeError('--consCores must be <= --maxCores')
    else:
        if not options.consCores:
            raise RuntimeError('--consCores required for non single_machine batch systems')
    if options.maxCores is not None and options.maxCores < 2**20 and options.consCores > int(options.maxCores):
        raise RuntimeError('--consCores must be <= --maxCores')
    
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    start_time = timeit.default_timer()
    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))

    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            hal_id = toil.restart()
        else:

            # load up the seqfile and figure out the outgroups and schedule
            config_node = ET.parse(options.configFile).getroot()
            config_wrapper = ConfigWrapper(config_node)
            config_wrapper.substituteAllPredefinedConstantsWithLiterals(options)
            config_wrapper.setSystemMemory(options)
            config_wrapper.applySlurmChunkScaling(options)
            if options.maxOutgroups:
                config_wrapper.setMaxNumOutgroups(options.maxOutgroups)

            # apply gpu override
            config_wrapper.initLastz(options)
            mc_tree, input_seq_map, og_candidates = parse_seqfile(options.seqFile, config_wrapper, root_name = options.root)
            logger.info('Tree: {}'.format(NXNewick().writeString(mc_tree)))
            og_map = compute_outgroups(mc_tree, config_wrapper, set(og_candidates), options.root,
                                       chrom_info_file = options.chromInfo)
            og_log = ''
            for source, ogs_dists in og_map.items():
                og_log += ' {}: {}\n'.format(source, ogs_dists)
            logger.info('Outgroups:\n{}'.format(og_log))
            event_set = get_event_set(mc_tree, config_wrapper, og_map, options.root, subtree=False)
            # infer default root
            if not options.root:
                options.root = mc_tree.getRootName()
            # take presence of root as sign we want to include it
            if options.root in input_seq_map and input_seq_map[options.root]:
                event_set.add(options.root)
                        
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
                
            # Make sure we have the dna-brnn model in the filestore if we need it
            loadDnaBrnnModel(toil, config_node)

            # run the whole workflow
            hal_id = toil.start(Job.wrapJobFn(progressive_workflow, options, config_node, mc_tree, og_map, input_seq_id_map))

        toil.exportFile(hal_id, makeURL(options.outputHal))
    
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("Cactus has finished after {} seconds".format(run_time))

if __name__ == '__main__':
    main()
