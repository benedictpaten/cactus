#!/usr/bin/env python3

"""
Generate the alignment file that needs to be input to cactus_consolidated

Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

from toil.job import Job
from toil.statsAndLogging import logger
from toil.lib.bioio import getLogLevelString
from sonLib.bioio import newickTreeParser
import os
import shutil
import math
from cactus.paf.paf import get_event_pairs, get_leaves, get_node, get_distances
from cactus.shared.common import cactus_call, getOptionalAttrib
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from cactus.shared.common import cactus_clamp_memory

def run_lastz(job, name_A, genome_A, name_B, genome_B, distance, params):
    # Create a local temporary file to put the alignments in.
    work_dir = job.fileStore.getLocalTempDir()
    alignment_file = os.path.join(work_dir, '{}_{}.paf'.format(name_A, name_B))

    # Get the input files
    genome_a_file = os.path.join(work_dir, '{}.fa'.format(name_A))
    genome_b_file = os.path.join(work_dir, '{}.fa'.format(name_B))
    job.fileStore.readGlobalFile(genome_A, genome_a_file)
    job.fileStore.readGlobalFile(genome_B, genome_b_file)

    # Get the params to do the alignment
    lastz_params_node = params.find("blast")
    lastz_divergence_node = lastz_params_node.find("lastzArguments")
    divergences = params.find("constants").find("divergences")
    for i in "one", "two", "three", "four", "five":
        if distance <= float(divergences.attrib[i]):
            lastz_params = lastz_divergence_node.attrib[i]
            break
    else:
        lastz_params = lastz_divergence_node.attrib["default"]
    logger.info("For distance {} for genomes {}, {} using {} lastz parameters".format(distance, genome_A,
                                                                                      genome_B, lastz_params))

    gpu = getOptionalAttrib(lastz_params_node, 'gpu', typeFn=int, default=0)
    cpu = getOptionalAttrib(lastz_params_node, 'cpu', typeFn=int, default=None)    
    if gpu:
        lastz_bin = 'run_segalign'
        suffix_a, suffix_b = '', ''
        assert gpu > 0
        lastz_params += ' --num_gpu {}'.format(gpu)
    else:
        lastz_bin = 'lastz'
        suffix_a = '[multiple][nameparse=darkspace]'
        suffix_b = '[nameparse=darkspace]'

    # Generate the alignment
    lastz_cmd = [lastz_bin,
                 '{}{}'.format(os.path.basename(genome_a_file), suffix_a),
                 '{}{}'.format(os.path.basename(genome_b_file), suffix_b),
                 '--format=paf:minimap2'] + lastz_params.split(' ')
        
    # note: it's very important to set the work_dir here, because cactus_call is not able to
    # sort out the mount directory by itself, presumably due to square brackets...
    segalign_messages = cactus_call(parameters=lastz_cmd, outfile=alignment_file, work_dir=work_dir, returnStdErr=gpu>0, gpus=gpu,
                                    cpus=cpu, job_memory=job.memory)

    if gpu:
        # run_segalign can crash and still exit 0, so it's worth taking a moment to check the log for errors
        segalign_messages = segalign_messages.lower()
        for line in segalign_messages.split("\n"):
            if not line.startswith("signals delivered"):
                for keyword in ['terminate', 'error', 'fail', 'assert', 'signal', 'abort', 'segmentation', 'sigsegv', 'kill']:
                    if keyword in line and 'signals' not in line:
                        job.fileStore.logToMaster("Segalign offending line: " + line)  # Log the messages
                        raise RuntimeError('{} exited 0 but keyword "{}" found in stderr'.format(lastz_cmd, keyword))

    # Add in mismatches so we can trim alignments later
    alignment_file_with_mismatches = os.path.join(work_dir, '{}_{}_mismatches.paf'.format(name_A, name_B))
    cactus_call(parameters=['paffy', 'add_mismatches', "-i", alignment_file, genome_a_file, genome_b_file ],
                outfile=alignment_file_with_mismatches)

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file_with_mismatches)


def run_minimap2(job, name_A, genome_A, name_B, genome_B, distance, params):
    # Create a local temporary file to put the alignments in.
    work_dir = job.fileStore.getLocalTempDir()
    alignment_file = os.path.join(work_dir, '{}_{}.paf'.format(name_A, name_B))

    # Get the input files
    genome_a_file = os.path.join(work_dir, '{}.fa'.format(name_A))
    genome_b_file = os.path.join(work_dir, '{}.fa'.format(name_B))
    job.fileStore.readGlobalFile(genome_A, genome_a_file)
    job.fileStore.readGlobalFile(genome_B, genome_b_file)

    minimap2_params = params.find("blast").attrib["minimap2_params"]
    
    # Generate the alignment
    minimap2_cmd = ['minimap2', '-c', genome_a_file, genome_b_file] + minimap2_params.split(' ')
    cactus_call(parameters=minimap2_cmd, outfile=alignment_file, job_memory=job.memory)

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file)

def run_fastga(job, name_A, genome_A, name_B, genome_B, distance, params):
    # Create a local temporary file to put the alignments in.
    work_dir = job.fileStore.getLocalTempDir()
    fastga_tempdir = os.path.join(work_dir, 'fgatmp')
    os.mkdir(fastga_tempdir)
    
    alignment_file = os.path.join(work_dir, '{}_{}.paf'.format(name_A, name_B))

    # Get the input files
    genome_a_file = os.path.join(work_dir, '{}.fa'.format(name_A))
    genome_b_file = os.path.join(work_dir, '{}.fa'.format(name_B))
    job.fileStore.readGlobalFile(genome_A, genome_a_file)
    job.fileStore.readGlobalFile(genome_B, genome_b_file)

    fastga_params = params.find("blast").attrib["fastga_params"]
    
    # Generate the alignment
    # note: FastGA very picky about spacking and ordering of parameters
    fastga_cmd = ['FastGA']
    if fastga_params:
        fastga_cmd += fastga_params.split(' ')
    fastga_cmd += ['-P{}'.format(fastga_tempdir), '-T{}'.format(job.cores), genome_b_file, genome_a_file]
    cactus_call(parameters=fastga_cmd, outfile=alignment_file, job_memory=job.memory)

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file)


def combine_chunks(job, chunked_alignment_files, batch_size):
    if len(chunked_alignment_files) >= 2 * batch_size:
        # run combine_chunks in batches
        root_job = Job()
        job.addChild(root_job)
        batch_results = []
        for chunk_idx in range(math.ceil(len(chunked_alignment_files) / batch_size)):
            batch = chunked_alignment_files[chunk_idx * batch_size : chunk_idx * batch_size + batch_size]
            batch_results.append(root_job.addChildJobFn(combine_chunks, batch, batch_size,
                                                        disk=sum([f.size for f in batch])).rv())
        return root_job.addFollowOnJobFn(merge_combined_chunks, batch_results,
                                         disk=2*sum([f.size for f in chunked_alignment_files])).rv()
    else:
        # Make combined alignments file
        alignment_file = job.fileStore.getLocalTempFile()
        for i, chunk in enumerate(chunked_alignment_files): # Append each of the chunked alignment files into one file
            cactus_call(parameters=['paffy', 'dechunk', '-i', job.fileStore.readGlobalFile(chunk)],
                        outfile=alignment_file, outappend=True, rt_log_cmd=(i<5 or i>=len(chunked_alignment_files)-5))
            job.fileStore.deleteGlobalFile(chunk) # Cleanup the old files
        return job.fileStore.writeGlobalFile(alignment_file)  # Return the final alignments file copied to the jobstore


def merge_combined_chunks(job, combined_chunks):
    # mege up and return some chunks, deleting them too
    output_path = job.fileStore.getLocalTempFile()
    with open(output_path, 'a') as output_file:
        for i,chunk in enumerate(combined_chunks):
            chunk_path = job.fileStore.readGlobalFile(chunk, mutable=True)
            with open(chunk_path, 'r') as chunk_file:
                shutil.copyfileobj(chunk_file, output_file)
            job.fileStore.deleteGlobalFile(chunk)
    return job.fileStore.writeGlobalFile(output_path)


def make_chunked_alignments(job, event_a, genome_a, event_b, genome_b, distance, params):
    lastz_params_node = params.find("blast")
    gpu = getOptionalAttrib(lastz_params_node, 'gpu', typeFn=int, default=0)
    fastga = getOptionalAttrib(lastz_params_node, 'mapper', typeFn=str) == 'fastga'
    if gpu or fastga:
        # wga-gpu has a 6G limit, so we always override
        lastz_params_node.attrib['chunkSize'] = '6000000000'
    lastz_cores = getOptionalAttrib(lastz_params_node, 'cpu', typeFn=int, default=None)
    lastz_memory = getOptionalAttrib(lastz_params_node, 'lastz_memory', typeFn=int, default=None)
    
    def make_chunks(genome):
        output_chunks_dir = job.fileStore.getLocalTempDir()
        fasta_chunk_cmd = ['faffy', 'chunk',
                           '-c', params.find("blast").attrib["chunkSize"],
                           '-o', params.find("blast").attrib["overlapSize"],
                           '--dir', output_chunks_dir,
                           job.fileStore.readGlobalFile(genome)]
        cactus_call(parameters=fasta_chunk_cmd)
        return [job.fileStore.writeGlobalFile(os.path.join(output_chunks_dir, chunk), cleanup=True)
                for chunk in os.listdir(output_chunks_dir)]
    # Chunk each input genome
    chunks_a = make_chunks(genome_a)
    chunks_b = make_chunks(genome_b)

    # Align all chunks from genome_A against all chunks from genome_B
    accelerators = 'cuda:{}'.format(gpu) if gpu else None
    chunked_alignment_files = []
    for i, chunk_a in enumerate(chunks_a):
        for j, chunk_b in enumerate(chunks_b):
            mappers = { "lastz":run_lastz, "minimap2":run_minimap2, "fastga":run_fastga}
            mappingFn = mappers[params.find("blast").attrib["mapper"]]
            memory = lastz_memory if lastz_memory else max(200000000, 15*(chunk_a.size+chunk_b.size))
            chunked_alignment_files.append(job.addChildJobFn(mappingFn, '{}_{}'.format(event_a, i), chunk_a,
                                                             '{}_{}'.format(event_b, j), chunk_b, distance, params,
                                                             cores=lastz_cores,
                                                             disk=max(4*(chunk_a.size+chunk_b.size), memory),
                                                             memory=cactus_clamp_memory(memory),
                                                             accelerators=accelerators).rv())

    dechunk_batch_size = getOptionalAttrib(lastz_params_node, 'dechunkBatchSize', typeFn=int, default=1e9)
    return job.addFollowOnJobFn(combine_chunks, chunked_alignment_files, dechunk_batch_size).rv()  # Combine the chunked alignment files


def invert_alignments(job, alignment_file):
    """ Invert the pafs in the alignment_file """
    alignment_file_local = job.fileStore.readGlobalFile(alignment_file)
    inverted_alignment_file_local = job.fileStore.getLocalTempFile()  # Get a temporary file to store the alignments in
    cactus_call(parameters=['paffy', 'invert', "-i", alignment_file_local], outfile=inverted_alignment_file_local, outappend=True,
                job_memory=job.memory)
    job.fileStore.deleteGlobalFile(alignment_file)
    return job.fileStore.writeGlobalFile(inverted_alignment_file_local)


def make_ingroup_to_outgroup_alignments_0(job, ingroup_event, outgroup_events, event_names_to_sequences, distances, params):
    # Generate the alignments fle
    alignment_file = job.addChildJobFn(make_ingroup_to_outgroup_alignments_1, ingroup_event, outgroup_events,
                                            event_names_to_sequences, distances, params).rv()

    # Invert the final alignment so that the query is the outgroup and the target is the ingroup
    return job.addFollowOnJobFn(invert_alignments, alignment_file).rv()


def make_ingroup_to_outgroup_alignments_1(job, ingroup_event, outgroup_events, event_names_to_sequences, distances, params):
    #  a job should never set its own follow-on, so we hang everything off root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)

    #  align ingroup to first outgroup to produce paf alignments
    outgroup = outgroup_events[0] # The first outgroup
    logger.info("Building alignment between ingroup event: {} and outgroup event: {}".format(ingroup_event.iD, outgroup.iD))
    alignment = root_job.addChildJobFn(make_chunked_alignments,
                                       outgroup.iD, event_names_to_sequences[outgroup.iD],
                                       ingroup_event.iD, event_names_to_sequences[ingroup_event.iD], distances[ingroup_event, outgroup], params,
                                       disk=4*(event_names_to_sequences[ingroup_event.iD].size+event_names_to_sequences[outgroup.iD].size)).rv()

    #  post process the alignments and recursively generate alignments to remaining outgroups
    return root_job.addFollowOnJobFn(make_ingroup_to_outgroup_alignments_2, alignment, ingroup_event, outgroup_events[1:],
                                     event_names_to_sequences, distances, params,
                                     disk=8*(event_names_to_sequences[ingroup_event.iD].size+event_names_to_sequences[outgroup.iD].size),
                                     memory=cactus_clamp_memory(8*(event_names_to_sequences[ingroup_event.iD].size+event_names_to_sequences[outgroup.iD].size))).rv() if len(outgroup_events) > 1 else alignment


def make_ingroup_to_outgroup_alignments_2(job, alignments, ingroup_event, outgroup_events,
                                          event_names_to_sequences, distances, params):
    
    # a job should never set its own follow-on, so we hang everything off root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)
        
    # identify all ingroup sub-sequences that remain unaligned longer than a threshold as follows:

    # run paffy to_bed to create a bed of aligned coverage
    work_dir = job.fileStore.getLocalTempDir()
    alignments_file = os.path.join(work_dir, '{}.paf'.format(ingroup_event.iD))
    job.fileStore.readGlobalFile(alignments, alignments_file)
    bed_file = os.path.join(work_dir, '{}.bed'.format(ingroup_event.iD))  # Get a temporary file to store the bed file in
    messages = cactus_call(parameters=['paffy', 'to_bed', "--excludeAligned", "--binary",
                                       "--minSize", params.find("blast").attrib["trimMinSize"],
                                       "-i", alignments_file,
                                       "--logLevel", getLogLevelString()],
                                       outfile=bed_file, returnStdErr=True, job_memory=job.memory)
    logger.info("paffy to_bed event:{}\n{}".format(ingroup_event.iD, messages[:-1]))  # Log paffy to_bed

    # run faffy extract to extract unaligned sequences longer than a threshold creating a reduced subset of A
    seq_file = os.path.join(work_dir, '{}_subseq.fa'.format(ingroup_event.iD))  # Get a temporary file to store the subsequences in
    ingroup_seq_file = os.path.join(work_dir, '{}.fa'.format(ingroup_event.iD))
    job.fileStore.readGlobalFile(event_names_to_sequences[ingroup_event.iD], ingroup_seq_file)  # The ingroup sequences
    messages = cactus_call(parameters=['faffy', 'extract', "-i", bed_file, ingroup_seq_file,
                                       "--flank", params.find("blast").attrib["trimFlanking"],
                                       "--logLevel", getLogLevelString()],
                           outfile=seq_file, returnStdErr=True, job_memory=job.memory)
    logger.info("faffy extract event:{}\n{}".format(ingroup_event.iD, messages[:-1]))  # Log faffy extract

    # replace the ingroup sequences with remaining sequences
    event_names_to_sequences[ingroup_event.iD] = job.fileStore.writeGlobalFile(seq_file)

    # recursively make alignments with the remaining outgroups
    alignments2 = root_job.addChildJobFn(make_ingroup_to_outgroup_alignments_1, ingroup_event, outgroup_events,
                                         event_names_to_sequences, distances, params).rv()

    return root_job.addFollowOnJobFn(make_ingroup_to_outgroup_alignments_3, ingroup_event, event_names_to_sequences[ingroup_event.iD],
                                     alignments, alignments2).rv()


def make_ingroup_to_outgroup_alignments_3(job, ingroup_event, ingroup_seq_file, alignments, alignments2, has_resources=False):
    # Merge the together two alignment files

    if not has_resources:
        # unpack promises for disk requirement
        return job.addChildJobFn(make_ingroup_to_outgroup_alignments_3, ingroup_event, ingroup_seq_file, alignments,
                                 alignments2, has_resources=True, disk=3*(alignments.size + alignments2.size)).rv()

    alignments = job.fileStore.readGlobalFile(alignments)  # Copy the global alignment files locally
    alignments2 = job.fileStore.readGlobalFile(alignments2)

    # use paffy dechunk to correct the subsequence coordinates of alignments2
    alignments2_corrected = job.fileStore.getLocalTempFile()
    messages = cactus_call(parameters=['paffy', 'dechunk', "-i", alignments2, "--query", "--logLevel", getLogLevelString()],
                           outfile=alignments2_corrected, returnStdErr=True, job_memory=job.memory)
    logger.info("paffy dechunk event:{}\n{}".format(ingroup_event.iD, messages[:-1]))  # Log paffy dechunk

    # merge the two alignment files together
    merged_alignments = job.fileStore.getLocalTempFile()
    cactus_call(parameters=['cat', alignments, alignments2_corrected], outfile=merged_alignments, job_memory=job.memory)  # Don't bother to log

    # Delete the file containing the subset of ingroup sequences as we don't need it any longer and it takes up space
    job.fileStore.deleteGlobalFile(ingroup_seq_file)

    return job.fileStore.writeGlobalFile(merged_alignments)


def merge_alignments(job, alignment_file1, alignment_file2):
    """" Merge together two alignment files """
    # Get a temporary directory to work in
    work_dir = job.fileStore.getLocalTempDir()

    # The merged alignment file
    merged_alignments_file = os.path.join(work_dir, 'output_alignments.paf')

    # Read the files to local disk and merge them together, don't bother to log
    cactus_call(parameters=['cat',
                            job.fileStore.readGlobalFile(alignment_file1, os.path.join(work_dir, 'alignment_1.paf')),
                            job.fileStore.readGlobalFile(alignment_file2, os.path.join(work_dir, 'alignment_2.paf'))],
                outfile=merged_alignments_file, job_memory=job.memory)

    # Write the merged file back to the job store and return its path
    return job.fileStore.writeGlobalFile(merged_alignments_file)


def chain_alignments_splitting_ingroups_and_outgroups(job, ingroup_alignment_files, ingroup_alignment_names,
                                                      outgroup_alignment_files, outgroup_alignment_names,
                                                      reference_event_name, params):
    """ Chains/tiles/etc. the ingroup alignments to each other so that every ingroup sequence has a primary alignment to
    another ingroup sequence, and separately each outgroup sequence has a primary alignment to an ingroup."""

    # Check we have the expected number of alignment files
    assert len(ingroup_alignment_files) == len(ingroup_alignment_names)
    assert len(outgroup_alignment_files) == len(outgroup_alignment_names)

    # Chain and pick the primary alignments of the ingroups to each other
    chained_ingroup_alignments = job.addChildJobFn(chain_alignments, ingroup_alignment_files,
                                                   ingroup_alignment_names, reference_event_name, params).rv()

    # Separately pick the primary of the outgroups to the ingroups. By setting include_inverted_alignments=False
    # we only get outgroup-to-ingroup alignments and not imgroup-to-ouygroup alignments and therefore primary
    # alignments along the outgroup sequences
    chained_outgroup_alignments = job.addChildJobFn(chain_alignments, outgroup_alignment_files,
                                                    outgroup_alignment_names, reference_event_name, params,
                                                    include_inverted_alignments=False).rv()

    # Calculate approximately total alignment file size
    total_file_size = sum(alignment_file.size for alignment_file in ingroup_alignment_files + outgroup_alignment_files)

    # Merge the resulting two alignment files into a single set of alignments
    return job.addFollowOnJobFn(merge_alignments, chained_ingroup_alignments, chained_outgroup_alignments,
                                disk=total_file_size*10).rv()


def chain_alignments(job, alignment_files, alignment_names, reference_event_name, params,
                     include_inverted_alignments=True):
    """The following recapitulates the pipeline showing in paffy/tests/paf_pipeline_test.sh"""

    root_job = Job()
    job.addChild(root_job)

    assert len(alignment_files) == len(alignment_names)
    
    # do the chaining
    chained_alignment_files = []
    for alignment_file, alignment_name in zip(alignment_files, alignment_names):
        chained_alignment_files.append(root_job.addChildJobFn(chain_one_alignment, alignment_file, alignment_name,
                                                              params, include_inverted_alignments,
                                                              disk=6*alignment_file.size,
                                                              memory=cactus_clamp_memory(32*alignment_file.size)).rv())
        
    # do the tiling and filtering
    return root_job.addFollowOnJobFn(tile_alignments, chained_alignment_files, reference_event_name, params,
                                     disk=6*sum([alignment_file.size for alignment_file in alignment_files]),
                                     memory=cactus_clamp_memory(32*sum([alignment_file.size for alignment_file in alignment_files]))).rv()


def chain_one_alignment(job, alignment_file, alignment_name, params, include_inverted_alignments):
    """ run paffy chain on one PAF. include_inverted_alignemnts is a flag to control if we additionally include
    the inverted paf records for chaining.
    """

    work_dir = job.fileStore.getLocalTempDir()
    alignment_path = os.path.join(work_dir, alignment_name + '.paf')
    alignment_inv_path = os.path.join(work_dir, alignment_name + '.inv.paf')
    output_path = os.path.join(work_dir, alignment_name + '.chained.paf')

    # Copy the alignments from the job store
    job.fileStore.readGlobalFile(alignment_file, alignment_path)

    # Get the forward and reverse versions of each alignment for symmetry with chaining if include_inverted_alignments
    # is set
    if include_inverted_alignments:
        shutil.copyfile(alignment_path, alignment_inv_path)
        cactus_call(parameters=['paffy', 'invert', "-i", alignment_inv_path], outfile=alignment_path, outappend=True,
                    job_memory=job.memory)

    # Now chain the alignments
    cactus_call(parameters=['paffy', 'chain', "-i", alignment_path,
                            "--maxGapLength", params.find("blast").attrib["chainMaxGapLength"],
                            "--chainGapOpen", params.find("blast").attrib["chainGapOpen"],
                            "--chainGapExtend", params.find("blast").attrib["chainGapExtend"],
                            "--trimFraction", params.find("blast").attrib["chainTrimFraction"],
                            "--logLevel", getLogLevelString()],
                outfile=output_path, job_memory=job.memory)

    job.fileStore.deleteGlobalFile(alignment_file)

    return job.fileStore.writeGlobalFile(output_path)
    
    
def tile_alignments(job, alignment_files, reference_event_name, params):
    # do everything post-chaining

    work_dir = job.fileStore.getLocalTempDir()

    # concatenate the input into one big paf
    local_paths = [os.path.join(work_dir, 'chained_{}_{}.paf'.format(reference_event_name, i)) for i in range(len(alignment_files))]
    for local_path, alignment_file in zip(local_paths, alignment_files):
        job.fileStore.readGlobalFile(alignment_file, local_path)
    chained_paf_path = os.path.join(work_dir, 'chained_{}'.format(reference_event_name))
    with open(chained_paf_path, 'a') as chained_paf_file:
        for local_path in local_paths:
            with open(local_path, 'r') as local_file:
                shutil.copyfileobj(local_file, chained_paf_file)
    for alignment_file in alignment_files:
        job.fileStore.deleteGlobalFile(alignment_file)

    # Now tile to select the primary alignments
    tiled_paf_path = os.path.join(work_dir, 'tiled_{}.paf'.format(reference_event_name))
    cactus_call(parameters=['paffy', 'tile', "-i", chained_paf_path, "--logLevel", getLogLevelString()],
                           outfile=tiled_paf_path, job_memory=job.memory)

    os.remove(chained_paf_path)
    trimmed_paf_path = os.path.join(work_dir, 'trim_{}.paf'.format(reference_event_name))

    # Trim the poorly aligned tails off the ends of the alignments after tiling - doing so before tiling
    # can create gaps in alignment chains which allows in spurious chains
    cactus_call(parameters=['paffy', 'trim', "-i", tiled_paf_path,
                            "--trimIdentity", params.find("blast").attrib["pafTrimIdentity"]],
                outfile=trimmed_paf_path, job_memory=job.memory)

    os.remove(tiled_paf_path)
    filter_paf_path = os.path.join(work_dir, 'filter_{}.paf'.format(reference_event_name))
                                    
    # Filter to primary alignments
    cactus_call(parameters=['paffy', 'filter', "-i", trimmed_paf_path, "--maxTileLevel", "1"],
                           outfile=filter_paf_path, job_memory=job.memory)

    os.remove(trimmed_paf_path)
    output_alignments_file = os.path.join(work_dir, 'output_alignments.paf')

    # Do we want to include the secondary alignments
    use_secondary_alignments = int(params.find("blast").attrib["outputSecondaryAlignments"])  # We should really switch to
    # the empty string being false instead of 0

    if use_secondary_alignments:
        # Filter to secondary alignments and put in the final output file
        cactus_call(parameters=['paffy', 'filter', "-i", filter_paf_path, "--maxTileLevel", "1", '-x'],
                    outfile=output_alignments_file, job_memory=job.memory)

    primary_chain_paf_path = os.path.join(work_dir, 'primary_chain_{}.paf'.format(reference_event_name))
    
    # Rechain the "primary" alignments so we can see how good the chains of the primary alignments are
    cactus_call(parameters=['paffy', 'chain', "-i", filter_paf_path,
                            "--maxGapLength", params.find("blast").attrib["chainMaxGapLength"],
                            "--chainGapOpen", params.find("blast").attrib["chainGapOpen"],
                            "--chainGapExtend", params.find("blast").attrib["chainGapExtend"],
                            "--trimFraction", params.find("blast").attrib["chainTrimFraction"],
                            "--logLevel", getLogLevelString()],
                outfile=primary_chain_paf_path, job_memory=job.memory)

    os.remove(filter_paf_path)

    # Filter primary alignments not in good chains
    cactus_call(parameters=['paffy', 'filter', "-i", primary_chain_paf_path,
                            "--minChainScore", params.find("blast").attrib["minPrimaryChainScore"]],
                outfile=output_alignments_file, outappend=True, job_memory=job.memory)

    if use_secondary_alignments:
        # Get the primaries we've filtered and switch them to secondaries in the final output
        prim_filter_cmd = [['paffy', 'filter', "-i", primary_chain_paf_path, "-x",
                            "--minChainScore", params.find("blast").attrib["minPrimaryChainScore"]]]
        prim_filter_cmd += [['sed', 's/tp:A:P/tp:A:S/']]
        # Switch low score primaries to secondaries
        prim_filter_cmd += [['sed', 's/tl:i:1/tl:i:2/']]
        cactus_call(parameters=prim_filter_cmd, outfile=output_alignments_file, outappend=True)

    return job.fileStore.writeGlobalFile(output_alignments_file)  # Copy back
        

def sanitize_then_make_paf_alignments(job, event_tree_string, event_names_to_sequences, ancestor_event_string, params):
    sanitize_job = job.addChildJobFn(sanitize_fasta_headers, event_names_to_sequences)
    paf_job = sanitize_job.addFollowOnJobFn(make_paf_alignments, event_tree_string, sanitize_job.rv(), ancestor_event_string, params)
    return paf_job.rv()


def make_paf_alignments(job, event_tree_string, event_names_to_sequences, ancestor_event_string, params):
    # a job should never set its own follow-on, so we hang everything off the root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)

    logger.info("Parsing species tree: {}".format(event_tree_string))
    event_tree = newickTreeParser(event_tree_string)

    ancestor_event = get_node(event_tree, ancestor_event_string)
    ingroup_events = get_leaves(ancestor_event) # Get the set of ingroup events
    # to be consistent with pre-refactor (and work with updating tests), we include the root when its id is input
    # and just treat it as an ingroup
    if ancestor_event.iD in event_names_to_sequences and event_names_to_sequences[ancestor_event.iD]:
        ingroup_events.append(ancestor_event)
    outgroup_events = [event for event in get_leaves(event_tree) if event not in ingroup_events]  # Set of outgroups
    logger.info("Got ingroup events: {} for ancestor event: {}".format(" ".join([i.iD for i in ingroup_events]),
                                                                       ancestor_event_string))

    # Calculate the total sequence size
    total_sequence_size = sum(event_names_to_sequences[event.iD].size for event in get_leaves(event_tree))

    # for each pair of ingroups make alignments
    ingroup_alignments = []
    # for better logs
    ingroup_alignment_names = []
    for ingroup, ingroup2, distance_a_b in get_event_pairs(ancestor_event, ingroup_events):
        logger.info("Building alignment between event: {} (ingroup) and event: {} (ingroup)".format(ingroup.iD, ingroup2.iD))
        ingroup_alignments.append(root_job.addChildJobFn(make_chunked_alignments,
                                                         ingroup.iD, event_names_to_sequences[ingroup.iD],
                                                         ingroup2.iD, event_names_to_sequences[ingroup2.iD], distance_a_b, params,
                                                         disk=2*total_sequence_size).rv())
        ingroup_alignment_names.append('{}-{}_vs_{}'.format(ancestor_event_string, ingroup.iD, ingroup2.iD))

    # Get the outgroup events
    distances = get_distances(event_tree)  # Distances between all pairs of nodes
    outgroup_events.sort(key=lambda outgroup: distances[ancestor_event, outgroup])  # Sort from closest to furthest
    logger.info("Got outgroup events: {} for ancestor event: {}".format(" ".join([i.iD for i in outgroup_events]),
                                                                        ancestor_event.iD))

    # for each ingroup make alignments to the outgroups
    if int(params.find("blast").attrib["trimIngroups"]):  # Trim the ingroup sequences
        outgroup_alignments = [root_job.addChildJobFn(make_ingroup_to_outgroup_alignments_0, ingroup, outgroup_events,
                                                      dict(event_names_to_sequences), distances, params).rv()
                                for ingroup in ingroup_events] if len(outgroup_events) > 0 else []
    else:
        outgroup_alignments = [root_job.addChildJobFn(make_chunked_alignments,
                                                      # Ingroup will be the target, outgroup the query
                                                      ingroup.iD, event_names_to_sequences[ingroup.iD],
                                                      outgroup.iD, event_names_to_sequences[outgroup.iD],
                                                      distances[ingroup, outgroup], params,
                                                      disk=2*total_sequence_size).rv()
                               for ingroup in ingroup_events for outgroup in outgroup_events]
    # for better logs
    outgroup_alignment_names = ['{}-og_{}'.format(ancestor_event_string, i) for i in range(len(outgroup_alignments))]

    # Now do the chaining, either getting ingroup primary alignments separately to primary ingroup to outgroup
    # alignments if the following option is true, otherwise get every sequence to pick its primary alignment without
    # regard to whether the other sequence is an ingroup or an outgroup
    if int(params.find("blast").attrib["pickIngroupPrimaryAlignmentsSeparatelyToOutgroups"]):
        return root_job.addFollowOnJobFn(chain_alignments_splitting_ingroups_and_outgroups,
                                         ingroup_alignments, ingroup_alignment_names,
                                         outgroup_alignments, outgroup_alignment_names,
                                         ancestor_event_string, params).rv()

    return job.addFollowOnJobFn(chain_alignments, ingroup_alignments + outgroup_alignments,
                                ingroup_alignment_names + outgroup_alignment_names,
                                ancestor_event_string, params).rv()


def trim_unaligned_sequences(job, sequences, alignments, params, has_resources=False):

    if not has_resources:
        return job.addChildJobFn(trim_unaligned_sequences, sequences, alignments, params, has_resources=True,
                                 disk=8*sum([seq.size for seq in sequences]) + 32*alignments.size,
                                 memory=cactus_clamp_memory(8*sum([seq.size for seq in sequences]) + 32*alignments.size)).rv()

    work_dir = job.fileStore.getLocalTempDir()
    alignments_file = os.path.join(work_dir, 'alignments.paf')    
    job.fileStore.readGlobalFile(alignments, alignments_file)  # Download the alignments

    # Make a bed of aligned sequence
    # run paffy to_bed to create a bed of aligned coverage
    bed_file = alignments_file + '.bed' # Get a temporary file to store the bed file in
    messages = cactus_call(parameters=['paffy', 'to_bed', "--binary", "--excludeUnaligned", "--includeInverted",
                                       '-i', alignments_file, "--logLevel", getLogLevelString()],
                           outfile=bed_file, returnStdErr=True, job_memory=job.memory)
    # Log paffy to_bed
    logger.info("paffy to_bed:\n{}".format(messages[:-1]))

    trimmed_sequence_files = []
    for i, sequence in enumerate(sequences):
        # run faffy extract to extract unaligned sequences longer than a threshold creating a reduced subset of A
        seq_file = os.path.join(work_dir, '{}.fa'.format(i))
        job.fileStore.readGlobalFile(sequence, seq_file)  # Load original sequence file
        trimmed_seq_file = seq_file + '.trim'  # Get a temporary file to store the extracted subsequences
        messages = cactus_call(parameters=['faffy', 'extract', "-i", bed_file, seq_file, "--skipMissing", "--minSize", "1",
                                           "--flank", params.find("blast").attrib["trimOutgroupFlanking"],
                                           "--logLevel", getLogLevelString()],
                               outfile=trimmed_seq_file, returnStdErr=True, job_memory=job.memory)
        logger.info("faffy extract \n{}".format(messages[:-1]))  # Log faffy extract
        trimmed_sequence_files.append(trimmed_seq_file)

    # Now convert the alignments to refer to the reduced sequences
    trimmed_alignments = alignments_file + '.trim'  # Get a temporary file to store the "trimmed" alignments
    messages = cactus_call(parameters=['paffy', 'upconvert', "-i", alignments_file, "--logLevel", getLogLevelString()] +
                           trimmed_sequence_files, outfile=trimmed_alignments, returnStdErr=True)
    logger.info("paffy upconvert\n{}".format(messages[:-1]))  # Log

    return [job.fileStore.writeGlobalFile(i) for i in trimmed_sequence_files], \
        job.fileStore.writeGlobalFile(trimmed_alignments)  # Return the trimmed sequence files and trimmed alignments

