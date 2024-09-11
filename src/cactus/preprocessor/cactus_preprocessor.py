#!/usr/bin/env python3

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""

"""
import os
import sys
import math
import copy
from argparse import ArgumentParser
import xml.etree.ElementTree as ET

from toil.statsAndLogging import logger

from sonLib.bioio import getTempDirectory
from toil.common import Toil
from toil.job import Job
from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import runGetChunks
from cactus.shared.common import makeURL
from cactus.shared.common import readGlobalFileWithoutCache
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import enableDumpStack
from cactus.shared.common import unzip_gzs
from cactus.shared.common import zip_gzs
from cactus.shared.version import cactus_commit
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger

from cactus.shared.common import cactus_override_toil_options
from cactus.preprocessor.checkUniqueHeaders import checkUniqueHeaders
from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMask import LastzRepeatMaskJob
from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMask import RepeatMaskOptions
from cactus.preprocessor.dnabrnnMasking import DnabrnnMaskJob, loadDnaBrnnModel
from cactus.preprocessor.redMasking import RedMaskJob
from cactus.preprocessor.cutHeaders import CutHeadersJob
from cactus.preprocessor.fileMasking import maskJobOverride, FileMaskingJob
from cactus.progressive.cactus_prepare import human2bytesN

class PreprocessorOptions:
    def __init__(self, chunkSize, memory, cpu, check, proportionToSample, unmask,
                 preprocessJob, checkAssemblyHub=None, lastzOptions=None, minPeriod=None,
                 gpu=0, lastz_memory=None, dnabrnnOpts=None,
                 dnabrnnAction=None, redOpts=None, redPrefilterOpts=None, eventName=None, minLength=None,
                 cutBefore=None, cutBeforeOcc=None, cutAfter=None, inputBedID=None):
        self.chunkSize = chunkSize
        self.memory = memory
        self.cpu = cpu
        self.check = check
        self.proportionToSample=proportionToSample
        self.unmask = unmask
        self.preprocessJob = preprocessJob
        self.checkAssemblyHub = checkAssemblyHub
        self.lastzOptions = lastzOptions
        self.minPeriod = minPeriod
        self.gpu = gpu
        self.gpuLastzInterval = self.chunkSize
        if self.gpu:
            self.chunkSize = 0
        self.lastz_memory= lastz_memory
        self.dnabrnnOpts = dnabrnnOpts
        self.dnabrnnAction = dnabrnnAction
        assert dnabrnnAction in ('softmask', 'hardmask', 'clip')
        self.redOpts = redOpts
        self.redPrefilterOpts = redPrefilterOpts
        self.eventName = eventName
        self.minLength = minLength        
        self.cutBefore = cutBefore
        self.cutBeforeOcc = cutBeforeOcc
        self.cutAfter = cutAfter
        self.inputBedID = inputBedID

class CheckUniqueHeaders(RoundedJob):
    """
    Check that the headers of the input file meet certain naming requirements.
    """
    def __init__(self, prepOptions, inChunkID):
        disk = 2*inChunkID.size
        RoundedJob.__init__(self, memory=prepOptions.memory, cores=prepOptions.cpu, disk=disk,
                     preemptable=True)
        self.prepOptions = prepOptions
        self.inChunkID = inChunkID

    def run(self, fileStore):
        inChunk = fileStore.readGlobalFile(self.inChunkID)
        outChunk = fileStore.getLocalTempFile()
        with open(inChunk) as inFile, open(outChunk, 'w') as outFile:
            checkUniqueHeaders(inFile, outFile, self.prepOptions.eventName, checkAssemblyHub=self.prepOptions.checkAssemblyHub)
        # We re-write the file here so that the output's lifecycle
        # matches the other chunked jobs, which usually write a new
        # chunk. 
        return fileStore.writeGlobalFile(outChunk)

class MergeChunks(RoundedJob):
    def __init__(self, prepOptions, chunkIDList):
        RoundedJob.__init__(self, preemptable=True)
        self.prepOptions = prepOptions
        self.chunkIDList = chunkIDList

    def run(self, fileStore):
        return self.addFollowOn(MergeChunks2(self.prepOptions, self.chunkIDList)).rv()

class MergeChunks2(RoundedJob):
    """merge a list of chunks into a fasta file"""
    def __init__(self, prepOptions, chunkIDList):
        disk = 2*sum([chunkID.size for chunkID in chunkIDList])
        RoundedJob.__init__(self, cores=prepOptions.cpu, memory=prepOptions.memory, disk=disk,
                     preemptable=True)
        self.prepOptions = prepOptions
        self.chunkIDList = chunkIDList

    def run(self, fileStore):
        chunkList = [readGlobalFileWithoutCache(fileStore, fileID) for fileID in self.chunkIDList]

        #Docker expects paths relative to the work dir
        chunkList = [os.path.basename(chunk) for chunk in chunkList]
        outSequencePath = fileStore.getLocalTempFile()
        cactus_call(outfile=outSequencePath, stdin_string=" ".join(chunkList),
                    parameters=["faffy", "merge"])
        for chunkID in self.chunkIDList:
            fileStore.deleteGlobalFile(chunkID)
        return fileStore.writeGlobalFile(outSequencePath)

class PreprocessSequence(RoundedJob):
    """Cut a sequence into chunks, process, then merge
    """
    def __init__(self, prepOptions, inSequenceID, chunksToCompute=None):
        disk = 3*inSequenceID.size if hasattr(inSequenceID, "size") else None
        RoundedJob.__init__(self, memory=prepOptions.memory, disk=disk,
                     preemptable=True)
        self.prepOptions = prepOptions
        self.inSequenceID = inSequenceID
        self.chunksToCompute = chunksToCompute

    def getChunkedJobForCurrentStage(self, seqIDs, proportionSampled, inChunkID, chunk_i):
        """
        Give the chunked work to the appropriate job.
        """
        if self.prepOptions.preprocessJob == "checkUniqueHeaders":
            return CheckUniqueHeaders(self.prepOptions, inChunkID)
        elif self.prepOptions.preprocessJob == "lastzRepeatMask":
            repeatMaskOptions = RepeatMaskOptions(proportionSampled=proportionSampled,
                                                  minPeriod=self.prepOptions.minPeriod,
                                                  lastzOpts=self.prepOptions.lastzOptions,
                                                  gpu=self.prepOptions.gpu,
                                                  cpu=self.prepOptions.cpu,
                                                  lastz_memory=self.prepOptions.lastz_memory,
                                                  gpuLastzInterval=self.prepOptions.gpuLastzInterval,
                                                  eventName='{}_{}'.format(self.prepOptions.eventName, chunk_i))
            return LastzRepeatMaskJob(repeatMaskOptions=repeatMaskOptions,
                                      queryID=inChunkID,
                                      targetIDs=seqIDs)
        elif self.prepOptions.preprocessJob == "dna-brnn":
            return DnabrnnMaskJob(inChunkID,
                                  dnabrnnOpts=self.prepOptions.dnabrnnOpts,
                                  minLength=self.prepOptions.minLength,
                                  action=self.prepOptions.dnabrnnAction,
                                  eventName=self.prepOptions.eventName,
                                  cpu=self.prepOptions.cpu)
        elif self.prepOptions.preprocessJob == "red":
            return RedMaskJob(inChunkID,
                              redOpts=self.prepOptions.redOpts,
                              redPrefilterOpts=self.prepOptions.redPrefilterOpts,
                              eventName=self.prepOptions.eventName,
                              unmask=self.prepOptions.unmask)  
        elif self.prepOptions.preprocessJob == "cutHeaders":
            return CutHeadersJob(inChunkID,
                                 cutBefore=self.prepOptions.cutBefore,
                                 cutBeforeOcc=self.prepOptions.cutBeforeOcc,
                                 cutAfter=self.prepOptions.cutAfter)
        elif self.prepOptions.preprocessJob == 'maskFile':
            return FileMaskingJob(inChunkID,
                                  minLength=self.prepOptions.minLength,
                                  action=self.prepOptions.dnabrnnAction,
                                  inputBedID=self.prepOptions.inputBedID,
                                  eventName=self.prepOptions.eventName)
        else:
            raise RuntimeError("Unknown preprocess job %s" % self.prepOptions.preprocessJob)

    def run(self, fileStore):
        logger.info("Preparing sequence for preprocessing")

        inSequence = fileStore.readGlobalFile(self.inSequenceID)

        if self.prepOptions.chunkSize <= 0:
            # In this first case we don't need to break up the sequence
            chunked = False
            inChunkList = [inSequence]
        else:
            # chunk it up
            chunked = True
            inChunkDirectory = getTempDirectory(rootDir=fileStore.getLocalTempDir())
            inChunkList = runGetChunks(sequenceFiles=[inSequence], chunksDir=inChunkDirectory,
                                       chunkSize=self.prepOptions.chunkSize,
                                       overlapSize=0)
            inChunkList = [os.path.abspath(path) for path in inChunkList]
        logger.info("Chunks = %s" % inChunkList)

        inChunkIDList = [fileStore.writeGlobalFile(chunk, cleanup=True) for chunk in inChunkList]
        outChunkIDList = []
        #For each input chunk we create an output chunk, it is the output chunks that get concatenated together.
        if not self.chunksToCompute:
            self.chunksToCompute = list(range(len(inChunkList)))
        for i in self.chunksToCompute:
            #Calculate the number of chunks to use
            inChunkNumber = int(max(1, math.ceil(len(inChunkList) * self.prepOptions.proportionToSample)))
            assert inChunkNumber <= len(inChunkList) and inChunkNumber > 0
            #Now get the list of chunks flanking and including the current chunk
            j = max(0, i - inChunkNumber//2)
            inChunkIDs = inChunkIDList[j:j+inChunkNumber]
            if len(inChunkIDs) < inChunkNumber: #This logic is like making the list circular
                inChunkIDs += inChunkIDList[:inChunkNumber-len(inChunkIDs)]
            assert len(inChunkIDs) == inChunkNumber
            if self.prepOptions.gpu:
                # when using gpu lastz, we pass through the proportion directly to segalign
                proportionSampled = self.prepOptions.proportionToSample
            else:
                # otherwise, it's taken from the ratio of chunks
                proportionSampled = float(inChunkNumber)/len(inChunkIDList)
            outChunkIDList.append(self.addChild(self.getChunkedJobForCurrentStage(inChunkIDs, proportionSampled, inChunkIDList[i], i)).rv())

        if chunked:
            # Merge results of the chunking process back into a genome-wide file
            return self.addFollowOn(MergeChunks(self.prepOptions, outChunkIDList)).rv()
        else:
            # Didn't chunk--we have a genome-wide fasta file
            return outChunkIDList[0]

def unmaskFasta(inFasta, outFasta):
    """Uppercase a fasta file (removing the soft-masking)."""
    with open(outFasta, 'w') as out:
        for line in open(inFasta):
            if len(line) > 0 and line[0] != ">":
                out.write(line.upper())
            else:
                out.write(line)

class BatchPreprocessor(RoundedJob):
    def __init__(self, prepXmlElems, inSequenceID, iteration = 0):
        self.prepXmlElems = prepXmlElems
        self.inSequenceID = inSequenceID
        self.iteration = iteration
        RoundedJob.__init__(self, preemptable=True)

    def run(self, fileStore):
        # Parse the "preprocessor" config xml element
        assert self.iteration < len(self.prepXmlElems)

        lastIteration = self.iteration == len(self.prepXmlElems) - 1

        prepNode = self.prepXmlElems[self.iteration]
        if getOptionalAttrib(prepNode, "active", typeFn = bool, default=True):
            prepOptions = PreprocessorOptions(chunkSize = int(prepNode.get("chunkSize", default="-1")),
                                              preprocessJob=prepNode.attrib["preprocessJob"],
                                              memory = int(prepNode.get("memory", default=0)),
                                              cpu = int(prepNode.get("cpu", default=1)),
                                              check = bool(int(prepNode.get("check", default="0"))),
                                              proportionToSample = getOptionalAttrib(prepNode, "proportionToSample", typeFn=float, default=1.0),
                                              unmask = getOptionalAttrib(prepNode, "unmask", typeFn=bool, default=False),
                                              lastzOptions = getOptionalAttrib(prepNode, "lastzOpts", default=""),
                                              minPeriod = getOptionalAttrib(prepNode, "minPeriod", typeFn=int, default=0),
                                              checkAssemblyHub = getOptionalAttrib(prepNode, "checkAssemblyHub", typeFn=bool, default=False),
                                              gpu = getOptionalAttrib(prepNode, "gpu", typeFn=int, default=0),
                                              lastz_memory = getOptionalAttrib(prepNode, "lastz_memory", typeFn=int, default=None),
                                              dnabrnnOpts = getOptionalAttrib(prepNode, "dna-brnnOpts", default=""),
                                              dnabrnnAction = getOptionalAttrib(prepNode, "action", typeFn=str, default="softmask"),
                                              redOpts = getOptionalAttrib(prepNode, "redOpts", default=""),
                                              redPrefilterOpts = getOptionalAttrib(prepNode, "redPrefilterOpts", default=""),
                                              eventName = getOptionalAttrib(prepNode, "eventName", typeFn=str, default=None),
                                              minLength = getOptionalAttrib(prepNode, "minLength", typeFn=int, default=1),
                                              cutBefore = getOptionalAttrib(prepNode, "cutBefore", typeFn=str, default=None),
                                              cutBeforeOcc = getOptionalAttrib(prepNode, "cutBeforeOcc", typeFn=int, default=None),
                                              cutAfter = getOptionalAttrib(prepNode, "cutAfter", typeFn=str, default=None),
                                              inputBedID = getOptionalAttrib(prepNode, "inputBedID", typeFn=str, default=None))

            if prepOptions.unmask:
                inSequence = fileStore.readGlobalFile(self.inSequenceID)
                unmaskedInputFile = fileStore.getLocalTempFile()
                unmaskFasta(inSequence, unmaskedInputFile)
                fileStore.deleteGlobalFile(self.inSequenceID)
                self.inSequenceID = fileStore.writeGlobalFile(unmaskedInputFile)

            ppJob = self.addChild(PreprocessSequence(prepOptions, self.inSequenceID))
            outSeqID = ppJob.rv()
            self.addFollowOnJobFn(clean_if_different, self.inSequenceID, outSeqID)
        else:
            logger.info("Skipping inactive preprocessor {}".format(prepNode.attrib["preprocessJob"]))
            outSeqID = self.inSequenceID

        if lastIteration == False:
            return self.addFollowOn(BatchPreprocessor(self.prepXmlElems, outSeqID, self.iteration + 1)).rv()
        else:
            return outSeqID

def clean_if_different(job, file_id, other_file_id):
    """ remove file_id from jobstore if its differetn from other_file_id"""
    if file_id != other_file_id:
        job.fileStore.deleteGlobalFile(file_id)
        
############################################################
############################################################
############################################################
##The preprocessor phase, which modifies the input sequences
############################################################
############################################################
############################################################

class CactusPreprocessor(RoundedJob):
    """Modifies the input genomes, doing things like masking/checking, etc.
    """
    def __init__(self, inputSequenceIDs, configNode, eventNames=[]):
        RoundedJob.__init__(self, disk=sum([id.size for id in inputSequenceIDs if hasattr(id, 'size')]), preemptable=True)
        self.inputSequenceIDs = inputSequenceIDs
        self.configNode = configNode
        self.eventNames = eventNames

    def run(self, fileStore):
        outputSequenceIDs = []
        if self.eventNames:
            assert len(self.eventNames) == len(self.inputSequenceIDs)
            configs = []
            for eventName in self.eventNames:
                conf = copy.deepcopy(self.configNode)
                for node in conf.findall("preprocessor"):
                    node.attrib["eventName"] = eventName
                # if we don't make different configs, the same reference somehow gets passed to mulitple childs below
                configs.append(conf)

        for i, inputSequenceID in enumerate(self.inputSequenceIDs):
            confNode = configs[i] if self.eventNames else self.configNode
            outputSequenceIDs.append(self.addChild(CactusPreprocessor2(inputSequenceID, confNode)).rv())
        return outputSequenceIDs

    @staticmethod
    def getOutputSequenceFiles(inputSequences, outputSequenceDir):
        """Function to get unambiguous file names for each input sequence in the output sequence dir.
        """
        if not outputSequenceDir.startswith("s3://") and not os.path.isdir(outputSequenceDir):
            os.mkdir(outputSequenceDir)
        return [ os.path.join(outputSequenceDir, inputSequences[i].split("/")[-1] + "_%i" % i) for i in range(len(inputSequences)) ]

class CactusPreprocessor2(RoundedJob):
    def __init__(self, inputSequenceID, configNode):
        RoundedJob.__init__(self, preemptable=True)
        self.inputSequenceID = inputSequenceID
        self.configNode = configNode

    def run(self, fileStore):
        prepXmlElems = self.configNode.findall("preprocessor")

        if len(prepXmlElems) == 0: #Just cp the file to the output file
            return self.inputSequenceID
        else:
            logger.info("Adding child batch_preprocessor target")
            return self.addChild(BatchPreprocessor(prepXmlElems, self.inputSequenceID, 0)).rv()

def stageWorkflow(outputSequenceDir, configNode, inputSequences, toil, restart=False,
                  outputSequences = [], maskMode=None, maskAction=None,
                  maskFile=None, minLength=None, inputEventNames=None, brnnCores=None,
                  gpu_override=False, options=None):
    #Replace any constants
    if not outputSequences:
        outputSequences = CactusPreprocessor.getOutputSequenceFiles(inputSequences, outputSequenceDir)
    else:
        assert len(outputSequences) == len(inputSequences)

    # Make sure we have the dna-brnn model in the filestore if we need it
    loadDnaBrnnModel(toil, configNode, maskAlpha = maskMode == 'brnn')
        
    if configNode.find("constants") != None:
        ConfigWrapper(configNode).substituteAllPredefinedConstantsWithLiterals(options)
    if maskMode:
        assert maskMode in ['red', 'lastz', 'brnn', 'none']
        red = maskMode == 'red'
        lastz = maskMode == 'lastz'
        brnn = maskMode == 'brnn'
        ConfigWrapper(configNode).setPreprocessorActive("red", red)
        ConfigWrapper(configNode).setPreprocessorActive("lastzRepeatMask", lastz)
        ConfigWrapper(configNode).setPreprocessorActive("dna-brnn", brnn)
        if brnn and maskAction == 'clip':
            for node in configNode.findall("preprocessor"):
                if getOptionalAttrib(node, "preprocessJob") == 'dna-brnn':
                    node.attrib["action"] = "clip"
    if brnnCores is not None:
        for node in configNode.findall("preprocessor"):
            if getOptionalAttrib(node, "preprocessJob") == 'dna-brnn':
                node.attrib["cpu"] = brnnCores
    if minLength is not None:
        for node in configNode.findall("preprocessor"):
            node.attrib["minLength"] = minLength
    if not restart:
        inputSequenceIDs = []
        for seq in inputSequences:
            logger.info("Importing {}".format(seq))
            inputSequenceIDs.append(toil.importFile(makeURL(seq)))
        maskFileID = toil.importFile(makeURL(maskFile)) if maskFile else None
        unzip_job = Job.wrapJobFn(unzip_then_pp, configNode, inputSequences, inputSequenceIDs, inputEventNames,
                                  maskFile, maskFileID, maskAction, minLength)
        outputSequenceIDs = toil.start(unzip_job)
    else:
        outputSequenceIDs = toil.restart()
    for seqID, path in zip(outputSequenceIDs, outputSequences):
        try:
            iter(seqID)
            # dna-brnn and maskFile will output a couple of bed files.  we scrape those out here
            toil.exportFile(seqID[0], makeURL(path))
            toil.exportFile(seqID[1], makeURL(path) + '.bed')
            toil.exportFile(seqID[2], makeURL(path) + '.mask.bed')
        except:
            toil.exportFile(seqID, makeURL(path))

def unzip_then_pp(job, config_node, input_fa_paths, input_fa_ids, input_event_names, mask_file_path, mask_file_id, mask_file_action, min_length):
    """ unzip then preprocess """
    unzip_job = job.addChildJobFn(unzip_gzs, input_fa_paths, input_fa_ids)
    if mask_file_id is not None:
        mask_unzip_job = unzip_job.addChildJobFn(unzip_gzs, [mask_file_path], [mask_file_id])
        config_node = mask_unzip_job.addFollowOnJobFn(maskJobOverride, config_node, mask_file_path, mask_unzip_job.rv(0), mask_file_action, min_length,
                                                      disk=mask_file_id.size*20).rv()
    pp_job = unzip_job.addFollowOn(CactusPreprocessor([unzip_job.rv(i) for i in range(len(input_fa_ids))], config_node, eventNames=input_event_names))
    zip_job = pp_job.addFollowOnJobFn(zip_gzs, input_fa_paths,  pp_job.rv(), list_elems = [0])
    return zip_job.rv()
    
def runCactusPreprocessor(outputSequenceDir, configFile, inputSequences, toilDir):
    toilOptions = Job.Runner.getDefaultOptions(toilDir)
    toilOptions.logLevel = "INFO"
    toilOptions.disableCaching = True
    with Toil(toilOptions) as toil:
        stageWorkflow(outputSequenceDir, ET.parse(options.configFile).getroot(), inputSequences, toil, options=toilOptions)

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument("inSeqFile", type=str, nargs='?', default=None, help = "Input Seq file")
    parser.add_argument("outSeqFile", type=str, nargs='?', default=None, help = "Output Seq file (ex generated with cactus-prepare)")
    parser.add_argument("--configFile", default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--inputNames", nargs='*', help='input genome names (not paths) to preprocess (all leaves from Input Seq file if none specified)')
    parser.add_argument("--inPaths", nargs='*', help='Space-separated list of input fasta paths (to be used in place of --inSeqFile')
    parser.add_argument("--outPaths", nargs='*', help='Space-separated list of output fasta paths (one for each inPath, used in place of --outSeqFile)')
    parser.add_argument("--maskMode", type=str, help='Masking mode, one of {"lastz", "brnn", "red", "none", "config"}. Default="config".', default='config')
    parser.add_argument("--maskAction", type=str, help='Masking action, one of {"softmask", "hardmask", "clip"}. Default="softmask"', default='softmask')
    parser.add_argument("--minLength", type=int, help='Minimum interval threshold for masking.  Overrides config')
    parser.add_argument("--maskFile", type=str, help='Add masking from BED or PAF file to sequences, **ignoring all other preprocessors**')
    parser.add_argument("--ignore", nargs='*', help='Space-separate list of genomes from inSeqFile to ignore', default=[])
    parser.add_argument("--brnnCores", type=int, help='Specify number of cores for each dna-brnn job (overriding default value from the config)')
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

    parser.add_argument("--pangenome", action="store_true", help='Do not mask. Just add Cactus-style unique prefixes and strip anything up to and including last #')

    options = parser.parse_args()
    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    # we have two modes: operate directly on paths or rely on the seqfiles.  they cannot be mixed
    if options.inSeqFile or options.outSeqFile:
        if not options.inSeqFile or not options.outSeqFile or options.inPaths or options.outPaths:
            raise RuntimeError('--inSeqFile must be used in conjunction with --outSeqFile and not with --inPaths nor --outPaths')
    elif options.inPaths or options.outPaths:
        if not options.inPaths or not options.outPaths or options.inSeqFile or options.outSeqFile:
            raise RuntimeError('--inPaths must be used in conjunction with --outPaths and not with --inSeqFile nor --outSeqFile')
        if len(options.inPaths) != len(options.outPaths):
            raise RuntimeError('--inPaths and --outPaths must have the same number of arguments')
        if not options.inputNames or len(options.inputNames) != len(options.inPaths):
            raise RuntimeError('--inputNames must be used in conjunction with --inputPaths to specify (exactly) one event name for each input')
    else:
        raise RuntimeError('--inSeqFile/--outSeqFile/--inputNames or --inPaths/--outPaths required to specify input')
    if options.maskMode not in ['lastz', 'brnn', 'red', 'none', 'config']:
        raise RuntimeError('--maskMode must be one of {"lastz", "brnn", "red", "none", "config"}')    
    if options.maskAction not in ['softmask', 'hardmask', 'clip']:
        raise RuntimeError('--maskAction must be one of {"softmask", "hardmask", "clip"}')
    if options.maskMode == 'lastz' and options.maskAction != 'softmask':
        raise RuntimeError('only softmasking is supported with lastz')
    if options.maskMode == 'red' and options.maskAction != 'softmask':
        raise RuntimeError('only softmasking is supported with red')    
    if options.maskFile and options.maskFile.endswith('.paf') and not options.inputNames and not options.inSeqFile:
        raise RuntimeError('paf masking requires event names specified wither with an input seqfile or with --inputNames')
    if options.maskFile and options.minLength is None:
        raise RuntimeError('--minLength must be used with --maskFile')
    
    inSeqPaths = []
    outSeqPaths = []
    inNames = options.inputNames

    #load cactus config
    configNode = ET.parse(options.configFile).getroot()
    #we never want to preprocess minigraph sequences
    graph_event = getOptionalAttrib(findRequiredNode(configNode, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    options.ignore.append(graph_event)

    # toggle on the gpu
    config_wrapper = ConfigWrapper(configNode)
    config_wrapper.initLastz(options)

    # apply pangenome overrides
    if options.pangenome:
        # don't mask
        options.maskMode = 'none'
        # set cutBefore to #
        # so HG002#0#chrY would turn into id=HG002.0|chrY (so long as event name is HG002.0)
        for node in configNode.findall("preprocessor"):
            if getOptionalAttrib(node, "preprocessJob") == 'cutHeaders':
                node.attrib["cutBefore"] = "#"
    
    # mine the paths out of the seqfiles
    if options.inSeqFile:
        inSeqFile = SeqFile(options.inSeqFile, defaultBranchLen=config_wrapper.getDefaultBranchLen(options.pangenome))
        outSeqFile = SeqFile(options.outSeqFile, defaultBranchLen=config_wrapper.getDefaultBranchLen(options.pangenome))

        if not inNames:
            inNames = [inSeqFile.tree.getName(node) for node in inSeqFile.tree.getLeaves()]

        for inName in inNames:
            if inName in options.ignore:
                # "convenience" functionality: we let the --ignore option update the output seqfile
                # to reflect the fact that we're not touching the original input
                outSeqFile.pathMap[inName] = inSeqFile.pathMap[inName]
                continue
            if inName not in inSeqFile.pathMap or inName not in outSeqFile.pathMap:
                raise RuntimeError('{} not present in input and output Seq files'.format(inName))
            inPath = inSeqFile.pathMap[inName]
            outPath = outSeqFile.pathMap[inName]
            if os.path.isdir(inPath):
                try:
                    os.makedirs(outPath)
                except:
                    pass
                assert os.path.isdir(inPath) == os.path.isdir(outPath)
                inSeqPaths += [os.path.join(inPath, seqPath) for seqPath in os.listdir(inPath)]
                outSeqPaths += [os.path.join(outPath, seqPath) for seqPath in os.listdir(inPath)]
            else:
                inSeqPaths += [inPath]
                outSeqPaths += [outPath]

        if options.ignore:
            # see comment above
            with open(options.outSeqFile, 'w') as outSF:
                outSF.write(str(outSeqFile))

    # we got path names directly from the command line
    else:
        inSeqPaths = options.inPaths
        outSeqPaths = options.outPaths

    assert outSeqPaths

    if options.ignore and inNames:
        for ignore_event in options.ignore:
            if ignore_event in inNames:
                del inNames[inNames.index(ignore_event)]

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    
    with Toil(options) as toil:
        stageWorkflow(outputSequenceDir=None,
                      configNode=configNode,
                      inputSequences=inSeqPaths,
                      toil=toil,
                      restart=options.restart,
                      outputSequences=outSeqPaths,
                      maskMode=options.maskMode if options.maskMode != 'config' else None,
                      maskAction=options.maskAction,
                      minLength=options.minLength,
                      inputEventNames=inNames,
                      brnnCores=options.brnnCores,
                      options=options)

if __name__ == '__main__':
    main()
