#!/usr/bin/env python3
"""Uses dna-brnn to mask alpha satellites with a given length threshold
"""

import os
import re
import sys
import shutil

from toil.lib.threading import cpu_count

from sonLib.bioio import catFiles

from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from cactus.shared.common import cactusRootPath
from toil.realtimeLogger import RealtimeLogger

class DnabrnnMaskJob(RoundedJob):
    def __init__(self, fastaID, dnabrnnOpts, hardmask, cpu, minLength=None):
        memory = 4*1024*1024*1024
        disk = 2*(fastaID.size)
        cores = min(cpu_count(), cpu)
        RoundedJob.__init__(self, memory=memory, disk=disk, cores=cores, preemptable=True)
        self.fastaID = fastaID
        self.minLength = minLength
        self.dnabrnnOpts = dnabrnnOpts
        self.hardmask = hardmask

    def run(self, fileStore):
        """
        mask alpha satellites with dna-brnn
        """
        work_dir = fileStore.getLocalTempDir()
        fastaFile = os.path.join(work_dir, 'seq.fa')
        fileStore.readGlobalFile(self.fastaID, fastaFile)

        cmd = ['dna-brnn', fastaFile] + self.dnabrnnOpts.split()
        
        if '-i' not in self.dnabrnnOpts:
            # pull up the model
            # todo: is there are more robust way?
            model_path = os.path.join(work_dir, 'model.knm')
            embedded_model = os.path.join(cactusRootPath(), 'attcc-alpha.knm')
            # we copy it over for container purposes
            shutil.copyfile(embedded_model, model_path)
            cmd += ['-i', model_path]        
        
        if self.cores:
            cmd += ['-t', str(self.cores)]

        bedFile = fileStore.getLocalTempFile()

        # run dna-brnn to make a bed file
        cactus_call(outfile=bedFile, parameters=cmd)

        maskedFile = fileStore.getLocalTempFile()

        mask_cmd = ['cactus_fasta_softmask_intervals.py', '--origin=zero', bedFile]
        if self.minLength:
            mask_cmd += '--minLength={}'.format(self.minLength)

        if self.hardmask:
            mask_cmd += ['--mask=N']

        # do the softmasking
        cactus_call(infile=fastaFile, outfile=maskedFile, parameters=mask_cmd)

        return fileStore.writeGlobalFile(maskedFile)


