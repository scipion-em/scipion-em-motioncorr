# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [4]
# *
# * [1] SciLifeLab, Stockholm University
# * [4] MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# ******************************************************************************

import os
import time
from math import ceil, sqrt
from threading import Thread
from uuid import uuid4
from collections import OrderedDict
from datetime import datetime

from emtools.utils import Pipeline

import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.gui.plotter import Plotter
from pyworkflow.protocol import STEPS_SERIAL
from pwem.emlib.image import ImageHandler
from pwem.objects import Image, Float, SetOfMovies
from pwem.protocols import ProtAlignMovies

from .. import Plugin
from ..constants import *
from ..convert import *
from .protocol_motioncorr import ProtMotionCorr


class ProtMotionCorrTasks(ProtMotionCorr):
    """ This protocol wraps motioncor2 movie alignment program developed at UCSF.

    Motioncor2 performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'motioncor2 tasks'

    def __init__(self, **kwargs):
        ProtMotionCorr.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL

    # We are not using the steps mechanism for parallelism from Scipion
    def _stepsCheck(self):
        pass

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):
        ProtMotionCorr._defineAlignmentParams(self, form)

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self.samplingRate = self.inputMovies.get().getSamplingRate()
        self._insertFunctionStep('_convertInputStep')
        self._insertFunctionStep('_processAllMoviesStep')

    def _loadInputMovies(self):
        inputMovies = OrderedDict()
        moviesFile = self.getInputMovies().getFileName()
        movieSet = SetOfMovies(filename=moviesFile)
        movieSet.loadAllProperties()
        for m in movieSet:
            mid = m.getObjId()
            if mid not in inputMovies:
                inputMovies[mid] = m.clone()
        movieSet.close()
        return inputMovies

    def _processAllMoviesStep(self):
        self.program = Plugin.getProgram()
        self.command = self._getCmd()

        mc = Pipeline()
        g = mc.addGenerator(self._generateBatches)
        gpus = self.getGpuList()
        mc1 = mc.addProcessor(g.outputQueue, self._getMcProcessor(gpus[0]))
        for gpu in gpus[1:]:
            mc.addProcessor(g.outputQueue, self._getMcProcessor(gpu),
                            outputQueue=mc1.outputQueue)
        moveQueue = mc1.outputQueue

        mc.run()

    def _generateBatches(self):
        """ Check for new input movies and generate new tasks. """
        # FIXME: Take streaming into account
        def _createBatch(movies):
            batch_id = str(uuid4())
            batch_path = self._getTmpPath(batch_id)
            pwutils.cleanPath(batch_path)
            os.mkdir(batch_path)

            for movie in movies:
                fn = movie.getFileName()
                baseName = os.path.basename(fn)
                self._getPath
                os.symlink(os.path.abspath(fn),
                           os.path.join(batch_path, baseName))
            return {
                'movies': movies,
                'id': batch_id,
                'path': batch_path
            }

        inputMovies = self._loadInputMovies()
        movies = []
        for movie in inputMovies.values():
            movies.append(movie)

            if len(movies) == 16:
                yield _createBatch(movies)
                movies = []

        if movies:
            yield _createBatch(movies)

    def _getMcProcessor(self, gpu):
        def _processBatch(batch):
            try:
                cmd = self.command.replace('-Gpu #', f'-Gpu {gpu}')
                self.runJob(self.program, cmd, cwd=batch['path'])

            except Exception as e:
                self.error("ERROR: Motioncor2 has failed for batch %s. --> %s\n"
                           % (batch['id'], str(e)))
                import traceback
                traceback.print_exc()

        return _processBatch

    # --------------------------- INFO functions ------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _getCmd(self):
        """ Set return a command string that will be used for each batch. """
        inputMovies = self.getInputMovies()
        argsDict = self._getMcArgs()
        print("gpu: ", argsDict['-Gpu'])
        argsDict['-Gpu'] = '#'
        argsDict['-Serial'] = 1

        # FIXME: Duplicated from ProtMotionCorr._processMovie
        if self.isEER:
            argsDict['-EerSampling'] = self.eerSampling.get() + 1
            argsDict['-FmIntFile'] = "../../extra/FmIntFile.txt"
        elif self.doApplyDoseFilter:
            preExp, dose = self._getCorrectedDose(inputMovies)
            argsDict.update({'-FmDose': dose,
                             '-InitDose': preExp if preExp > 0.001 else 0})

        if inputMovies.getGain():
            argsDict.update({'-Gain': f'"{inputMovies.getGain()}"',
                             '-RotGain': self.gainRot.get(),
                             '-FlipGain': self.gainFlip.get()})

        if inputMovies.getDark():
            argsDict['-Dark'] = inputMovies.getDark()

        # Get input format, but for the batch
        firstMovie = inputMovies.getFirstItem()
        ext = pwutils.getExt(firstMovie.getFileName()).lower()
        if ext in ['.mrc', '.mrcs']:
            inprefix = '-InMrc'
        elif ext in ['.tif', '.tiff']:
            inprefix = '-InTiff'
        else:
            raise Exception(f"Unsupported format '{ext}' for batch processing "
                            f"in Motioncor protocol. ")

        argsDict[inprefix] = './'
        argsDict['-OutMrc'] = 'output_'

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += self.extraParams2.get()

        return cmd
