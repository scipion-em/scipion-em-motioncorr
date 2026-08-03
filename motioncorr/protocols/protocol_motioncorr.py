# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Vahid Abrishami (vabrishami@cnb.csic.es) [2]
# *              Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [3]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [4]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [3] Department of Anatomy and Cell Biology, McGill University
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
from os.path import exists
from threading import Thread

import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.gui.plotter import Plotter
from pwem.objects import Image, Float
from pwem.protocols import ProtAlignMovies
from pyworkflow.utils import weakImport

from .. import Plugin, DEFECTS_FILE_EER
from .protocol_base import ProtMotionCorrBase


class ProtMotionCorr(ProtMotionCorrBase, ProtAlignMovies):
    """
    Corrects beam-induced motion in cryo-EM movies and produces aligned
    micrographs suitable for downstream image processing. The protocol
    improves image quality by compensating for frame-to-frame specimen
    drift and by preparing consistent micrograph outputs for further
    structural analysis.

    AI Generated:

    Motion Correction (ProtMotionCorr) - User Manual
        Overview

        The Motion Correction protocol is designed to process cryo-EM movie
        stacks acquired as sequences of individual frames. During electron
        exposure, particles often experience beam-induced movement that
        causes blurring and reduces high-resolution information. The main
        objective of this protocol is to compensate for that movement so
        that the final micrographs better represent the underlying
        biological specimen.

        For biological users, this step is one of the earliest and most
        important stages of single-particle and tomography workflows.
        Correcting motion before contrast transfer function estimation,
        particle picking, or reconstruction improves the reliability of all
        downstream measurements. In practical terms, motion correction often
        determines how much structural detail can ultimately be recovered.

        Inputs and Biological Context

        The protocol operates on imported movies rather than on already
        averaged micrographs. Each movie contains multiple frames that
        capture the progressive electron exposure of the sample. The input
        data are expected to include acquisition metadata such as sampling
        rate, electron dose, and frame organization, since these values are
        biologically important for accurate correction and proper exposure
        interpretation.

        In many experiments, movies may also include detector gain
        references or defect information. These auxiliary inputs help
        preserve image fidelity by compensating for detector-specific
        imperfections. For biological interpretation, reliable detector
        correction is essential because fixed artifacts can otherwise mimic
        structural features or bias downstream analysis.

        Dose Weighting and Radiation Damage

        A central biological feature of this protocol is dose weighting.
        During electron exposure, later frames accumulate more radiation
        damage than earlier ones. Dose weighting accounts for this effect
        by giving stronger emphasis to the less damaged signal.

        In practical cryo-EM workflows, dose weighting usually improves the
        preservation of high-resolution structural information. This is
        particularly valuable when studying macromolecular assemblies where
        fine conformational details matter. When dose information is not
        available or not reliable, this feature should be used with caution
        because biologically meaningful weighting depends directly on
        accurate exposure metadata.

        Global and Local Motion Correction

        Biological specimens do not always move as rigid bodies. Some
        micrographs exhibit approximately uniform drift, while others show
        local deformations caused by beam-induced bending of ice or sample
        support instability.

        This protocol supports both global and local correction strategies.
        Global correction is generally sufficient for stable samples with
        limited deformation. Local correction becomes especially important
        for larger fields of view, thinner ice, or highly flexible samples
        where different regions may move differently during exposure.

        For biological users, local correction often improves the quality
        of particles located near the edges of the image or in regions
        affected by anisotropic drift. However, overly aggressive local
        correction can become unstable when signal is weak, so moderate
        settings are often preferable for exploratory analyses.

        Frame Grouping and EER Data

        The protocol can combine frames into larger fractions before
        alignment. This is often useful when very high frame rates produce
        extremely low signal per frame. Grouping can stabilize alignment
        while keeping the overall exposure structure biologically
        meaningful.

        Modern detectors may also produce EER movies, which represent a
        different acquisition strategy based on very fine temporal
        sampling. For biological applications, EER handling is especially
        useful when precise dose fractionation is required, such as in
        high-resolution single-particle work where exposure modeling
        strongly influences the final map quality.

        Output Products and Interpretation

        The primary result of the protocol is a set of aligned micrographs
        that can be used directly for downstream cryo-EM processing.
        Depending on the selected options, the protocol may also generate
        dose-weighted micrographs, non-weighted micrographs, aligned movie
        stacks, and separated odd-even outputs.

        These alternative outputs serve different biological purposes.
        Dose-weighted micrographs are usually preferred for reconstruction,
        while non-weighted micrographs can be useful for contrast transfer
        function estimation. Odd-even outputs are valuable when validating
        data consistency or preparing independent half-data analyses.

        Motion Diagnostics and Quality Assessment

        Beyond producing corrected micrographs, the protocol also provides
        motion-related diagnostic information. This includes global motion
        trajectories, accumulated displacement estimates, and optional
        power spectrum or thumbnail visualizations.

        For biological interpretation, these diagnostics are highly
        informative. Excessive motion may indicate poor ice quality,
        unstable support films, charging effects, or problematic exposure
        conditions. Reviewing motion behavior early can prevent large
        downstream processing efforts on low-quality datasets.

        Practical Recommendations

        In routine cryo-EM practice, it is usually advisable to begin with
        dose weighting enabled and moderate local correction. For stable
        samples, default settings are often sufficient. For very noisy
        datasets or unusually high frame-rate acquisitions, modest frame
        grouping can improve robustness.

        If local motion appears unstable, reducing correction complexity
        often provides more reliable biological results than attempting
        highly aggressive refinement. When preparing data for high-
        resolution reconstruction, visual inspection of corrected
        micrographs and motion trajectories remains an essential quality
        control step.

        Final Perspective

        For most cryo-EM users, motion correction is not simply a technical
        preprocessing operation. It is the stage where raw detector movies
        are transformed into biologically interpretable images. Accurate
        correction of drift, careful treatment of dose effects, and
        thoughtful evaluation of motion behavior together establish the
        foundation for reliable structural analysis.
    """

    _label = 'movie alignment'
    evenOddCapable = True

    def __init__(self, **kwargs):
        ProtAlignMovies.__init__(self, **kwargs)
        ProtMotionCorrBase.__init__(self, **kwargs)

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif', '.eer', '.gain'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):
        form.addHidden('doSaveAveMic', params.BooleanParam,
                       default=True)
        form.addHidden('useAlignToSum', params.BooleanParam,
                       default=True)

        group = form.addGroup('Alignment')
        line = group.addLine('Frames to ALIGN and SUM',
                             help='Frames range to ALIGN and SUM on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to align, it means that you will '
                                  'align until the last frame of the movie. '
                                  'When using EER, this option is IGNORED!')
        line.addParam('alignFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('alignFrameN', params.IntParam, default=0,
                      label='to')

        group.addParam('binFactor', params.FloatParam, default=1.,
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')

        line = group.addLine('Crop offsets (px)',
                             expertLevel=cons.LEVEL_ADVANCED)
        line.addParam('cropOffsetX', params.IntParam, default=0, label='X',
                      expertLevel=cons.LEVEL_ADVANCED)
        line.addParam('cropOffsetY', params.IntParam, default=0, label='Y',
                      expertLevel=cons.LEVEL_ADVANCED)

        line = group.addLine('Crop dimensions (px)',
                             help='How many pixels to crop from offset\n'
                                  'If equal to 0, use maximum size.',
                             expertLevel=cons.LEVEL_ADVANCED)
        line.addParam('cropDimX', params.IntParam, default=0, label='X',
                      expertLevel=cons.LEVEL_ADVANCED)
        line.addParam('cropDimY', params.IntParam, default=0, label='Y',
                      expertLevel=cons.LEVEL_ADVANCED)

        if self.evenOddCapable:
            form.addParam('splitEvenOdd', params.BooleanParam,
                          default=False,
                          label='Split & sum odd/even frames?',
                          help='Generate odd and even sums using odd and even frames '
                               'respectively when this option is enabled.')

        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Save aligned movie?")

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Compute PSD?",
                      help="If Yes, the protocol will compute for each "
                           "aligned micrograph the PSD using EMAN2.")

        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      expertLevel=cons.LEVEL_ADVANCED,
                      default=False,
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail with EMAN2 and keep it with the '
                           'micrograph object for visualization purposes. ')

        form.addParam('extraProtocolParams', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional protocol parameters',
                      help="Here you can provide some extra parameters for the "
                           "protocol, not the underlying motioncor program."
                           "You can provide many options separated by space. "
                           "\n\n*Options:* \n\n"
                           "--dont_use_worker_thread \n"
                           " Now by default we use a separate thread to compute"
                           " PSD and thumbnail (if is required). This allows "
                           " more effective use of the GPU card, but requires "
                           " an extra CPU. Use this option (NOT RECOMMENDED) if "
                           " you want to prevent this behaviour")

        self._defineCommonParams(form)

    # --------------------------- STEPS functions -----------------------------
    def _processMovie(self, movie):
        inputMovies = self.getInputMovies()
        movieFolder = self._getOutputMovieFolder(movie)
        outputMicFn = self._getOutputMicName(movie)

        argsDict = self._getMcArgs()
        argsDict['-OutMrc'] = f'"{outputMicFn}"'
        argsDict['-LogDir'] = './'
        if self.isEER:
            argsDict.update({'-FmIntFile': "../../extra/FmIntFile.txt"})
        if self.defectFile.get():
            argsDict['-DefectFile'] = self.defectFile.get()
        elif self.defectMap.get():
            argsDict['-DefectMap'] = self.defectMap.get()
        elif exists(self._getExtraPath(DEFECTS_FILE_EER)):
            argsDict['-DefectFile'] = f"../../extra/{DEFECTS_FILE_EER}"
        args = self._getInputFormat(movie.getFileName())
        args += ' '.join(['%s %s' % (k, v)
                          for k, v in argsDict.items()])
        args += ' ' + self.extraParams2.get()

        try:
            self.runJob(Plugin.getProgram(), args, cwd=movieFolder,
                        env=Plugin.getEnviron())
            self._moveOutput(movie)

            def _extraWork():
                # we need to move shifts log to extra dir before parsing
                try:
                    self._saveAlignmentPlots(movie, inputMovies.getSamplingRate())
                    outMicFn = self._getExtraPath(self._getMicFn(movie))

                    if self.doComputePSD:
                        self._computePSD(outMicFn, outputFn=self._getPsdCorr(movie))

                    if self._doComputeMicThumbnail():
                        self.computeThumbnail(
                            outMicFn, outputFn=self._getOutputMicThumbnail(movie))
                except:
                    self.error(f"ERROR: Extra work (i.e plots, PSD, thumbnail) "
                               f"has failed for {movie.getFileName()}\n")

            if self._useWorkerThread():
                thread = Thread(target=_extraWork)
                thread.start()
            else:
                _extraWork()

        except Exception as e:
            self.error(f"ERROR: Motioncor has failed for {movie.getFileName()} --> {str(e)}\n")
            import traceback
            traceback.print_exc()
        
    def _insertFinalSteps(self, deps):
        stepsId = []
        if self._useWorkerThread():
            stepsId.append(self._insertFunctionStep('waitForThreadStep',
                                                    prerequisites=deps,
                                                    needsGPU=False))
        return stepsId

    def waitForThreadStep(self):
        # Quick and dirty (maybe desperate) way to wait
        # if the PSD and thumbnail were computed in a thread
        # If running in streaming this will not be necessary
        time.sleep(10)  # wait 10 sec for the thread to finish

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, 'outputMicrographs') or \
                hasattr(self, 'outputMicrographsDoseWeighted'):
            summary.append('Aligned %d movies using Motioncor.'
                           % self.getInputMovies().getSize())
            if self.splitEvenOdd and self._createOutputWeightedMicrographs():
                summary.append('Even/odd outputs are dose-weighted!')
        else:
            summary.append('Output is not ready')

        return summary

    def _validate(self):
        errors = ProtMotionCorrBase._validate(self)

        # check eman2 plugin
        if self.doComputeMicThumbnail or self.doComputePSD:
            try:
                from pwem import Domain
                _ = Domain.importFromPlugin('eman2', doRaise=True)
            except:
                errors.append("EMAN2 plugin not found!\nComputing thumbnails "
                              "or PSD requires EMAN2 plugin and binaries installed.")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getCwdPath(self, movie, path):
        return os.path.join(self._getOutputMovieFolder(movie), path)

    def _getMovieLogFile(self, movie):
        usePatches = self.patchX != 0 or self.patchY != 0
        return '%s%s-Full.log' % (self._getMovieRoot(movie),
                                  '-Patch' if usePatches else '')

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_psd', 'png', extra=True)

    def getInputMovies(self):
        return self.inputMovies.get()

    def _computePSD(self, inputFn, outputFn, scaleFactor=6):
        """ Generate a thumbnail of the PSD with EMAN2"""
        args = f"{inputFn} {outputFn} "
        args += f"--process=math.realtofft --meanshrink {scaleFactor} "
        args += "--fixintscaling=sane"

        from pwem import Domain
        eman2 = Domain.importFromPlugin('eman2')
        from pyworkflow.utils.process import runJob
        runJob(self._log, eman2.Plugin.getProgram('e2proc2d.py'), args,
               env=eman2.Plugin.getEnviron())

        return outputFn

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setPlotInfo(movie, mic)

    def _getMicFn(self, movie):
        if self.doApplyDoseFilter:
            return self._getOutputMicWtName(movie)
        else:
            return self._getOutputMicName(movie)

    def _setPlotInfo(self, movie, mic):
        mic.plotGlobal = Image(location=self._getPlotGlobal(movie))
        if self.doComputePSD:
            mic.psdCorr = Image(location=self._getPsdCorr(movie))
        if self._doComputeMicThumbnail():
            mic.thumbnail = Image(location=self._getOutputMicThumbnail(movie))
        if self.doApplyDoseFilter:
            total, early, late = self.calcFrameMotion(movie)
            mic._rlnAccumMotionTotal = Float(total)
            mic._rlnAccumMotionEarly = Float(early)
            mic._rlnAccumMotionLate = Float(late)

    def _saveAlignmentPlots(self, movie, pixSize):
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movie)
        first, _ = self._getFramesRange()
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first, pixSize)
        plotter.savefig(self._getPlotGlobal(movie))
        plotter.close()

    def _moveOutput(self, movie):
        """ Move output from tmp to extra folder. """
        def _moveToExtra(fn):
            """ Move file from movies tmp folder to extra """
            if os.path.exists(self._getCwdPath(movie, fn)):
                pwutils.moveFile(self._getCwdPath(movie, fn),
                                 self._getExtraPath(fn))

        _moveToExtra(self._getMovieLogFile(movie))
        _moveToExtra(self._getMicFn(movie))
        _moveToExtra(pwutils.replaceExt(movie.getBaseName(), "star"))

        if self._doSaveUnweightedMic():
            _moveToExtra(self._getOutputMicName(movie))

        if self.splitEvenOdd:
            _moveToExtra(self._getOutputMicEvenName(movie))
            _moveToExtra(self._getOutputMicOddName(movie))

        if self.doSaveMovie:
            outputMicFn = self._getCwdPath(movie, self._getOutputMicName(movie))
            outputMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
            movieFn = outputMicFn.replace('_aligned_mic.mrc',
                                          '_aligned_mic_Stk.mrc')
            pwutils.moveFile(movieFn, outputMovieFn)

    def _updateOutputSet(self, outputName, outputSet,
                         state=pwobj.Set.STREAM_OPEN):
        """ Redefine this method to set EER attrs. """
        first = getattr(self, '_firstUpdate', True)

        if first and outputName == 'outputMovies':
            with weakImport('relion'):
                from relion.convert.convert31 import OpticsGroups
                og = OpticsGroups.fromImages(outputSet)
                gain = self.getInputMovies().getGain()
                ogDict = {'rlnMicrographStartFrame': self.alignFrame0.get()}
                if self.isEER:
                    ogDict.update({'rlnEERGrouping': self.eerGroup.get(),
                                   'rlnEERUpsampling': self.eerSampling.get() + 1})
                if gain:
                    ogDict['rlnMicrographGainName'] = gain
                og.updateAll(**ogDict)
                og.toImages(outputSet)

        ProtAlignMovies._updateOutputSet(self, outputName, outputSet,
                                         state=state)
        self._firstUpdate = False

    def _createOutputMicrographs(self):
        createWeighted = self._createOutputWeightedMicrographs()
        # To create the unweighted average micrographs
        # we only consider the 'doSaveUnweightedMic' flag if the
        # weighted ones should be created.
        return not createWeighted or self.doSaveUnweightedMic

    def _doSaveUnweightedMic(self):
        """ Wraps the logic for saving unweighted mics that needs to consider the doApplyDoseFilter to be true"""
        return self._createOutputWeightedMicrographs() and self.doSaveUnweightedMic

    def _createOutputWeightedMicrographs(self):
        return self.doApplyDoseFilter

    def _doComputeMicThumbnail(self):
        return self.doComputeMicThumbnail

    def _useWorkerThread(self):
        return '--dont_use_worker_thread' not in self.extraProtocolParams.get()

    def getSamplingRate(self):
        return self.getInputMovies().getSamplingRate()

    def calcFrameMotion(self, movie):
        # based on relion 3.1 motioncorr_runner.cpp
        shiftsX, shiftsY = self._getMovieShifts(movie)
        a0, aN = self._getFramesRange()
        nframes = aN - a0 + 1
        preExp, dose = self._getCorrectedDose()
        # when using EER, the hardware frames are grouped
        if self.isEER:
            dose *= self.eerGroup.get()
        cutoff = (4 - preExp) // dose  # early is <= 4e/A^2
        total, early, late = 0., 0., 0.
        x, y, xOld, yOld = 0., 0., 0., 0.
        pix = self.getSamplingRate()
        try:
            for frame in range(2, nframes + 1):  # start from the 2nd frame
                x, y = shiftsX[frame - 1], shiftsY[frame - 1]
                d = sqrt((x - xOld) * (x - xOld) + (y - yOld) * (y - yOld))
                total += d
                if frame <= cutoff:
                    early += d
                else:
                    late += d
                xOld = x
                yOld = y
            return list(map(lambda x: pix * x, [total, early, late]))
        except IndexError:
            self.error(f"Expected {nframes} frames, found less. "
                       f"Check movie {movie.getFileName()}")


def createGlobalAlignmentPlot(meanX, meanY, first, pixSize):
    """ Create a plotter with the shift per frame. """
    sumMeanX = []
    sumMeanY = []

    def px_to_ang(px):
        y1, y2 = px.get_ylim()
        x1, x2 = px.get_xlim()
        ax_ang2.set_ylim(y1 * pixSize, y2 * pixSize)
        ax_ang.set_xlim(x1 * pixSize, x2 * pixSize)
        ax_ang.figure.canvas.draw()
        ax_ang2.figure.canvas.draw()

    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax_px = figure.add_subplot(111)
    ax_px.grid()
    ax_px.set_xlabel('Shift x (px)')
    ax_px.set_ylabel('Shift y (px)')

    ax_ang = ax_px.twiny()
    ax_ang.set_xlabel('Shift x (A)')
    ax_ang2 = ax_px.twinx()
    ax_ang2.set_ylabel('Shift y (A)')

    i = first
    # The output and log files list the shifts relative to the first frame.
    # ROB unit seems to be pixels since sampling rate is only asked
    # by the program if dose filtering is required
    skipLabels = ceil(len(meanX) / 10.0)
    labelTick = 1

    for x, y in zip(meanX, meanY):
        sumMeanX.append(x)
        sumMeanY.append(y)
        if labelTick == 1:
            ax_px.text(x - 0.02, y + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    # automatically update lim of ax_ang when lim of ax_px changes.
    ax_px.callbacks.connect("ylim_changed", px_to_ang)
    ax_px.callbacks.connect("xlim_changed", px_to_ang)

    ax_px.plot(sumMeanX, sumMeanY, color='b')
    ax_px.plot(sumMeanX, sumMeanY, 'yo')
    ax_px.plot(sumMeanX[0], sumMeanY[0], 'ro', markersize=10, linewidth=0.5)
    ax_px.set_title('Global frame alignment')

    plotter.tightLayout()

    return plotter
