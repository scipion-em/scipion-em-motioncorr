# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# **************************************************************************

import os

import pyworkflow.em as pwem
import pyworkflow.utils as pwutils

from .constants import *


_references = ['Li2013', 'Zheng2017']


class Plugin(pwem.Plugin):
    _homeVar = MOTIONCOR2_HOME
    _pathVars = [MOTIONCOR2_HOME]
    _versions = {
        MOTIONCOR2: ['01302017', '1.0.2', '1.0.5', '1.1.0'],
        MOTIONCORR: ['2.1']}

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(MOTIONCORR_HOME, 'motioncorr-2.1')
        cls._defineVar(MOTIONCORR_BIN, 'dosefgpu_driftcorr_7.5')
        cls._defineEmVar(MOTIONCOR2_HOME, 'motioncor2-1.1.0')
        cls._defineVar(MOTIONCOR2_BIN, 'MotionCor2_1.1.0-Cuda80')

    @classmethod
    def getProgram(cls, mcVar=MOTIONCOR2):
        """ Re-implement this method because this plugin really deals
        with two different programs that can be the "home".
        """
        return os.path.join(cls.getVar('%s_HOME' % mcVar), 'bin',
                            os.path.basename(cls.getVar('%s_BIN' % mcVar)))

    @classmethod
    def getEnviron(cls, mcVar=MOTIONCOR2):
        """ Return the environment to run motioncor2 or motioncorr. """
        environ = pwutils.Environ(os.environ)
        cudaLib = cls.getCudaLib(mcVar, environ)
        environ.addLibrary(cudaLib)

        return environ

    @classmethod
    def getActiveVersion(cls, mcVar=MOTIONCOR2):
        return pwem.Plugin.getActiveVersion(home=cls.getProgram(mcVar),
                                            versions=cls._versions[mcVar])

    @classmethod
    def getSupportedVersions(cls, mcVar=MOTIONCOR2):
        return cls._versions[mcVar]

    @classmethod
    def getCudaLib(cls, mcVar=MOTIONCOR2, environ=None):
        e = environ or pwutils.Environ(os.environ)
        return e.getFirst(['%s_CUDA_LIB' % mcVar, CUDA_LIB])

    @classmethod
    def defineBinaries(cls, env):
        env.addPackage('motioncorr', version='2.1',
                       tar='motioncorr_v2.1.tgz',
                       default=True)

        env.addPackage('motioncor2', version='17.01.30',
                       tar='motioncor2_01302017.tgz')

        env.addPackage('motioncor2', version='1.0.2',
                       tar='motioncor2-1.0.2.tgz')

        env.addPackage('motioncor2', version='1.0.5',
                       tar='motioncor2-1.0.5.tgz')

        env.addPackage('motioncor2', version='1.1.0',
                       tar='motioncor2-1.1.0.tgz',
                       default=True)


pwem.Domain.registerPlugin(__name__)
