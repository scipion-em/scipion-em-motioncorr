# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# **************************************************************************

import os

import pwem
import pyworkflow.utils as pwutils

from .constants import *


__version__ = '3.1.0'
_references = ['Zheng2017']


class Plugin(pwem.Plugin):
    _homeVar = MOTIONCOR2_HOME
    _pathVars = [MOTIONCOR2_HOME]
    _supportedVersions = ['1.2.6', '1.3.0', '1.3.1', '1.3.2', '1.4.0']
    _url = "https://github.com/scipion-em/scipion-em-motioncorr"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(MOTIONCOR2_HOME, 'motioncor2-1.4.0')
        cls._defineVar(MOTIONCOR2_BIN, 'MotionCor2_1.4.0_Cuda101')
        cls._defineVar(MOTIONCOR2_CUDA_LIB, pwem.Config.CUDA_LIB)

    @classmethod
    def getProgram(cls):
        return os.path.join(cls.getHome('bin'),
                            os.path.basename(cls.getVar(MOTIONCOR2_BIN)))

    @classmethod
    def getEnviron(cls):
        """ Return the environment to run motioncor2. """
        environ = pwutils.Environ(os.environ)
        # Get motioncor2 CUDA library path if defined
        cudaLib = cls.getVar(MOTIONCOR2_CUDA_LIB, pwem.Config.CUDA_LIB)
        environ.addLibrary(cudaLib)

        return environ

    @classmethod
    def defineBinaries(cls, env):
        env.addPackage('motioncor2', version='1.2.6',
                       tar='motioncor2-1.2.6.tgz')

        env.addPackage('motioncor2', version='1.3.0',
                       tar='motioncor2-1.3.0.tgz')

        env.addPackage('motioncor2', version='1.3.1',
                       tar='motioncor2-1.3.1.tgz')

        env.addPackage('motioncor2', version='1.3.2',
                       tar='motioncor2-1.3.2.tgz')

        env.addPackage('motioncor2', version='1.4.0',
                       tar='motioncor2-1.4.0.tgz',
                       default=True)
