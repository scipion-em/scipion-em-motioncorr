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


__version__ = '3.9'
_references = ['Zheng2017']


class Plugin(pwem.Plugin):
    _homeVar = MOTIONCOR2_HOME
    _pathVars = [MOTIONCOR2_CUDA_LIB]
    _supportedVersions = ['1.4.0', '1.4.2', '1.4.4', '1.4.5', '1.4.7', '1.5.0']
    _url = "https://github.com/scipion-em/scipion-em-motioncorr"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(MOTIONCOR2_HOME, 'motioncor2-1.5.0')
        cls._defineVar(MOTIONCOR2_CUDA_LIB, pwem.Config.CUDA_LIB)

        # Define the variable default value based on the guessed cuda version
        cudaVersion = cls.guessCudaVersion(MOTIONCOR2_CUDA_LIB)
        cls._defineVar(MOTIONCOR2_BIN, 'MotionCor2_1.5.0_Cuda%s%s_05-31-2022' % (
            cudaVersion.major, cudaVersion.minor))


    @classmethod
    def getProgram(cls):
        return os.path.join(cls.getHome('bin'),
                            os.path.basename(cls.getVar(MOTIONCOR2_BIN)))

    @classmethod
    def validateInstallation(cls):
        """ Check if the binaries are properly installed. """
        try:
            if not os.path.exists(cls.getProgram()):
                return [f"{cls.getProgram()} does not exist, please verify "
                        "the following variables or edit the config file:\n\n"
                        f"{MOTIONCOR2_HOME}: {cls.getVar(MOTIONCOR2_HOME)}\n"
                        f"{MOTIONCOR2_BIN}: {cls.getVar(MOTIONCOR2_BIN)}"]

            if not os.path.exists(cls.getVar(MOTIONCOR2_CUDA_LIB)):
                return [f"{cls.getVar(MOTIONCOR2_CUDA_LIB)} does not exist, "
                        f"please verify {MOTIONCOR2_CUDA_LIB} variable."]

        except Exception as e:
            return ["validateInstallation fails: %s" % e]

    @classmethod
    def versionGE(cls, version):
        """ Return True if current version of motioncor2 is greater
         or equal than the input argument.
         Params:
            version: string version (semantic version, e.g 1.0.1)
        """
        v1 = int(Plugin.getActiveVersion().replace('.', ''))
        v2 = int(version.replace('.', ''))

        if v1 < v2:
            return False
        return True

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
        for v in cls._supportedVersions:
            env.addPackage('motioncor2', version=v,
                           tar='motioncor2-%s.tgz' % v,
                           default=v == '1.5.0')
