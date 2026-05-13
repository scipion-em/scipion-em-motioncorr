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
from os.path import dirname

import pwem
import pyworkflow.utils as pwutils
from pyworkflow import SPA, TOMO

from .constants import *


__version__ = '3.17.2'
_references = ['Zheng2017']


class Plugin(pwem.Plugin):
    _homeVar = MOTIONCOR_HOME
    _pathVars = [MOTIONCOR_CUDA_LIB]
    _supportedVersions = [V1_2_4]
    _url = "https://github.com/scipion-em/scipion-em-motioncorr"
    _processingField = [SPA, TOMO]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(MOTIONCOR_HOME, f'motioncor3-{V1_2_4}')
        cls._defineVar(MOTIONCOR_CUDA_LIB, pwem.Config.CUDA_LIB)
        # Define the variable default value based on the guessed cuda version
        cls._defineVar(MOTIONCOR_BIN, 'MotionCor3')

    @classmethod
    def getProgram(cls):
        return os.path.join(cls.getHome('bin'),
                            os.path.basename(cls.getVar(MOTIONCOR_BIN)))

    @classmethod
    def validateInstallation(cls):
        """ Check if the binaries are properly installed. """
        try:
            if not os.path.exists(cls.getProgram()):
                return [f"{cls.getProgram()} does not exist, please verify "
                        "the following variables or edit the config file:\n\n"
                        f"{MOTIONCOR_HOME}: {cls.getVar(MOTIONCOR_HOME)}\n"
                        f"{MOTIONCOR_BIN}: {cls.getVar(MOTIONCOR_BIN)}"]

            if not os.path.exists(cls.getVar(MOTIONCOR_CUDA_LIB)):
                return [f"{cls.getVar(MOTIONCOR_CUDA_LIB)} does not exist, "
                        f"please verify {MOTIONCOR_CUDA_LIB} variable."]

        except Exception as e:
            return [f"validateInstallation fails: {str(e)}"]

    @classmethod
    def getEnviron(cls):
        """ Return the environment to run motioncor. """
        environ = pwutils.Environ(os.environ)
        # Get motioncor CUDA library path if defined
        cudaLib = cls.getVar(MOTIONCOR_CUDA_LIB, pwem.Config.CUDA_LIB)
        environ.addLibrary(cudaLib)

        return environ

    @classmethod
    def defineBinaries(cls, env):

        for v in cls._supportedVersions:
            if v == V1_2_4:
                commit = "dd8b683"
            else:
                commit = "main"

        binary = cls.getVar(MOTIONCOR_BIN)
        MOTIONCOR_INSTALLED = f'{MOTIONCOR_BIN}_installed'

        cmd = [
            f'cd .. && rm -rf motioncor3-{v} && '
            f'git clone https://github.com/CZImagingInstitute/MotionCor3.git motioncor3-{v} && '
            f'cd motioncor3-{v} && '
            f'git checkout {commit} && '
            'mkdir bin && '
            r"sed -i '/^CUFLAG = -Xptxas -dlcm=ca -O2 \\/,/code=sm_70$/c\CUFLAG = -Xptxas -dlcm=ca -O2 -arch=all' makefile11 && "       
            r"sed -i '/-L\/usr\/lib64 \\/a\    -Xcompiler -no-pie \\' makefile11 && "
            r"sed -i '/^PRJLIB =/a BINARY ?= MotionCor3' makefile11 && "
            r"sed -i 's/-o MotionCor3/-o $(BINARY)/' makefile11 && "
            f' make exe -f makefile11 BINARY={binary} CUDAHOME={dirname(pwem.Config.CUDA_LIB)} -j 16 && '
            f'cp {binary} bin/ && '
            f'touch {MOTIONCOR_INSTALLED}'
            ]
        
        installationCmds = [
            (cmd, MOTIONCOR_INSTALLED)
        ]

        env.addPackage('motioncor3', version=v,
                        tar = 'void.tgz',
                        neededProgs = ['git', 'gcc', 'g++', 'make', 'cmake'],
                        commands = installationCmds,
                        updateCuda = True,
                        default = True)
