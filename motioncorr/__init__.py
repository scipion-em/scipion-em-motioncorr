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

import pyworkflow.em
from pyworkflow.utils import Environ, getEnvVariable

from .constants import *


_references = ['Li2013', 'Zheng2017']

# The following class is required for Scipion to detect this Python module
# as a Scipion Plugin. It needs to specify the PluginMeta __metaclass__
# Some function related to the underlying package binaries need to be
# implemented
class Plugin:
    __metaclass__ = pyworkflow.em.PluginMeta

    @classmethod
    def getEnviron(cls, useMC2=False):
        """ Return the environ settings to run motioncorr programs. """
        environ = Environ(os.environ)

        if os.path.exists(MOTIONCORR_PATH) and not useMC2:
            environ.update({'PATH': os.path.join(os.environ['MOTIONCORR_HOME'], 'bin')},
                           position=Environ.BEGIN)

        if os.path.exists(MOTIONCOR2_PATH):
            environ.update({'PATH': os.path.join(os.environ['MOTIONCOR2_HOME'], 'bin')},
                           position=Environ.BEGIN)

        cudaLib = cls.getCudaLib(environ, useMC2)
        environ.addLibrary(cudaLib)

        return environ

    @classmethod
    def getActiveVersion(cls, var):
        varHome = var + '_HOME'
        path = os.environ[varHome]
        for v in cls.getSupportedVersions(var):
            if v in path or v in os.path.realpath(path):
                return v
        return ''

    @classmethod
    def getSupportedVersions(cls, var):
        if var == 'MOTIONCORR':
            return ['2.1']
        elif var == 'MOTIONCOR2':
            return ['01302017', '1.0.0', '1.0.2', '1.0.4', '1.0.5']

    @classmethod
    def getCudaLib(cls, environ=None, useMC2=False):

        e = environ or Environ(os.environ)
        cudaLib = MOTIONCOR2_CUDA_LIB if useMC2 else MOTIONCORR_CUDA_LIB
        return e.getFirst(cudaLib, CUDA_LIB)

    @classmethod
    def validateInstallation(cls):
        """ This function will be used to check if package is
        properly installed."""
        environ = cls.getEnviron()

        missingPaths = ["%s: %s" % (var, environ[var])
                        for var in [GAUTOMATCH_HOME]
                        if not os.path.os.path.exists(environ[var])]

        return (["Missing variables:"] + missingPaths) if missingPaths else []

