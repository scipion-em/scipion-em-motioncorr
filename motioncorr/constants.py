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
import pyworkflow.utils as pwutils

CUDA_LIB = 'CUDA_LIB'

MOTIONCORR = 'MOTIONCORR'
MOTIONCORR_HOME = 'MOTIONCORR_HOME'
MOTIONCORR_CUDA_LIB = 'MOTIONCORR_CUDA_LIB'
MOTIONCORR_BIN = 'MOTIONCORR_BIN'

MOTIONCOR2 = 'MOTIONCOR2'
MOTIONCOR2_HOME = 'MOTIONCOR2_HOME'
MOTIONCOR2_CUDA_LIB = 'MOTIONCOR2_CUDA_LIB'
MOTIONCOR2_BIN = 'MOTIONCOR2_BIN'


# FIXME: Remove the following block when done
# def _getHome(key, default):
#     """ Get the required home path, if not present..
#     the default value will be used from EM_ROOT.
#     """
#     return os.environ.get(key, os.path.join(os.environ['EM_ROOT'], default))
#
# # For some reason all variables end up with SCIPION HOME prepended %&%#!
# MOTIONCORR = pwutils.getEnvVariable('MOTIONCORR')
# MOTIONCORR = os.path.basename(MOTIONCORR)
#
# MOTIONCOR2_BIN = pwutils.getEnvVariable('MOTIONCOR2_BIN', 'motioncor2')
# MOTIONCOR2_BIN = os.path.basename(MOTIONCOR2_BIN)
#
# MOTIONCORR_PATH = os.path.join(_getHome('MOTIONCORR_HOME', 'motioncorr'),
#                                'bin', MOTIONCORR)
#
# MOTIONCOR2_PATH = os.path.join(_getHome('MOTIONCOR2_HOME', 'motioncor2'),
#                                'bin', MOTIONCOR2_BIN)


