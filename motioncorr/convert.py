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

from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md


def parseMovieAlignment2(logFile):
    """ Get global frame shifts relative to the first frame
    (for the plots)
    """
    f = open(logFile, 'r')
    first = None
    xshifts = []
    yshifts = []

    for line in f:
        l = line.strip()
        if '#' not in l and len(l) > 0:
            parts = l.split()
            if first is None:  # read the first frame number
                first = int(parts[0])  # take the id from first column
            # take the shifts from the last two columns of the line
            xshifts.append(float(parts[1]))
            yshifts.append(float(parts[2]))
    f.close()
    return xshifts, yshifts


def getMovieFileName(movie):
    """ Add the :mrcs or :ems extensions to movie files to be
    recognized by Xmipp as proper stack files.
    Note: Copied from xmipp3.convert
    """
    fn = movie.getFileName()
    if fn.endswith('.mrc'):
        fn += ':mrcs'
    elif fn.endswith('.em'):
        fn += ':ems'

    return fn


def writeShiftsMovieAlignment(movie, xmdFn, s0, sN):
    """Note: Function copied from xmipp3.convert to avoid dependency.
    """
    movieAlignment = movie.getAlignment()
    shiftListX, shiftListY = movieAlignment.getShifts()
    # Generating metadata for global shifts
    a0, aN = movieAlignment.getRange()
    globalShiftsMD = md.MetaData()
    alFrame = a0

    if s0 < a0:
        for i in range(s0, a0):
            objId = globalShiftsMD.addObject()
            imgFn = ImageHandler.locationToXmipp((i, getMovieFileName(movie)))
            globalShiftsMD.setValue(md.MDL_IMAGE, imgFn, objId)
            globalShiftsMD.setValue(md.MDL_SHIFT_X, 0.0, objId)
            globalShiftsMD.setValue(md.MDL_SHIFT_Y, 0.0, objId)

    for shiftX, shiftY in zip(shiftListX, shiftListY):
        if s0 <= alFrame <= sN:
            objId = globalShiftsMD.addObject()

            imgFn = ImageHandler.locationToXmipp((alFrame, getMovieFileName(movie)))
            globalShiftsMD.setValue(md.MDL_IMAGE, imgFn, objId)
            globalShiftsMD.setValue(md.MDL_SHIFT_X, shiftX, objId)
            globalShiftsMD.setValue(md.MDL_SHIFT_Y, shiftY, objId)
        alFrame += 1

    if sN > aN:
        for j in range(aN, sN):
            objId = globalShiftsMD.addObject()
            imgFn = ImageHandler.locationToXmipp((j+1, getMovieFileName(movie)))
            globalShiftsMD.setValue(md.MDL_IMAGE, imgFn, objId)
            globalShiftsMD.setValue(md.MDL_SHIFT_X, 0.0, objId)
            globalShiftsMD.setValue(md.MDL_SHIFT_Y, 0.0, objId)

    globalShiftsMD.write(xmdFn)
