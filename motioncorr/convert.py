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

from tifffile import TiffFile
import xml.etree.ElementTree as ET
import logging
logger = logging.getLogger(__name__)


def parseMovieAlignment2(logFile):
    """ Get global frame shifts relative to the first frame. """
    first = None
    xshifts = []
    yshifts = []
    with open(logFile, 'r') as f:
        for line in f:
            l = line.strip()
            if '#' not in l and len(l) > 0:
                parts = l.split()
                if first is None:  # read the first frame number
                    first = int(parts[0])  # take the id from first column
                # take the shifts from the last two columns of the line
                xshifts.append(float(parts[1]))
                yshifts.append(float(parts[2]))

    xoff, yoff = -xshifts[0], -yshifts[0]
    xshifts = [x + xoff for x in xshifts]
    yshifts = [y + yoff for y in yshifts]

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


def parseEERDefects(fn):
    """ Extract defects coords from XML metadata inside EER *.gain file. """
    defects = []  # x y w h

    if not fn.endswith(".gain"):
        return defects

    logger.info("Parsing defects from EER gain file..")

    with TiffFile(fn) as tif:
        for page in tif.pages:
            for tag in page.tags:
                if tag.code == 65100:  # TFS EER gain Metadata
                    gainStr = tag.value
                    break
            break

    for item in ET.fromstring(gainStr.decode('utf-8')):
        if item.tag == "point":
            point = item.text.split(",")
            defects.append((point[0], point[1], 1, 1))
        elif item.tag == "area":
            area = item.text.split(",")
            defects.append((area[0],
                            area[1],
                            int(area[2])-int(area[0])+1,
                            int(area[3])-int(area[1])+1))
        elif item.tag == "col":
            area = item.text.split("-")
            defects.append((area[0],
                            0,
                            int(area[1])-int(area[0])+1,
                            4096))
        elif item.tag == "row":
            area = item.text.split("-")
            defects.append((0, area[0], 0, 4096,
                            int(area[1])-int(area[0])+1))

    return defects
