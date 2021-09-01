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

import io
from collections import OrderedDict
from emtable import Table
from tifffile import TiffFile


class OpticsGroups:
    """ Store information about optics groups in an indexable way.
    Existing groups can be accessed by number of name.
    """

    def __init__(self, opticsTable):
        self.__fromTable(opticsTable)

    def __fromTable(self, opticsTable):
        self._dict = OrderedDict()
        # Also allow indexing by name
        self._dictName = OrderedDict()
        # Map optics rows both by name and by number
        for og in opticsTable:
            self.__store(og)

    def __store(self, og):
        self._dict[og.rlnOpticsGroup] = og
        self._dictName[og.rlnOpticsGroupName] = og

    def __getitem__(self, item):
        if isinstance(item, int):
            return self._dict[item]
        elif isinstance(item, str):
            return self._dictName[item]
        raise TypeError("Unsupported type '%s' of item '%s'"
                        % (type(item), item))

    def __contains__(self, item):
        return item in self._dict or item in self._dictName

    def __iter__(self):
        """ Iterate over all optics groups. """
        return iter(self._dict.values())

    def __len__(self):
        return len(self._dict)

    def __str__(self):
        return self.toString()

    def first(self):
        """ Return first optics group. """
        return next(iter(self._dict.values()))

    def update(self, ogId, **kwargs):
        og = self.__getitem__(ogId)
        newOg = og._replace(**kwargs)
        self.__store(newOg)
        return newOg

    def updateAll(self, **kwargs):
        """ Update all Optics Groups with these values. """
        missing = {k: v for k, v in kwargs.items() if not self.hasColumn(k)}
        existing = {k: v for k, v in kwargs.items() if self.hasColumn(k)}

        self.addColumns(**missing)

        for og in self:
            self.update(og.rlnOpticsGroup, **existing)

    def add(self, newOg):
        self.__store(newOg)

    def hasColumn(self, colName):
        return hasattr(self.first(), colName)

    def addColumns(self, **kwargs):
        """ Add new columns with default values (type inferred from it). """
        items = self.first()._asdict().items()
        cols = [Table.Column(k, type(v)) for k, v in items]

        for k, v in kwargs.items():
            cols.append(Table.Column(k, type(v)))

        t = Table(columns=cols)

        for og in self._dict.values():
            values = og._asdict()
            values.update(kwargs)
            t.addRow(**values)

        self.__fromTable(t)

    @staticmethod
    def fromStar(starFilePath):
        """ Create an OpticsGroups from a given STAR file.
        """
        return OpticsGroups(Table(fileName=starFilePath, tableName='optics'))

    @staticmethod
    def fromString(stringValue):
        """ Create an OpticsGroups from string content (STAR format)
        """
        f = io.StringIO(stringValue)
        t = Table()
        t.readStar(f, tableName='optics')
        return OpticsGroups(t)

    @staticmethod
    def fromImages(imageSet):
        acq = imageSet.getAcquisition()
        params = {'rlnImageSize': imageSet.getXDim(),
                  'rlnMicrographPixelSize': imageSet.getSamplingRate(),
                  'rlnMicrographOriginalPixelSize': imageSet.getSamplingRate()}
        try:
            og = OpticsGroups.fromString(acq.opticsGroupInfo.get())
            # always update sampling and image size from the set
            og.updateAll(**params)
            return og
        except:
            params.update({
                'rlnVoltage': acq.getVoltage(),
                'rlnSphericalAberration': acq.getSphericalAberration(),
                'rlnAmplitudeContrast': acq.getAmplitudeContrast(),
            })
            return OpticsGroups.create(**params)

    @staticmethod
    def create(**kwargs):
        opticsString1 = """

# version 30001

data_optics

loop_
_rlnOpticsGroupName #1
_rlnOpticsGroup #2
_rlnMicrographOriginalPixelSize #3
_rlnVoltage #4
_rlnSphericalAberration #5
_rlnAmplitudeContrast #6
_rlnImageSize #7
_rlnImageDimensionality #8
opticsGroup1            1      1.000000   300.000000     2.700000     0.100000     256            2
        """

        og = OpticsGroups.fromString(opticsString1)
        fog = og.first()
        newColumns = {k: v for k, v in kwargs.items() if not hasattr(fog, k)}
        og.addColumns(**newColumns)
        og.update(1, **kwargs)
        return og

    def _write(self, f):
        # Create columns from the first row
        items = self.first()._asdict().items()
        cols = [Table.Column(k, type(v)) for k, v in items]
        t = Table(columns=cols)
        for og in self._dict.values():
            t.addRow(*og)
        t.writeStar(f, tableName='optics')

    def toString(self):
        """ Return a string (STAR format) with the current optics groups.
        """
        f = io.StringIO()
        self._write(f)
        result = f.getvalue()
        f.close()

        return result

    def toStar(self, starFile):
        """ Write current optics groups to a given file.
        """
        self._write(starFile)

    def toImages(self, imageSet):
        """ Store the optics groups information in the image acquisition.
        """
        imageSet.getAcquisition().opticsGroupInfo.set(self.toString())


def parseMovieAlignment2(logFile):
    """ Get global frame shifts relative to the first frame
    (for the plots)
    """
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
    import xml.etree.ElementTree as ET
    defects = []  # x y w h

    if not fn.endswith(".gain"):
        return defects

    print("Parsing defects from EER gain file..")

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
