from typing import Union
import os
import glob
from lxml import etree


def __init__(self,
             path: Union[str, None] = None,
             prefix: Union[str, None] = None,
             version: Union[float, None] = None):
    # ------------------------------------------------------------
    # check and set the value of path
    if path is None:
        # default to using $ESPRESSO_TMPDIR
        if 'ESPRESSO_TMPDIR' in os.environ:
            path = os.environ['ESPRESSO_TMPDIR']
        else:
            path = '.'
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"The provided path does not exist: {path}")
    self.path = path
    # ------------------------------------------------------------
    # check and set the value of the prefix
    if prefix is None:
        savefolders = glob.glob(os.path.join(path, '*.save'))
        if len(savefolders) > 1:
            raise ValueError(
                "More than one save folder under path, please specify a prefix")
        elif len(savefolders) == 0:
            raise ValueError("No save folder detected under path")
        savefolder = savefolders[0]
        prefix = savefolder[:-len('.save')]
    else:
        savefolder = os.path.join(path, f'{prefix}.save')
        if not os.path.exists(self.savefolder):
            raise FileNotFoundError(
                f"Unable to locate save folder: {savefolder}")
    self.savefolder = savefolder
    self.prefix = prefix
    # ------------------------------------------------------------
    # find the data xml file and determine version number
    xml_found = False
    xmlfile_enums = enumerate(
        (os.path.join(self.savefolder, 'data-file-schema.xml'),
         os.path.join(self.savefolder, 'data-file.xml'))
    )
    for i, xmlfile in xmlfile_enums:
        if os.path.exists(xmlfile):
            xml_found = True
            break
    if not xml_found:
        raise ValueError("Unable to locate xml file")
    self.root = etree.parse(xmlfile).getroot()
    if i == 0:
        child = self.root.find('general_info')
        version = child.find('creator').attrib['VERSION']
    elif i == 1:
        child = self.root.find('HEADER')
        version = child.find('CREATOR').attrib['VERSION']
    if version is None:
        raise ValueError("Unable to read version from xmlfile")
    self.version = version
    self.xmlfile = xmlfile
