from galaxy.datatypes import data
from galaxy.datatypes.metadata import MetadataElement

import os
import logging

log = logging.getLogger(__name__)

class HDR(data.Text):
    file_ext = "hdr"
