from galaxy.datatypes import data
from galaxy.datatypes.metadata import MetadataElement

import os
import logging

log = logging.getLogger(__name__)

class BIL(data.Text):
    file_ext = ".bil"


