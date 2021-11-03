"""
Tools to process Paoli et al. 2019 data
https://doi.org/10.1101/2021.03.24.436479
"""

import os
import re
import shutil
import tarfile
import pyfastx
from phyloplacement.utils import (listTarDir, setDefaultOutputPath,
                                  terminalExecute,
                                  createTemporaryFilePath) 


def getGenomeIDs(paoli_tar: str) -> list:
    """
    Get Paoli database genome IDs from tar file
    """
    return [
        file.split('aPaoli/')[1].split('.fasta')[0]
        for file in listTarDir(paoli_tar)[1:]
        ]