#!/usr/bin/env python
# conda activate traits

import re
import pyfastx
from phyloplacement.utils import readFromPickleFile, setDefaultOutputPath 
from phyloplacement.mardb import getMarDBentryCode, filterMarDBbyEntryCodes


if __name__ == '__main__':

    # Retrieve mardb genomes
    label_dict = readFromPickleFile('/home/robaina/Documents/TRAITS/tests/ref_reduced_clean_id_dict.pickle')
    nxr_entry_codes = {getMarDBentryCode(v) for v in label_dict.values()}

    filterMarDBbyEntryCodes(
        input_fasta='/home/robaina/Documents/MAR_database/mardb_assembly_V6.fa',
        entry_codes=nxr_entry_codes,
        output_fasta='/home/robaina/Documents/TRAITS/nxr_genomes.fasta'
    )
    
    print('Done!')