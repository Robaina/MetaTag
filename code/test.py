#!/usr/bin/env python
# conda activate traits

import re
import pyfastx
from phyloplacement.utils import readFromPickleFile, setDefaultOutputPath 
from phyloplacement.mardb import getMarDBentryCode, getMARdbGenomeByEntryCode


if __name__ == '__main__':

    # Retrieve mardb genomes
    label_dict = readFromPickleFile('/home/robaina/Documents/TRAITS/tests/ref_reduced_clean_id_dict.pickle')
    nxr_entry_codes = {getMarDBentryCode(v) for v in label_dict.values()}

    for entry_code in nxr_entry_codes:
        print(f'Running entry code: {entry_code}')
        getMARdbGenomeByEntryCode(
            input_fasta='/home/robaina/Documents/MAR_database/mardb_assembly_V6.fa',
            entry_code=entry_code,
            output_fasta=f'/home/robaina/Documents/TRAITS/tests/nxr_genomes/{entry_code}.fa'
        )
    
    print('Done!')