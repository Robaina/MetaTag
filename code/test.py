#!/usr/bin/env python
# conda activate traits

import re
import pyfastx
from phyloplacement.utils import readFromPickleFile, setDefaultOutputPath 


db_entry = re.compile('\[mmp_id=(.*)\] ')

def getMarDBentryCode(label: str) -> str:
    return re.search(db_entry, label).group(1)


def getMARDBsequencesByIDs(input_fasta: str, entry_codes: set,
                        output_fasta: str = None) -> None:
    """
    Filter records in mardb fasta file matching provided entry codes
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, '_fitered')
    
    fasta = pyfastx.Fasta(input_fasta, build_index=False, full_name=True)
    with open(output_fasta, 'w') as outfile:
        for record_name, record_seq in fasta:
            entry_code = getMarDBentryCode(record_name)
            if entry_code in entry_codes:
                outfile.write(f'>{record_name}\n{record_seq}\n')

if __name__ == '__main__':


    # Retrieve mardb nucleotide sequences
    label_dict = readFromPickleFile('/home/robaina/Documents/TRAITS/tests/ref_reduced_clean_id_dict.pickle')
    nxr_entry_codes = {getMarDBentryCode(v) for v in label_dict.values()}

    getMARDBsequencesByIDs(
        input_fasta='/home/robaina/Documents/MAR_database/mardb_nucleotides_V6.fna',
        entry_codes=nxr_entry_codes,
        output_fasta='/home/robaina/Documents/TRAITS/nxr_nt.fasta'
    )
    
    print('Done!')


    """
    I see after running the script that entry ids are not unique, there are several nt sequences with
    the same id. How can one then get nt sequences for given mardb protein sequence?
    """