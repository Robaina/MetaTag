"""
Tools to process MARdb data
"""

import re
import os
import shutil
import tempfile
import pyfastx
from phyloplacement.utils import setDefaultOutputPath, terminalExecute 

db_entry = re.compile('\[mmp_id=(.*)\] ')
only_letters = re.compile('[^A-Z]')

def getMarDBentryCode(label: str) -> str:
    return re.search(db_entry, label).group(1)

def filterMarDBrecordsbyEntryCodes(input_fasta: str, entry_codes: set,
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

def getMARdbGenomeByEntryCode(entry_code: str, input_fasta: str,
                              output_fasta: str = None,
                              clean: bool = True) -> None:
    """
    Get full or partial genomes with given MARdb entry code.
    If clean = True, remove characters which are not letters
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta,
                                            tag=f'_genome_{entry_code}',
                                            extension='.fa')
    
    def is_empty(fasta_file):
        with open(fasta_file, 'r') as file:
            return '>' not in file.read()

    def cleanOutputFasta(output_fasta: str) -> None:
        with tempfile.TemporaryFile() as tfile, \
             open(output_fasta, 'r') as file:
            for line in file:
                if '>' in line:
                    line = only_letters.sub('', line)
                tfile.write(line)
        shutil.move(tfile, output_fasta)

    cmd_str = (
        f'grep -A1 {entry_code} {input_fasta} > {output_fasta}'
    )
    terminalExecute(cmd_str, suppress_output=False)
    if clean:
        cleanOutputFasta(output_fasta)
    if is_empty(output_fasta):
        os.remove(output_fasta)

def relabelMarDB(label_dict: dict) -> dict:
    """
    Convert mardb long labels into short labels 
    displaying mardb id and species (if present)
    """
    db_code_pattern = re.compile('\[mmp_(.*)\]')
    species_pattern = re.compile('\[(.*?)\]')

    def editMarDBlabel(label: str) -> str:
        try:
            species = re.search(
                species_pattern,
                re.sub(db_code_pattern, '', label)
                ).group(1)
        except:
            species = 'Undetermined'
        mar_id = label.split(' ')[0]
        return f'{mar_id}_{species}'

    return {
        k: editMarDBlabel(v)
        for k, v in label_dict.items()
    }

# def getMARdbGenomeByEntryCode(entry_code: str, input_fasta: str,
#                               output_fasta: str = None) -> None:
#     """
#     Get full or partial genomes with given MARdb entry code.
#     If more than one record found for entry code (e.g., when dealing
#     with contigs of partial genomes), then all records are merged together
#     into a single fasta file.
#     """
#     if output_fasta is None:
#         output_fasta = setDefaultOutputPath(input_fasta,
#                                             tag=f'_genome_{entry_code}',
#                                             extension='.fa')
#     cmd_str = (
#         f'grep -A1 {entry_code} {input_fasta} | grep -v {entry_code} '
#         '''| tr -d ''' + r"'\n'" + " "
#         f'> {output_fasta}'
#     )
#     terminalExecute(cmd_str, suppress_output=False)
#     with open(output_fasta, 'r+') as file:
#         content = file.read()
#         file.seek(0)
#         file.write(f'>{entry_code}\n' + content)