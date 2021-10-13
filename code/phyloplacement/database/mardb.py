"""
Tools to process MARdb data
"""

import re
import pyfastx
from phyloplacement.utils import setDefaultOutputPath, terminalExecute 


db_entry = re.compile('\[mmp_id=(.*)\] ')

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
                              output_fasta: str = None) -> None:
    """
    Get full or partial genomes with given MARdb entry code.
    If more than one record found for entry code (e.g., when dealing
    with contigs of partial genomes), then all records are merged together
    into a single fasta file.
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta,
                                            tag=f'_genome_{entry_code}',
                                            extension='.fa')
    cmd_str = (
        f'grep -A1 {entry_code} {input_fasta} | grep -v {entry_code} '
        '''| tr -d ''' + r"'\n'" + " "
        f'> {output_fasta}'
    )
    terminalExecute(cmd_str, suppress_output=False)
    with open(output_fasta, 'r+') as file:
        content = file.read()
        file.seek(0)
        file.write(f'>{entry_code}\n' + content)

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