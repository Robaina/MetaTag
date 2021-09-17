from phyloplacement.database import filterFASTAByHMM, removeDuplicatesFromFasta
from phyloplacement.phylotree import (runMuscle, runTrimal,
                                      convertFastaAlnToPhylip, runFastTree)

if __name__ == '__main__':
        
    # hmm = 'data/hmms/narGTIGR01580.1.HMM'
    # input_fasta = '/home/robaina/Documents/MAR_database/mardb_proteins_V6.faa'
    # input_fasta = '/usr/gonzalez/metagenomes/MarPeptides/mardb_proteins_V6.fasta'
    # input_fasta_no_dup = '/usr/gonzalez/metagenomes/MarPeptides/mardb_proteins_V6_no_duplicates.fasta'
    
    # removeDuplicatesFromFasta(input_fasta, output_fasta=input_fasta_no_dup)

    # filterFASTAByHMM(hmm_model=hmm,
    #                  input_fasta=input_fasta_no_dup,
    #                  output_fasta=None,
    #                  method='hmmsearch')
  
    # runMuscle(input_fasta='data/nxr/mardb_proteins_V6_TIGR015180.1.fasta',
    #           output_file='data/nxr/mardb_proteins_V6_TIGR015180.1.fasta.aln',
    #           maxiters=None)

    # runTrimal(input_aln='data/nxr/mardb_proteins_V6_TIGR015180.1.fasta.aln',
    #           output_aln=None)

    convertFastaAlnToPhylip(input_fasta_aln='data/nxr/mardb_proteins_V6_TIGR015180.1.fasta.aln',
                            output_file='data/nxr/mardb_proteins_V6_TIGR015180.1.phylip')
    
    runFastTree(input_phylip='data/nxr/mardb_proteins_V6_TIGR015180.1.phylip')



"""
NOTE: producing this error during phylip conversion

Traceback (most recent call last):
  File "/home/robaina/miniconda3/envs/traits/lib/python3.9/site-packages/Bio/File.py", line 72, in as_handle
    with open(handleish, mode, **kwargs) as fp:
TypeError: expected str, bytes or os.PathLike object, not TextIOWrapper

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/home/robaina/Documents/TRAITS/code/accesory/fasta_to_phylip.py", line 34, in <module>
    main()
  File "/home/robaina/Documents/TRAITS/code/accesory/fasta_to_phylip.py", line 30, in main
    AlignIO.write(records, output_handle, "phylip")
  File "/home/robaina/miniconda3/envs/traits/lib/python3.9/site-packages/Bio/AlignIO/__init__.py", line 215, in write
    count = writer_class(fp).write_file(alignments)
  File "/home/robaina/miniconda3/envs/traits/lib/python3.9/site-packages/Bio/AlignIO/Interfaces.py", line 129, in write_file
    self.write_alignment(alignment)
  File "/home/robaina/miniconda3/envs/traits/lib/python3.9/site-packages/Bio/AlignIO/PhylipIO.py", line 101, in write_alignment
    raise ValueError(
ValueError: Repeated name 'GCA_003230' (originally 'GCA_003230455.1_01774_MMP09240198'), possibly due to truncation
"""