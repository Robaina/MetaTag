#!/usr/bin/env python
# conda activate traits

from phyloplacement.database import filterFASTAByHMM, removeDuplicatesFromFasta
from phyloplacement.phylotree import (runMuscle, runTrimal,
                                      convertFastaAlnToPhylip, runFastTree,
                                      runIqTree)

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

    # convertFastaAlnToPhylip(input_fasta_aln='data/nxr/mardb_proteins_V6_TIGR015180.1.fasta.aln',
    #                         output_file='data/nxr/mardb_proteins_V6_TIGR015180.1.phylip')
    
    # runFastTree(input_algns='data/nxr/mardb_proteins_V6_TIGR015180.1.fasta_trimal.aln')

    runIqTree(input_algns='data/nxr/mardb_proteins_V6_TIGR015180.1.fasta_trimal.aln')