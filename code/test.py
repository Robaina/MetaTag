from phyloplacement.database import filterFASTAByHMM, removeDuplicatesFromFastaByID


if __name__ == '__main__':
    
    # hmm = '/home/robaina/hmmer/narGTIGR01580.1.HMM'
    # # ifile = '/home/robaina/cleangenomes/results/Marref_V6.fasta'
    # ifile = '/home/robaina/Documents/MAR_database/mardb_proteins_V6.faa'

    # filterFASTAByHMM(hmm_model=hmm,
    #                  input_fasta=ifile,
    #                  output_fasta=None,
    #                  method='hmmsearch')

    # input_fasta = '/home/robaina/Documents/MAR_database/mardb_proteins_V6.faa'
    
    hmm = '../data/hmms/narGTIGR01580.1.HMM'
    input_fasta = '/usr/gonzalez/metagenomes/MarPeptides/mardb_proteins_V6.fasta'
    input_fasta_no_dup = '/usr/gonzalez/metagenomes/MarPeptides/mardb_proteins_V6_no_duplicates.fasta'
    
    removeDuplicatesFromFastaByID(input_fasta, output_fasta=input_fasta_no_dup) 

    filterFASTAByHMM(hmm_model=hmm,
                     input_fasta=input_fasta_no_dup,
                     output_fasta=None,
                     method='hmmsearch')
