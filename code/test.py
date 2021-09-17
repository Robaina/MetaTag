from phyloplacement.database import filterFASTAByHMM, removeDuplicatesFromFasta


if __name__ == '__main__':
        
    hmm = 'data/hmms/narGTIGR01580.1.HMM'
    # input_fasta = '/home/robaina/Documents/MAR_database/mardb_proteins_V6.faa'
    input_fasta = '/usr/gonzalez/metagenomes/MarPeptides/mardb_proteins_V6.fasta'
    input_fasta_no_dup = '/usr/gonzalez/metagenomes/MarPeptides/mardb_proteins_V6_no_duplicates.fasta'
    
    # removeDuplicatesFromFastaByID(input_fasta, output_fasta=input_fasta_no_dup) 
    # removeDuplicatesFromFasta(input_fasta, output_fasta=input_fasta_no_dup)

    filterFASTAByHMM(hmm_model=hmm,
                     input_fasta=input_fasta_no_dup,
                     output_fasta=None,
                     method='hmmsearch')

# ValueError: The ID or alternative IDs of Hit 'WP_113256187.1_MMP09508687' exists in this QueryResult