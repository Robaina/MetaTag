from phyloplacement.database import filterFASTAByHMM


if __name__ == '__main__':
    
    hmm = '/home/robaina/hmmer/narGTIGR01580.1.HMM'
    ifile = '/home/robaina/cleangenomes/results/Marref_V6.fasta'

    filterFASTAByHMM(hmm_model=hmm,
                     input_fasta=ifile,
                     output_fasta=None,
                     method='hmmsearch')
