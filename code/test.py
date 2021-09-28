#!/usr/bin/env python
# conda activate traits

from phyloplacement.database import (filterFASTAByHMM, removeDuplicatesFromFasta,
                                     filterFastaBySequenceLength, runCDHIT)

from phyloplacement.alignment import (runMuscle, runMAFFT, runTrimal)

from phyloplacement.phylotree import (runFastTree, runEPAng,
                                      runIqTree, runTreeShrink)

from phyloplacement.visualization import plotTreeInBrowser

if __name__ == '__main__':
        
    # hmm = 'data/hmms/narGTIGR01580.1.HMM'
    # input_fasta = '/home/robaina/Documents/MAR_database/mardb_proteins_V6.faa'
    # input_fasta_no_dup = '/home/robaina/Documents/MAR_database/mardb_proteins_V6_no_duplicates.fasta'
    
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

    # runIqTree(input_algns='data/nxr/mardb_proteins_V6_TIGR015180.1.fasta_trimal.aln')
    
    """
    Large (36800) database with Molybdopterin.hmm. 
    Reduced to 17300 with cd-hit (default params)

    Very small (1268) database with TIGR015180 for narG
    """

    # filterFastaBySequenceLength(input_fasta=input_fasta_no_dup, minLength=100)
    # fa = pyfastx.Fasta('/home/robaina/Documents/TRAITS/data/nxr/mardb_proteins_V6_TIGR015180.1.fasta')
    # ids = fa.keys()
    # ids.filter
     
    # runMuscle(input_fasta='data/nxr/kitzinger2021/Nxr_kitzinger_2021.fasta')


    # runTrimal(input_aln='data/nxr/kitzinger2021/Nxr_kitzinger_2021.fasta.aln',
    #           output_aln='data/nxr/kitzinger2021/Nxr_kitzinger_2021.fasta.aln')

    # runIqTree(input_algns='data/nxr/kitzinger2021/Nxr_kitzinger_2021.fasta.aln',
    #           keep_recovery_files=False,
    #           output_dir='data/nxr/kitzinger2021',
    #           output_prefix=None)

    # ******************************************************************************************

    # runTreeShrink(input_tree='/home/robaina/Documents/TRAITS/data/nxr/iqtree_shrink/tree/mardb_proteins_V6_TIGR015180.1.fasta_trimal.aln.treefile',
    #               input_aln='/home/robaina/Documents/TRAITS/data/nxr/iqtree_shrink/tree/mardb_proteins_V6_TIGR015180.1.fasta_trimal.aln',
    #               output_dir='/home/robaina/Documents/TRAITS/data/nxr/iqtree_shrink_output',
    #               output_deleted_nodes=True,
    #               additional_args='--force --centroid -q 0.05')

    # plotTreeInBrowser(input_tree='/home/robaina/Documents/TRAITS/data/nxr/iqtree_shrink_output/mardb_proteins_V6_TIGR015180.1.fasta_trimal.aln_shrink.treefile',
    #                   output_dir='/home/robaina/Documents/TRAITS/data/nxr/tree-viz')
    
    
    runEPAng(input_aln='/home/robaina/Documents/TRAITS/data/nxr/iqtree_shrink_output/mardb_proteins_V6_TIGR015180.1.fasta_trimal_shrink.aln',
             input_tree='/home/robaina/Documents/TRAITS/data/nxr/iqtree_shrink_output/mardb_proteins_V6_TIGR015180.1.fasta_trimal.aln_shrink.treefile',
             input_query='/home/robaina/Documents/TRAITS/data/nxr/kitzinger2021/epang_test.fasta',
             output_dir='/home/robaina/Documents/TRAITS/data/nxr/epang',
             n_threads=None,
             additional_args='--redo')

# Run hmmsearch to get gene-specific sequences

# Run cd-hit to remove uninformative sequences

# Run MSA software: mafft (parallel), or muscle (single thread)

# Run fasttree or iqtree to get newick tree file

# Prun tree (and MSA if desired) of outlier branches: treshrink

# Visualize tree with empress tree-plot (html)

# Do placement of short reads on tree