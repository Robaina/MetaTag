#!/bin/bash
# chmod +x 0.sh
# Needs to have sequencestrimmed.phylip or alignment.phylip, and sequencesLongLabels.fasta.
# itolfasttree.py changes labels for longer labels in sequencesLongLabels.fasta.



rm *.txt

#fasttree sequencestrimmed.phylip > tree.txt

fasttree alignment.phylip > tree.txt

# nucleotides:
#fasttree -gtr -nt sequencestrimmed.phylip > tree_file 


python itol.py

