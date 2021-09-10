#!/bin/bash
# chmod +x iqtree.sh
# Needs to have sequencestrimmed.phylip or alignment.phylip, and sequencesLongLabels.fasta.
# itoliqtree.py changes labels for longer labels in sequencesLongLabels.fasta.



#Any gene
###iqtree -s alignment.phylip -st AA -nt AUTO -m TEST -bb 1000 -o 001_seq
#iqtree -s alignment.phylip -st AA -nt AUTO -m TEST -bb 1000

rm *.phylip.*

iqtree -s alignment.phylip -st AA -nt AUTO -m TEST -bb 1000

#iqtree -s sequencestrimmed.phylip -st AA -nt AUTO -m TEST -bb 1000 -o 001_seq

#iqtree -s alignment.phylip -st AA -nt AUTO -m TEST -bb 1000



python itol.py


