#!/usr/bin/env python
# Each of the peptides in sequencesLongLabels.fasta should have a label between brackets with the 
# format [Cluster X]. This way, one can automatically quantify each of the groups. 
# If it doesn't have any, then run this to add a generic label. Otherwise the papara script won't work.
# After you run it, do:
# mv sequencesLongLabels2.fasta sequencesLongLabels.fasta

from Bio import SeqIO

lbl=" [Cluster X]"

handle = open("sequencesLongLabels.fasta", "r")
filetowrite = open("sequencesLongLabels2.fasta", "w")

for xx in SeqIO.parse(handle, "fasta"):
	filetowrite.write(">"+str(xx.description)+lbl+"\n")
	filetowrite.write(str(xx.seq)+"\n")
filetowrite.close()
handle.close()


print (""); print ("I'm done!")
exit()