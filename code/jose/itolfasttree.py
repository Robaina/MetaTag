#!/usr/bin/env python

"""
NOTES (Semidán)

Run this script after running fasttree to make iTOL (tree visualization) - compatible input file

"""

import re
from Bio import SeqIO


treefiletoread="tree.txt"

regex = re.compile('[^a-zA-Z0-9]')

handle = open("sequencesLongLabels.fasta", "r")
lbl=[]; lbllong=[]
for loop in SeqIO.parse(handle, "fasta"):
	print (str(loop.description)[:str(loop.description).find("_")])
	lbl.append(str(loop.description)[:str(loop.description).find("_")])
	a=str(loop.description)
	a=regex.sub("_", a)
	while a.count("__")!=0:
		a=a.replace("__", "_")
	
	lbllong.append("'"+a+"'")
	
handle.close()

print (lbl)
print (lbllong)


tree = open(treefiletoread, 'r').read()

for i in range(len(lbl)):
	tree=tree.replace(lbl[i]+"_seq", lbllong[i])

print (tree)

filetowrite = open("0newtreelonger.txt", "w")
filetowrite.write(tree+"\n")
filetowrite.close()


print(); print ("I'm done!")




