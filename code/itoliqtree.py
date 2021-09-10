#!/usr/bin/env python
#

import re, os
from Bio import SeqIO



if os.path.exists("alignment.phylip.contree"):
	treefiletoread="alignment.phylip.contree"

if os.path.exists("sequencestrimmed.phylip.contree"):
	treefiletoread="sequencestrimmed.phylip.contree"

regex = re.compile('[^a-zA-Z0-9]')

handle = open("sequencesLongLabels.fasta", "r")
lbl=[]; lbllong=[]
sseq=[]
for loop in SeqIO.parse(handle, "fasta"):
	print (str(loop.description)[:str(loop.description).find("_")])
	lbl.append(str(loop.description)[:str(loop.description).find("_")])
	a=str(loop.description)
	a=regex.sub("_", a)
	while a.count("__")!=0:
		a=a.replace("__", "_")
	
	lbllong.append("'"+a+"'")
	sseq.append(str(loop.seq))
	
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
exit()

for i in range(len(lbllong)):
	lbllong[i]=lbllong[i].replace("'","").strip()

	if lbllong[i].find(" [Cluster")>-1:
		lbllong[i]=lbllong[i][:lbllong[i].find(" [Cluster")]

tree = open(treefiletoread, 'r').read()

for i in range(len(lbl)):
	tree=tree.replace(lbl[i]+"_seq", "'"+lbllong[i]+" [Cluster rubiscolarge]'")
    

print (tree)

filetowrite = open("0newtreelonger2.txt", "w")
filetowrite.write(tree+"\n")
filetowrite.close()

filetowrite = open("sequencesLongLabels2.fasta", "w")
for i in range(len(lbllong)):
	filetowrite.write(">"+lbllong[i]+" [Cluster rubiscolarge]\n")
	filetowrite.write(sseq[i]+"\n")
	
filetowrite.close()







