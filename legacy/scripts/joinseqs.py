#!/usr/bin/env python
# Reads fasta files in ./files and ignores sequences with repeated labels or sequences.
# It also removes sequences that don't start with an M. Removes numbers at the beginning of label.
# It sorts sequence files before anything. If the file label contains the word complete, it would
# insert that word before the label of the sequences.

import os, re
from Bio import SeqIO

ppath="./files/"

regex = re.compile('[^GPAVLIMCFYWHKRQNEDST]')

dirList=os.listdir(ppath)
mfname=[]
for fname in dirList:
	mfname.append(fname)
	
mfname.sort()
print (mfname)


lbl=[]; seqs=[]
filetowrite = open("seqs.fasta", "w")
i=0; nfile=0
for x in range(len(mfname)):
	nfile+=1
	print (mfname[x])
	handle = open(ppath+mfname[x], "r")
	j=0; duplicates=0
	for xx in SeqIO.parse(handle, "fasta"):
		j+=1
		#if not(str(xx.name) in lbl):
		if not(str(xx.name) in lbl) and not(str(xx.seq).upper().strip() in seqs) and len(regex.sub("", str(xx.seq).upper()))==len(str(xx.seq)) and str(xx.seq)[0]=="M":
		##if not(str(xx.name) in lbl) and not(str(xx.seq).upper().strip() in seqs):
		#if not(str(xx.seq).upper().strip() in seqs):
			lbl.append(str(xx.name))
			seqs.append(str(xx.seq).upper().strip())
			i+=1
			
			b=str(xx.description)
			if str(xx.description).find("_")>-1:
				for y in range(len(str(xx.description))):
					if ord(str(xx.description)[y])<48 or ord(str(xx.description)[y])>57:
						break
				if str(xx.description)[y]=="_":
					y+=1
					
				b=str(xx.description)[y:]

			if mfname[x].lower().find("complete")>-1:
				filetowrite.write(">complete"+b+"\n"+str(xx.seq)+"\n")
			else:

				filetowrite.write(">"+b+"\n"+str(xx.seq)+"\n")

		else:
			duplicates+=1
	handle.close()
	print (j)
	print ("Bad sequences: "+str(duplicates))

filetowrite.close()

print()
print ("Total sequences: "+str(i))
exit()
print()
input("Keep going? ")

filetoread = open("taxa.txt", "r")
genomes=[]
for line in filetoread:
	line=line.strip()
	if len(line)>0:
		genomes.append(line[:line.find(".")])
		
print (genomes)

filetowrite = open("results.txt", "w")
for i in range(len(genomes)):
	print (genomes[i])
	filetowrite.write(genomes[i]+"\t")
	j=0
	handle = open("seqs.fasta", "r")
	for xx in SeqIO.parse(handle, "fasta"):
		if str(xx.name)[:str(xx.name).find("__")]==genomes[i]:
			j+=1
	handle.close()
	print (j)
	filetowrite.write(str(j)+"\n")
filetowrite.close()
exit()

