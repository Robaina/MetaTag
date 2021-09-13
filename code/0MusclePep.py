#!/usr/bin/env python

"""
NOTES (Semidan)

Runs muscle: multiple sequence alignment
"""

import subprocess, re, os
from Bio import SeqIO
from Bio import AlignIO

cmd = ["rm alignment.*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm sequencestrimmed.*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

regex = re.compile('[^A-Z]')
regex2 = re.compile('[^GPAVLIMCFYWHKRQNEDST]')
regex3 = re.compile('[GPAVLIMCFYWHKRQNEDST]')

nameoffile="sequences"

records = list(SeqIO.parse(nameoffile, "fasta"))
numberofsequeces=len(records)
print()#; input(str(numberofsequeces)+" total sequences. ")
lennumberofsequeces=len(str(numberofsequeces))
print (lennumberofsequeces)

handle = open("sequences", "r")
filetowrite = open("sequences.fasta", "w")
filetowrite2 = open("sequencesLongLabels.fasta", "w")
filetowrite3 = open("leftovers.fasta", "w")
j=0
minseqlenght=1E6; maxseqlentgh=1; totallength=0
for i in SeqIO.parse(handle, "fasta"):

	#if len(str(i.seq))>=350-100 and len(str(i.seq))<=350+100 and len(regex2.sub("", str(i.seq).upper()))==len(str(i.seq)) and str(i.seq)[0]=="M":
	#if len(str(i.seq))>=400-100 and len(str(i.seq))<=400+100:
	#if len(regex2.sub("", str(i.seq).upper()))==len(str(i.seq)) and str(i.seq)[0]=="M":
	#if len(regex2.sub("", str(i.seq).upper()))==len(str(i.seq)):
	if True:
		j+=1
		print (i.name)
		
		filetowrite.write(">"+"0"*(lennumberofsequeces-len(str(j)))+str(j)+"_seq"+"\n")
		
		if str(i.description).find("_")>-1:
			for x in range(len(str(i.description))):
				if ord(str(i.description)[x])<48 or ord(str(i.description)[x])>57:
					break
			if str(i.description)[x]=="_":
				x+=1
				
			filetowrite2.write(">"+"0"*(lennumberofsequeces-len(str(j)))+str(j)+"_"+str(i.description)[x:]+"\n")
		else:
			filetowrite2.write(">"+"0"*(lennumberofsequeces-len(str(j)))+str(j)+"_"+str(i.description)+"\n")
			
		filetowrite.write(regex.sub("", str(i.seq).upper())+"\n")
		filetowrite2.write(regex.sub("", str(i.seq).upper())+"\n")
		if len(i.seq)>maxseqlentgh:
			maxseqlentgh=len(i.seq)
			lblmaxseqlentgh=str(i.description)
		if len(i.seq)<minseqlenght:
			minseqlenght =len(i.seq)
			lblminseqlenght=str(i.description)
		totallength+=len(i.seq)
	else:
		filetowrite3.write(">"+str(i.description)+"\n")
		filetowrite3.write(regex.sub("", str(i.seq).upper())+"\n")
filetowrite3.close()
filetowrite2.close()
filetowrite.close()
handle.close()

print()
print ("Minimum length: "+str(minseqlenght)+" ("+lblminseqlenght+")")
print ("Maximum length: "+str(maxseqlentgh)+" ("+lblmaxseqlentgh+")")
print ("Average length: "+str(totallength/j)); print()


records = list(SeqIO.parse("sequencesLongLabels.fasta", "fasta"))
numberofsequeces=len(records)
print()#; input(str(numberofsequeces)+" total sequences. ")
lennumberofsequeces=len(str(numberofsequeces))
print (lennumberofsequeces)

print(); q=input("Large dataset? ")
#q="no"
q=q.lower()[0]

cmd = ["rm *.html"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

if q=="y":
	print(); print ("muscle -in sequences.fasta -out sequences.fasta.aln -maxiters 2")
	os.system("muscle -in sequences.fasta -out sequences.fasta.aln -maxiters 2")
else:
	print(); print ("muscle -in sequences.fasta -out sequences.fasta.aln")
	os.system("muscle -in sequences.fasta -out sequences.fasta.aln")

input_handle = open("sequences.fasta.aln", "r")
output_handle = open("alignment.clustal", "w")
alignments = AlignIO.parse(input_handle, "fasta")
AlignIO.write(alignments, output_handle, "clustal")
output_handle.close()
input_handle.close()

input_handle = open("sequences.fasta.aln", "r")
output_handle = open("alignment.phylip", "w")
alignments = AlignIO.parse(input_handle, "fasta")
AlignIO.write(alignments, output_handle, "phylip")
output_handle.close()
input_handle.close()

regex = re.compile('[^-A-Z]')
minni=0
maxni=1E6
lbl=[]
sequences=[]
handle = open("sequences.fasta.aln", "r")
for i in SeqIO.parse(handle, "fasta"):
	lbl.append(str(i.name))
	a=regex.sub("", str(i.seq).upper())

	for x in range(len(a)):
		if a[x]!="-":
			break
		
	if x>minni:
		minni=x

	for x in range(1, len(a)):
		if a[-x]!="-":
			break
		
	if len(a)-x<maxni:
		print (len(a)-x)
		maxni=len(a)-x

	sequences.append(a)
		
		
print ("minni: "+str(minni))
print ("maxni: "+str(maxni))
	
handle.close()

filetowrite = open("sequencestrimmed.fasta", "w")
for j in range(len(lbl)):
	filetowrite.write(">"+lbl[j]+"\n")
	filetowrite.write(sequences[j][minni:maxni]+"\n")
filetowrite.close()

input_handle = open("sequencestrimmed.fasta", "r")
output_handle = open("sequencestrimmed.phylip", "w")
alignments = AlignIO.parse(input_handle, "fasta")
AlignIO.write(alignments, output_handle, "phylip")
output_handle.close()
input_handle.close()
	
print(); print ("I'm done!")

