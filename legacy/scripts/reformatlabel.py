#!/usr/bin/env python
# Run this one after loopparallel.py.
# It mofidied the labels to have the name of the file, followed by two "_", then the
# contig number, peptide number, gene start, gene end and orientation.

import os, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

work_dir = '/home/robaina'

pathresults = '/home/robaina/cleangenomes/results/'

pathtonucleotide = '/home/robaina/cleangenomes/cleanfiles/'

pathresultsformatted="./a/"

cmd = ["rm -r "+pathresultsformatted]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir "+pathresultsformatted]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

question=input("CDS too? ")
question=question.lower()[0]

cmd = [f"rm -r {work_dir}/prodigal/anucleotides"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
if question=="y":
	cmd = [f"mkdir {work_dir}/prodigal/anucleotides"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

dirList=os.listdir(pathresults)
filenumber=0
for fname in dirList:
	filenumber+=1
	print (str(filenumber)+"     "+fname)
	handle = open(pathresults+fname, "r")
	labl=[]
	mcontig=[]
	mpeptidename=[]
	mcdstart=[]
	mcdsend=[]
	morientation=[]
	msequence=[]
	for x in SeqIO.parse(handle, "fasta"):
		labl.append(str(x.description))
		mcontig.append(str(x.description)[:str(x.description).find("_")])
		peptidenumber=str(x.description)[1+str(x.description).find("_"):]
		peptidenumber=peptidenumber[:peptidenumber.find(" ")]
		mpeptidename.append(peptidenumber)
		cdsstart=str(x.description)[1+str(x.description).find("#"):]
		cdsstart=cdsstart[:cdsstart.find("#")].strip()
		mcdstart.append(cdsstart)
		cdsend=str(x.description)[1+str(x.description).find("#"):]
		cdsend=cdsend[1+cdsend.find("#"):]
		cdsend=cdsend[:cdsend.find("#")].strip()
		mcdsend.append(cdsend)
		orientation=str(x.description)[1+str(x.description).find("#"):]
		orientation=orientation[1+orientation.find("#"):]; orientation=orientation[1+orientation.find("#"):]
		orientation=orientation[:orientation.find("#")].strip()
		if orientation=="1":
			morientation.append("pos")
		if orientation=="-1":
			morientation.append("neg")
		msequence.append(str(x.seq).replace("*",""))
	handle.close()

	longestcontig=0
	longestpeptide=0
	longeststart=0
	longestend=0
	for i in range(len(labl)):
		if len(mcontig[i])> longestcontig:
			longestcontig=len(mcontig[i])
		if len(mpeptidename[i])> longestpeptide:
			longestpeptide =len(mpeptidename[i])
		if len(mcdstart[i])> longeststart:
			longeststart =len(mcdstart[i])
		if len(mcdsend[i])> longestend:
			longestend =len(mcdsend[i])

	filetowrite = open(pathresultsformatted+fname.replace(".aa","")+".fasta", "w")
	finallbl=[]
	for i in range(len(labl)):
		a=fname[:fname.find(".")]+"__"+"0"*(longestcontig-len(mcontig[i]))+mcontig[i]
		a+="_"+"0"*(longestpeptide-len(mpeptidename[i]))+mpeptidename[i]
		a+="_"+"0"*(longeststart-len(mcdstart[i]))+ mcdstart[i]
		a+="_"+"0"*(longestend-len(mcdsend[i]))+ mcdsend[i]
		a+="_"+morientation[i]
		filetowrite.write(">"+a+"\n"); finallbl.append(a)
		filetowrite.write(msequence[i]+"\n")
	filetowrite.close()

	if question=="y":
		filetowrite = open(f"{work_dir}/prodigal/anucleotides/"+fname.replace(".aa","")+".fasta", "w")
		for i in range(len(mcontig)):
			handle = open(pathtonucleotide+fname.replace(".aa",".fasta"), "r")
			for x in SeqIO.parse(handle, "fasta"):
				if mcontig[i]==str(x.name):
					filetowrite.write(">"+finallbl[i]+"\n")
					a=str(x.seq)[int(mcdstart[i])-1:int(mcdsend[i])]
					if morientation[i]=="pos":	
						filetowrite.write(a+"\n")
					if morientation[i]=="neg":
						filetowrite.write(str(Seq(a).reverse_complement())+"\n")
					break
			handle.close()
		filetowrite.close()

print()
print ("I'm done!")
