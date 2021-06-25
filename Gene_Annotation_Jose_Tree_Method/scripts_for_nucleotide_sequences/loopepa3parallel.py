#!/usr/bin/env python
# Have all sequences in ./hits.
# For each file, it will launch epa3parallel.py
# Launch 0MusclePep.py. Take alignment.phylip and sequencesLongLabels.fasta to ./iqtree.
# Go to ./iqtree and then:
# rm *.phylip.*
# iqtree -s alignment.phylip -st AA -nt AUTO -m TEST -bb 1000
# python itol.py
# Check file alignment.phylip.log to get best model in line:
# Best-fit model: LG+I+G4 chosen according to BIC
# and change the line in file epa3parallel.py

import re, os, subprocess, time
from Bio import SeqIO

cmd = ["rm -r ./resultclusterstemp/"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./resultclusterstemp/"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r ./removedseqtemp/"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./removedseqtemp/"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r ./resultsseq/"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./resultsseq/"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()


tosplit=500

dirList=os.listdir("/home/gonzalez/epa-ng/hits")
nfiles=0
for ffile in dirList:
	nfiles+=1
    
	cmd = ["rm -r ./hitsdirectory"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["mkdir ./hitsdirectory"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

	cmd = ["grep -c '^>' /home/gonzalez/epa-ng/hits/"+ffile]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (str(int(out.decode('ascii').strip()))+" sequences.")
	
	handle = open("/home/gonzalez/epa-ng/hits/"+ffile, "r")
	numberofseq=0
	j=0
	filetowrite = open("./hitsdirectory/"+ffile.replace(".fasta","")+"_"+str(j)+".fasta", "w")
	for xx in SeqIO.parse(handle, "fasta"):
		numberofseq+=1
		if numberofseq/tosplit==int(numberofseq/tosplit):
			print (numberofseq)
			j+=1
			filetowrite.close()
			filetowrite = open("./hitsdirectory/"+ffile.replace(".fasta","")+"_"+str(j)+".fasta", "w")
		filetowrite.write(">"+str(xx.name)+"\n")
		filetowrite.write(str(xx.seq)+"\n")

	handle.close()
	filetowrite.close()

	os.system("python epa3parallel2.py "+ffile)

	####################################################################################################
	filetoread = open("resultsclusters.txt")
	filetowrite = open("./resultclusterstemp/resultsclusters"+str(nfiles)+".txt", "w")
	
	mat=[]
	for line in filetoread:
		m=line.split()
		#print (m)
		if m[0]=="File":
			filetowrite.write(line)
			for i in range(len(m)):
				mat.append(0)	
			
		else:
			for i in range(1, len(m)):
				mat[i]+=int(m[i])
	
	filetowrite.write(ffile+"\t")
	
	towrite=""
	for i in range(1, len(mat)):
		towrite+=str(mat[i])+"\t"
	
	#print (mat)
	towrite=towrite[:-1]+"\n"
	filetowrite.write(towrite)
	
	filetoread.close()
	filetowrite.close()
	####################################################################################################
	filetoread = open("removedseq.txt")
	filetowrite = open("./removedseqtemp/removedseq"+str(nfiles)+".txt", "w")
	
	mat=[]; mat.append(0); mat.append(0); mat.append(0)
	for line in filetoread:
		m=line.split()
		#print (m)
		for i in range(1, len(m)):
			#print (m[i])
			mat[i-1]+=int(m[i])
	
	filetowrite.write(ffile+"\t")
	
	towrite=""
	for i in range(1, len(m)):
		towrite+=str(mat[i-1])+"\t"
	
	#print (mat)
	towrite=towrite[:-1]+"\n"
	filetowrite.write(towrite)
	
	filetoread.close()
	filetowrite.close()

	####################################################################################################

cmd = ["cat /home/gonzalez/epa-ng/removedseqtemp/*.txt > removedseq.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

filetowrite = open("resultsclusters.txt", "w")
act=1
dirList=os.listdir("./resultclusterstemp")
for ffile in dirList:
	filetoread = open("./resultclusterstemp/"+ffile)
	for line in filetoread:
		if line.startswith("File") and act==1:
			filetowrite.write(line)
			act=0
		if not(line.startswith("File")):
			filetowrite.write(line)
	filetoread.close()
filetowrite.close()



cmd = ["rm -r ./removedseqtemp"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm -r ./resultclusterstemp"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm -r ./hitsdirectory"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()


print()
print ("I'm done!")
exit()

