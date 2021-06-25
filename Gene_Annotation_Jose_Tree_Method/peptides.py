#!/usr/bin/env python
# Launch the hmm first at ./hmmer with pfamgenomesparallel.py.

import os, subprocess, re
from Bio import SeqIO

regex = re.compile('[^a-zA-Z0-9]')

cmd = ["rm -r ./DB"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./DB"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

print ("makeblastdb -in sequencesLongLabels.fasta -dbtype prot -out ./DB/temp")
cmd = ["makeblastdb -in sequencesLongLabels.fasta -dbtype prot -out ./DB/temp"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))

cmd = ["rm -r ./sequencestoblast"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./sequencestoblast"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

#print ("cat /home/gonzalez/hmmer/filesMegahit/*.fasta > temp.fasta")
#cmd = ["cat /home/gonzalez/hmmer/filesMegahit/*.fasta > temp.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print ("cat /home/gonzalez/hmmer/files/*.fasta > temp.fasta")
cmd = ["cat /home/gonzalez/hmmer/files/*.fasta > temp.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
#print ("cp K03520_aa.fasta temp.fasta")
#cmd = ["cp K03520_aa.fasta temp.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

print ("grep -c '^>' temp.fasta")
cmd = ["grep -c '^>' temp.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (regex.sub("", out.decode('ascii'))+" sequences.")
numberofseq=int(out)

nproc=int(open("nproc.txt", 'r').read())

tosplit=int(numberofseq/nproc)
print (str(int(1+numberofseq/tosplit))+" fasta files.")

handle = open("temp.fasta", "r")
j=1; numberofseq=0; numberofzeros=10+len(str(1+numberofseq/tosplit))
filetowrite = open("./sequencestoblast/"+(numberofzeros-len(str(j)))*"0"+str(j)+".fasta", "w")
print ("Writing "+(numberofzeros-len(str(j)))*"0"+str(j)+".fasta")
for xx in SeqIO.parse(handle, "fasta"):
	numberofseq+=1
	if numberofseq/tosplit==int(numberofseq/tosplit):
		j+=1
		filetowrite.close()
		print ("Writing "+(numberofzeros-len(str(j)))*"0"+str(j)+".fasta")
		filetowrite = open("./sequencestoblast/"+(numberofzeros-len(str(j)))*"0"+str(j)+".fasta", "w")
	filetowrite.write(">"+str(xx.description)+"\n")
	filetowrite.write(str(xx.seq)+"\n")
filetowrite.close()
handle.close()

print ("parallel -j "+str(nproc+2)+" 'FASTA={}; blastp -query $FASTA -db ./DB/temp -evalue 0.0001 -outfmt 6 -num_threads 1 -out ${FASTA%}.blast' ::: ./sequencestoblast/*.fasta")
cmd = ["parallel -j "+str(nproc+2)+" 'FASTA={}; blastp -query $FASTA -db ./DB/temp -evalue 0.0001 -outfmt 6 -num_threads 1 -out ${FASTA%}.blast' ::: ./sequencestoblast/*.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))

dirList=os.listdir("./sequencestoblast")
ffiles=[]
for fname in dirList:
	if not(fname.find(".blast")>0):
		ffiles.append(fname)
	
ffiles.sort()

filetowrite = open("results.txt", "w")
filetowrite.close()

for x in range(len(ffiles)):
	print ("cat ./sequencestoblast/"+ffiles[x]+".blast >> results.txt")
	cmd = ["cat ./sequencestoblast/"+ffiles[x]+".blast >> results.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r ./DB"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm -r ./sequencestoblast"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

filetoread = open("results.txt")
hits=[]
for line in filetoread:
	line=line.strip()
	items = line.split("\t")
	#if float(items[2])>=90 and float(items[-1])>=50:
	if float(items[2])>=80 and float(items[-1])>=50:
		if not(items[0] in hits):
			print (line)
			hits.append(items[0])
filetoread.close()

print (hits)

filetowrite = open("ttemp.fasta", "w")
handle = open("temp.fasta", "r")
lbl=[]; seqs=[]
for xx in SeqIO.parse(handle, "fasta"):
	if str(xx.name) in hits:
		filetowrite.write(">"+str(xx.name)+"\n")
		filetowrite.write(str(xx.seq)+"\n")
		
		lbl.append(str(xx.name))
		seqs.append(str(xx.seq).upper().strip())
		
filetowrite.close()

filetowrite = open("ttemp.fasta", "a")
handle = open("/home/gonzalez/hmmer/results.fasta", "r")
i=0
for xx in SeqIO.parse(handle, "fasta"):
	if not(str(xx.name) in lbl):
		filetowrite.write(">"+str(xx.name)+"\n")
		filetowrite.write(str(xx.seq)+"\n")
		i+=1
handle.close()
filetowrite.close()

input ("Extra sequences retrieved by hmm: "+str(i))

cmd = ["rm temp.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

tosplit=500

cmd = ["rm -r ./hitsdirectory"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./hitsdirectory"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["grep -c '^>' ttemp.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (str(int(out.decode('ascii').strip()))+" sequences.")
	
handle = open("ttemp.fasta", "r")
numberofseq=0
j=0
filetowrite = open("./hitsdirectory/"+str(j)+".fasta", "w")
for xx in SeqIO.parse(handle, "fasta"):
	numberofseq+=1
	if numberofseq/tosplit==int(numberofseq/tosplit):
		print (numberofseq)
		j+=1
		filetowrite.close()
		filetowrite = open("./hitsdirectory/"+str(j)+".fasta", "w")
	filetowrite.write(">"+str(xx.name)+"\n")
	filetowrite.write(str(xx.seq)+"\n")

handle.close()
filetowrite.close()

cmd = ["rm -r ./resultsseq"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./resultsseq"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

os.system("python peptidesepa3parallel.py ttemp.fasta")

cmd = ["rm ttemp.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm results.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()


print ()
print ("I'm done!")
exit()
