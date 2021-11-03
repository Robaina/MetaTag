#!/usr/bin/env python
# Launch the hmm first at ./hmmer with pfamgenomesparallel.py.
# Dont launch hmm first. Have the problem sequences in the folder /home/gonzalez/epa-ng/problemseqs.

import os, subprocess, re
from Bio import SeqIO

regex = re.compile('[^a-zA-Z0-9]')

print ("cat /home/gonzalez/epa-ng/problemseqs/*.fasta > ttemp.fasta")
cmd = ["cat /home/gonzalez/epa-ng/problemseqs/*.fasta > ttemp.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

tosplit=200

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

cmd = ["rm -r ./resultsclusterslbl"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./resultsclusterslbl"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

os.system("python peptidesepa3parallel.py ttemp.fasta")

#cmd = ["rm ttemp.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm results.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

print ("tar -zcvf resultsclusterslbl.tar.gz -C ./resultsclusterslbl .")
cmd = ["tar -zcvf resultsclusterslbl.tar.gz -C ./resultsclusterslbl ."]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))

print ("gunzip -t resultsclusterslbl.tar.gz")
cmd = ["gunzip -t resultsclusterslbl.tar.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))

print ("gunzip -c resultsclusterslbl.tar.gz | tar t > /dev/null")
cmd = ["gunzip -c resultsclusterslbl.tar.gz | tar t > /dev/null"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))

os.system("python /home/gonzalez/send_email.py")
print ()
print ("I'm done!")
exit()
