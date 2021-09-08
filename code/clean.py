#!/usr/bin/env python
# Cleans fasta sequences in folder.
# To run other programs, better if you choose option 3.

import os, re, subprocess, time

nmultiprocess=5

regexlbl = re.compile('[^a-zA-Z0-9]')
regexseq = re.compile('[^A-Z]')

from Bio import SeqIO

question = input("Should I change labels for file name followed by a numbers (1) or leave as is (2) or just numbers (3)? (enter if 3) ")
question=question.strip()
if question=="":
	question="3"
question=question[0]
question=int(question)

cmd = ["rm -r /usr/gonzalez/cleangenomes/newfiles"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir /usr/gonzalez/cleangenomes/newfiles"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

path="/usr/gonzalez/cleangenomes/files/"
dirList=os.listdir(path)
for fname in dirList:
	print (fname)
	newname=fname[:fname.rfind(".")]
	newname=regexlbl.sub("_", newname)
	#newname=regexlbl.sub("", newname) 	

	while newname[-1]=="_":
		newname=newname[:-1]
		
	while newname[0]=="_":
		newname=newname[1:]
		
	print (newname)

	#newname=newname.replace(".","_")
	while newname.find("__")!=-1:
		newname=newname.replace("__","_")
	newname=newname[0].upper()+newname[1:]

	print ("cp "+path+fname+" /usr/gonzalez/cleangenomes/newfiles/"+newname+".fasta")
	cmd = ["cp "+path+fname+" /usr/gonzalez/cleangenomes/newfiles/"+newname+".fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
time.sleep(2)

cmd = ["rm -r /usr/gonzalez/cleangenomes/cleanfiles"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir /usr/gonzalez/cleangenomes/cleanfiles"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

from ruffus import *

def reformatfile(fname, question):
	path="/usr/gonzalez/cleangenomes/newfiles/"

	filetowrite = open("/usr/gonzalez/cleangenomes/cleanfiles/"+fname, "w")
	handle = open(path+"/"+fname, "r")
	j=0
	for i in SeqIO.parse(handle, "fasta"):
		j+=1
		
		if question==1:
			filetowrite.write(">"+fname[:fname.rfind(".")]+"_"+str(j)+"\n")
		if question==2:
			filetowrite.write(">"+regexlbl.sub("_", str(i.description))+"\n")
		if question==3:
			filetowrite.write(">"+str(j)+"\n")
			
		filetowrite.write(regexseq.sub("", str(i.seq).upper())+"\n")
		if len(regexseq.sub("", str(i.seq).upper()))==0:
			print (j)
			print (fname)
			print ("No sequence data. I'm leaving.")
			exit()
	handle.close()
	filetowrite.close()

def fnames():
	path="/usr/gonzalez/cleangenomes/newfiles/"
	files=[]
	dirList=os.listdir(path)
	for fname in dirList:
		files.append(path+fname)

	parameters=[[0 for i in range(3)] for j in range(len(files))]
	for d1 in range(len(files)):
		parameters[d1][0]= files[d1]
		parameters[d1][1]= files[d1].replace(path,"/usr/gonzalez/cleangenomes/cleanfiles/")
		parameters[d1][2]= str(d1+1)

	for job_parameters in parameters:
		yield job_parameters

@files(fnames)
def parallel_task(input_file, output_file, i):
	print (i+" "+input_file.replace("/usr/gonzalez/cleangenomes/newfiles/",""))
	fname=input_file.replace("/usr/gonzalez/cleangenomes/newfiles/","")
	reformatfile(fname, question)
	

pipeline_run([parallel_task], verbose=2, multiprocess=nmultiprocess)

cmd = ["rm -r /usr/gonzalez/cleangenomes/newfiles"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

print("Files saved in /usr/gonzalez/cleangenomes/cleanfiles/")
print ("I'm done!")
