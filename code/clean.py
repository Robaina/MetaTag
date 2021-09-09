#!/usr/bin/env python
# Cleans fasta sequences in folder.
# To run other programs, better if you choose option 3.


"""
NOTES (Semidán):

{work_dir}/cleangenomes/files: required to parse input files

"""

import os, re, subprocess, time
from Bio import SeqIO
import ruffus

global work_dir

work_dir = '/home/robaina'

nmultiprocess = 5

regexlbl = re.compile('[^a-zA-Z0-9]')
regexseq = re.compile('[^A-Z]')

question = input("Should I change labels for file name followed by a numbers (1) or leave as is (2) or just numbers (3)? (enter if 3) ")
question = question.strip()
if question=="":
	question="3"
question=question[0]
question=int(question)

newfiles_dir = os.path.join(work_dir, 'cleangenomes', 'newfiles')
if os.path.isdir(newfiles_dir):
	cmd = [f"rm -r {newfiles_dir}"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = [f"mkdir {newfiles_dir}"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

path = os.path.join(work_dir, 'cleangenomes', 'files')
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
	while newname.find("__") != -1:
		newname=newname.replace("__","_")
	newname=newname[0].upper() + newname[1:]

	print ("cp " + path + fname + " " + newfiles_dir + newname + ".fasta")
	cmd = [f"cp {os.path.join(path, fname)} {os.path.join(newfiles_dir, newname)}.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
time.sleep(2)

cleanfiles_dir = os.path.join(work_dir, 'cleangenomes', 'cleanfiles')
# cmd = [f"rm -r {cleanfiles_dir}"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = [f"mkdir {cleanfiles_dir}"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()


def reformatfile(fname, question):
	# path="/usr/gonzalez/cleangenomes/newfiles/"
	cleanfiles_dir = '/home/robaina/cleangenomes/cleanfiles'
	print('path_clean_hola', os.path.join(cleanfiles_dir, fname))
	#filetowrite = open(os.path.join(cleanfiles_dir, fname), "w")
	filetowrite = open('/home/robaina/cleangenomes/newfiles' + fname)
	handle = open(os.path.join(newfiles_dir, fname), "r")
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
	# path="/usr/gonzalez/cleangenomes/newfiles/"
	files=[]
	dirList=os.listdir(newfiles_dir)
	for fname in dirList:
		files.append(os.path.join(newfiles_dir, fname))

	parameters=[[0 for i in range(3)] for j in range(len(files))]
	for d1 in range(len(files)):
		parameters[d1][0]= files[d1]
		parameters[d1][1]= files[d1].replace(newfiles_dir, cleanfiles_dir)
		parameters[d1][2]= str(d1+1)

	for job_parameters in parameters:
		yield job_parameters

@ruffus.files(fnames)
def parallel_task(input_file, output_file, i):
	print (i+" "+input_file.replace(newfiles_dir,""))
	fname=input_file.replace(newfiles_dir,"")
	reformatfile(fname, question)
	

ruffus.pipeline_run([parallel_task], verbose=2, multiprocess=nmultiprocess)

cmd = [f"rm -r {newfiles_dir}"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

print(f"Files saved in {cleanfiles_dir}")
print ("I'm done!")
