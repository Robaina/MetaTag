 #!/usr/bin/env python
# https://github.com/hyattpd/prodigal/wiki
# It runs prodigal on a set of nucleotide sequences in a folder. 
# Needs ruffus.
# Best if the sequences are cleaned with clean.py.
# Run this one first and then reformatlabel.py to modify the labels.

"""
NOTE (Semidan)
This script runs prodigal on a set of nucleotide sequences in a directory to translate them to
peptide sequences.

I understand that this script is run within cleangenomes directory, so it accesses 'files' within.

REQUIRES directory prodigal

ASKS:
-----

1. Default is genome sequence. Metagenome sequences? (yes, no)
2. How many processors?

OUTPUTS:
-------
1. Directory 'results'

2. Prodigal outputs GBK file by default


Dependencies:

1. Prodigal (installed in conda base: conda activate)
"""

import os, subprocess

global pathorigin, pathresults, nsequences


pathorigin="/home/robaina/cleangenomes/files/"

pathresults="/home/robaina/cleangenomes/results/"


cmd = ["ls -1 "+pathorigin+" | wc -l"]
pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p_status = pipe.wait()
out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))

nsequences = out.decode('ascii').strip()
print ("Number of sequences: " + nsequences)

cmd = ["rm *.gbk"]
pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p_status = pipe.wait()
out, err = pipe.communicate()

cmd = ["rm *.faa"]
pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p_status = pipe.wait()
out, err = pipe.communicate()

cmd = ["rm -r "+pathresults]
pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p_status = pipe.wait()
out, err = pipe.communicate()

cmd = ["mkdir "+pathresults]
pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p_status = pipe.wait()
out, err = pipe.communicate()

print ()
x = input("Default is genome sequence. Metagenome sequences? ")
if x=="":
	x="no"
x=x.lower()
x=x[0]

print()
xx = input("How many processors? ")
xx=int(xx)
print ()

from ruffus import *

def fnames():
	files=[]
	dirList=os.listdir(pathorigin)
	for fname in dirList:
		files.append(pathorigin+fname)

	parameters=[[0 for i in range(3)] for j in range(len(files))]
	for d1 in range(len(files)):
		parameters[d1][0]= files[d1]
		parameters[d1][1]= files[d1].replace(pathorigin,pathresults)
		parameters[d1][2]= str(d1+1)

	for job_parameters in parameters:
		yield job_parameters

@files(fnames)  # Each file gets its own process
def parallel_task(input_file, output_file, i):
	print (i+"/"+nsequences+" "+input_file.replace(pathorigin,"").replace(".fasta",""))
	fname=input_file.replace(pathorigin,"")
	if x=="y":
		cmd = ["prodigal -q -p meta -i "+pathorigin+fname+" -o /home/robaina/prodigal/"+fname.replace(".fasta","")+".gbk -a /home/robaina/prodigal/"+fname.replace(".fasta","")+".faa"]
		pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		p_status = pipe.wait()
		out, err = pipe.communicate()
		print (out.decode('ascii'), err.decode('ascii'))

	else:
		cmd = ["prodigal -q -i "+pathorigin+fname+" -o /home/robaina/prodigal/"+fname.replace(".fasta","")+".gbk -a /home/robaina/prodigal/"+fname.replace(".fasta","")+".faa"]
		pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		p_status = pipe.wait()
		out, err = pipe.communicate()

	cmd = ["mv /home/robaina/prodigal/"+fname.replace(".fasta","")+".faa "+pathresults+fname.replace(".fasta","")+".aa"]
	pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	p_status = pipe.wait()
	out, err = pipe.communicate()

	cmd = ["rm /home/robaina/prodigal/"+fname.replace(".fasta","")+".gbk"]
	pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	p_status = pipe.wait()
	out, err = pipe.communicate()
	

pipeline_run([parallel_task], verbose=1, multiprocess=xx)



print()
print ("I'm done!")


exit()
