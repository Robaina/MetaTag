#!/usr/bin/env python
# 



import re, os, subprocess, time, sys
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

##################################################################################################################################

nmultiprocess=20


nproc=1 # epa-ng. 0 for  maximum number of threads available.

# result of iqtree
aamodel="mtZOA+G4" # archaeal AmoA
aamodel="LG+F+I+G4" # rubisco large
aamodel="LG+F+I+G4" # 4-hydroxybutyryl-CoA dehydratase, https://doi.org/10.1073/pnas.1402028111 
aamodel="LG+I+G4" # RadA 
aamodel="LG+F+I+G4" # DmdA

cmd = ["grep 'Best-fit model:' alignment.phylip.log"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))
aamodel=out.decode('ascii')[16:]
aamodel=aamodel[:aamodel.find(" ")].strip()

print (aamodel)

readspath="/home/gonzalez/epa-ng/hitsdirectory/"

pathresults="/home/gonzalez/epa-ng/results/"

# Needs these files:
# sequencesLongLabels.fasta, alignment.phylip.contree, alignment.phylip 

##################################################################################################################################

ffile=sys.argv[1]

#filetowrite = open("./resultsseq/"+ffile, "w")
#filetowrite.close()

cmd = ["rm removedseq*.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm resultsclusters*.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm temp*.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm del?_*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm *.job*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm query*.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm reference*.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm -r outdir*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm -r temp*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r "+pathresults]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir "+pathresults]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r resultsseqtemp"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir resultsseqtemp"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["ls -1 "+readspath+" | wc -l"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))
nsequences=out.decode('ascii').strip()
print ("Number of sequences: "+nsequences)

from Bio import Phylo
'''
tree = Phylo.read('alignment.phylip.contree', 'newick')

#print (tree.total_branch_length())
termnames = [term.name for term in tree.get_terminals()]

print (len(termnames))
print (tree.total_branch_length()/len(termnames))

print (tree.get_terminals())
print (tree.get_nonterminals())
print (tree.depths())
print ()
print (tree.distance("33_seq"))

maxbranchlength=0
for i in range(len(termnames)):
	if tree.distance(termnames[i])>maxbranchlength:
		maxbranchlength=tree.distance(termnames[i])
		print (termnames[i], maxbranchlength)

exit()
	
print (branchlength)


maxbranchlength=max(branchlength)


print(maxbranchlength)

exit()


nonterm=tree.get_nonterminals()

for i in range(len(nonterm)):
	print (i, nonterm[i])
	input("")

exit()

termnames = [term.name for term in tree.get_terminals()]

leafnumber=tree.count_terminals()

filetowrite = open("temp.txt", "w")
filetowrite.write(str(tree))
filetowrite.close()

branchlength=[]
for i in range(len(termnames)):
	branchlength.append(0)

filetoread = open("temp.txt")
for line in filetoread:
	if line.find("name=")>-1:
		a=line[6+line.find("name="):]
		a=a[:a.find("'")]

		b=line[14+line.find("branch_length="):]
		b=b[:b.find(",")]

		branchlength[termnames.index(a)]=float(b)

filetoread.close()

print (termnames)
print (branchlength)

maxbranchlength=max(branchlength)
print(maxbranchlength)

#print (tree.total_branch_length()/len(termnames))
exit()


input(" ?")
totallength=0
for i in range(len(branchlength)):
	totallength+=branchlength[i]

averagelength=totallength/len(branchlength)

print (averagelength)
'''
regexlbl = re.compile('[^a-zA-Z0-9]')
regex2 = re.compile('[^GPAVLIMCFYWHKRQNEDST]')

filetoread = open("sequencesLongLabels.fasta")
lbl=[]; lbllong=[]; lbllongcomplete=[] 
for line in filetoread:
	if line.startswith(">"):
		lbllongcomplete.append(line[1:])
		line=line[1:1+line.find("]")]
		
		lbl.append(line[:line.find("_")]+"_seq")

		a=regexlbl.sub("_", line)
		while a.count("__")!=0:
			a=a.replace("__", "_")
	
		lbllong.append("'"+a+"'")

filetoread.close()

clusters=[]
for i in range(len(lbllongcomplete)):
	a=regexlbl.sub("_", lbllongcomplete[i].split("[")[-1])
	while a[-1]=="_":
		a=a[:-1]
		
	if not(a in clusters):
		clusters.append(a)
		
clusters.sort()

filetowrite = open("resultsclusters", "w")
towrite="File\t"
for i in range(len(clusters)):
	towrite+=clusters[i]+"\t"

towrite=towrite[:-1]+"\n"
filetowrite.write(towrite)
filetowrite.close()

print (clusters)

from Bio import AlignIO

def loop(input_file, output_file, filenumber, aamodel, nproc, clusters, ffile):
	fname=input_file.split("/")[-1]
	print (fname)

	filenumber=str(filenumber)
	
	cmd = ["rm *.job"+filenumber]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm -r outdir"+filenumber]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm query"+filenumber+".fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm reference"+filenumber+".fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm -r ./temp"+filenumber]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["mkdir ./temp"+filenumber]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

	filetowrite = open("del1_"+filenumber+".fasta", "w")
	handle = open(input_file, "r")
	metaseq1=[]
	metaseq2=[]
	jj=0
	for xx in SeqIO.parse(handle, "fasta"):
		jj+=1
		#filetowrite.write(">'META"+str(jj)+"'\n")
		filetowrite.write(">'META"+str(xx.name)+"'\n")
		filetowrite.write(regex2.sub("", str(xx.seq))+"\n")                           
		#metaseq1.append("'META"+str(jj)+"'")
		metaseq1.append("'META"+str(xx.name)+"'")
		metaseq2.append(str(xx.name))
	handle.close()
    
	filetowrite.close()

	print ("./papara_static_x86_64 -a -r -t alignment.phylip.contree -s alignment.phylip -q del1_"+filenumber+".fasta -n job"+filenumber)
	os.system("./papara_static_x86_64 -a -r -t alignment.phylip.contree -s alignment.phylip -q del1_"+filenumber+".fasta -n job"+filenumber)

	time.sleep(1)

	filetoread = open("papara_alignment.job"+filenumber)
	filetowrite1 = open("reference"+filenumber+".fasta", "w")
	filetowrite2 = open("query"+filenumber+".fasta", "w")
	for line in filetoread:
		if regexlbl.sub("", line).startswith("META"):
			filetowrite2.write(">"+line[:line.find(" ")].strip()+"\n")
			filetowrite2.write(line[line.find(" "):].strip()+"\n")
		if line.find("_seq")>-1:
			filetowrite1.write(">"+line[:line.find(" ")].strip()+"\n")
			filetowrite1.write(line[line.find(" "):].strip()+"\n")

	filetowrite1.close()
	filetowrite2.close()
	filetoread.close()

	print("epa-ng -T "+nproc+" --redo --tree alignment.phylip.contree --ref-msa reference"+filenumber+".fasta --query query"+filenumber+".fasta --out-dir ./temp"+filenumber+" --model "+aamodel)
	os.system("epa-ng -T "+nproc+" --redo --tree alignment.phylip.contree --ref-msa reference"+filenumber+".fasta --query query"+filenumber+".fasta --out-dir ./temp"+filenumber+" --model "+aamodel) 
	time.sleep(1)
	print ("gappa examine graft --fully-resolve --out-dir ./outdir"+filenumber+" --jplace-path /home/gonzalez/epa-ng/temp"+filenumber)
	os.system("gappa examine graft --fully-resolve --out-dir ./outdir"+filenumber+" --jplace-path /home/gonzalez/epa-ng/temp"+filenumber)
	time.sleep(1)

	#---------------------------------------------------------------------------------------------------------

	tree = Phylo.read("./outdir"+filenumber+"/epa_result.newick", "newick")
	
	termnames = [term.name for term in tree.get_terminals()]
	
	maxbranchlength=0
	for i in range(len(termnames)):
		if termnames[i].endswith("_seq"):
			if tree.distance(termnames[i])>maxbranchlength:
				maxbranchlength=tree.distance(termnames[i])
				#print (termnames[i], maxbranchlength)
				
	print (maxbranchlength)

	badsequences=[]	
	for i in range(len(termnames)):
		if termnames[i].startswith("META"):
			if tree.distance(termnames[i])>maxbranchlength*2:
				badsequences.append(termnames[i])

	print (badsequences)
	'''
	input("?")
	exit()
	
	badsequences=[]	
	filetoread = open("temp"+filenumber+".txt")
	for line in filetoread:
		if line.find("name=")>-1:
			a=line[6+line.find("name="):]
			a=a[:a.find("'")]
	
			b=line[14+line.find("branch_length="):]
			b=b[:b.find(",")]
	
			if float(b)>2*maxbranchlength:
				badsequences.append(regexlbl.sub("", a))
	filetoread.close()

	
	exit()
	

	
	#leafnumber=tree.count_terminals()
	
	filetowrite = open("temp"+filenumber+".txt", "w")
	filetowrite.write(str(tree))
	filetowrite.close()
	
	badsequences=[]	
	filetoread = open("temp"+filenumber+".txt")
	for line in filetoread:
		if line.find("name=")>-1:
			a=line[6+line.find("name="):]
			a=a[:a.find("'")]
	
			b=line[14+line.find("branch_length="):]
			b=b[:b.find(",")]
	
			if float(b)>2*maxbranchlength:
				badsequences.append(regexlbl.sub("", a))
	filetoread.close()
    
	print (badsequences)
	'''
	cmd = ["rm *.job"+filenumber]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm -r outdir"+filenumber]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm query"+filenumber+".fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm reference"+filenumber+".fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm -r ./temp"+filenumber]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["mkdir ./temp"+filenumber]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	#exit()
	filetowrite = open("del2_"+filenumber+".fasta", "w")
	handle = open("del1_"+filenumber+".fasta", "r")
	jj=0
	for xx in SeqIO.parse(handle, "fasta"):
		if not(regexlbl.sub("", str(xx.name)) in badsequences) and len(str(xx.seq))>5:
			jj+=1
			filetowrite.write(">"+str(xx.name)+"\n")
			filetowrite.write(str(xx.seq)+"\n")
	handle.close()
	filetowrite.close()
	
	tccounter=0

	if jj>0:

		print ("./papara_static_x86_64 -a -r -t alignment.phylip.contree -s alignment.phylip -q del2_"+filenumber+".fasta -n job"+filenumber)
		os.system("./papara_static_x86_64 -a -r -t alignment.phylip.contree -s alignment.phylip -q del2_"+filenumber+".fasta -n job"+filenumber)
		time.sleep(1)
	
		filetoread = open("papara_alignment.job"+filenumber)
		filetowrite1 = open("reference"+filenumber+".fasta", "w")
		filetowrite2 = open("query"+filenumber+".fasta", "w")
		for line in filetoread:
			if regexlbl.sub("", line).startswith("META"):
				filetowrite2.write(">"+line[:line.find(" ")].strip()+"\n")
				filetowrite2.write(line[line.find(" "):].strip()+"\n")
			if line.find("_seq")>-1:
				filetowrite1.write(">"+line[:line.find(" ")].strip()+"\n")
				filetowrite1.write(line[line.find(" "):].strip()+"\n")
	
		filetowrite1.close()
		filetowrite2.close()
		filetoread.close()
	
		print ("epa-ng -T "+nproc+" --tree alignment.phylip.contree --ref-msa reference"+filenumber+".fasta --query query"+filenumber+".fasta --out-dir ./temp"+filenumber+" --model "+aamodel)
		os.system("epa-ng -T "+nproc+" --tree alignment.phylip.contree --ref-msa reference"+filenumber+".fasta --query query"+filenumber+".fasta --out-dir ./temp"+filenumber+" --model "+aamodel) 
		time.sleep(1)
		print ("gappa examine graft --fully-resolve --out-dir ./outdir"+filenumber+" --jplace-path /home/gonzalez/epa-ng/temp"+filenumber)
		os.system("gappa examine graft --fully-resolve --out-dir ./outdir"+filenumber+" --jplace-path /home/gonzalez/epa-ng/temp"+filenumber)
		time.sleep(1)
		s = open("./outdir"+filenumber+"/epa_result.newick", "r").read()
		
		#---------------------------------------------------------------------------------------------------------
		print (s)
		
		for i in range(len(lbl)):
			s=s.replace(lbl[i], "'"+regexlbl.sub("", lbllong[i])+"'")
		print (s)
	
		filetowrite = open(output_file, "w")
		filetowrite.write(s)
		filetowrite.close()
	
		#---------------------------------------------------------------------------------------------------------
			
		tree = Phylo.read("./outdir"+filenumber+"/epa_result.newick", "newick")

		termnames = [term.name for term in tree.get_terminals()]
		
		maxbranchlength=0
		for i in range(len(termnames)):
			if termnames[i].endswith("_seq"):
				if tree.distance(termnames[i])>maxbranchlength:
					maxbranchlength=tree.distance(termnames[i])
					#print (termnames[i], maxbranchlength)
					
		print (maxbranchlength)

		tree = Phylo.read('alignment.phylip.contree', 'newick')
		termnames = [term.name for term in tree.get_terminals()]

	
		handle = open("del2_"+filenumber+".fasta", "r")
		closesthittree=[]
		metaseq=[]
		tree = Phylo.read("./outdir"+filenumber+"/epa_result.newick", "newick")
		for xx in SeqIO.parse(handle, "fasta"):
			print("-----------------"); print (str(xx.name))
			mindistance=1E6
			for i in range(len(termnames)):
				#print (termnames[i])
				#print (str(xx.name).replace("'",""))
	
				dist=tree.distance(termnames[i], str(xx.name).replace("'",""))
	
				if dist<mindistance:
					mindistance=dist
					a=termnames[i]
	
			#print (a)
			#print (mindistance)
			
			#print (tree.distance(str(xx.name).replace("'","")))
			#print (maxbranchlength)
			#print (tree.distance(str(xx.name).replace("'",""))<maxbranchlength)
			#input("? ")
			#if mindistance<2*maxbranchlength:
			if tree.distance(str(xx.name).replace("'",""))<maxbranchlength*2:
				a=regexlbl.sub("_", lbllongcomplete[lbl.index(a)].split("[")[-1])
				while a[-1]=="_":
					a=a[:-1]
					
				print (a)
				
				closesthittree.append(a)
				metaseq.append(str(xx.name))
				
				filetowriteresultsclusterslbl = open("./resultsclusterslbl/"+a+".fasta", "a")
				filetowriteresultsclusterslbl.write(">"+str(xx.name)[5:].replace("'", "")+"\n"+str(xx.seq)+"\n")
				filetowriteresultsclusterslbl.close()
				
		handle.close()
		#input("? ")
		#########################################################################################################
		
		print (metaseq)
		print (len(metaseq))
		
		metaseqoriginal=[]
		for i in range(len(metaseq)):
			for j in range(len(metaseq1)):
				if metaseq[i]==metaseq1[j]:
					metaseqoriginal.append(metaseq2[j])
					break
				
		print (metaseqoriginal)
		print (len(metaseqoriginal))

		#filetowrite = open("./resultsseq/"+ffile, "a")
		filetowrite = open("./resultsseqtemp/"+(10-len(str(filenumber)))*"0"+filenumber+".fasta", "w")
		handle = open(input_file, "r")
		
		for xx in SeqIO.parse(handle, "fasta"):
			#for i in range(len(metaseqoriginal)):
			if str(xx.name) in metaseqoriginal:
				#filetowrite.write(">"+str(xx.name)+" "+closesthittree[metaseq.index(str(xx.name))]+"\n")
				filetowrite.write(">"+str(xx.name)+"\n")
				filetowrite.write(str(xx.seq)+"\n")
				#break
		handle.close()
		filetowrite.close()

		#########################################################################################################
		
		
		filetowrite = open("resultsclusters"+filenumber+".txt", "w")
		
		towrite=fname+"\t"
		for i in range(len(clusters)):
			ccounter=0
			for j in range(len(closesthittree)):
				if closesthittree[j]==clusters[i]:
					ccounter+=1
					tccounter+=1
					
			towrite+=str(ccounter)+"\t"
					
		filetowrite.write(towrite[:-1]+"\n")
		filetowrite.close()
		
	if jj==0:
		filetowrite = open("resultsclusters"+filenumber+".txt", "w")
		
		towrite=fname+"\t"

		for i in range(len(clusters)):
			towrite+="0\t"
					
		filetowrite.write(towrite[:-1]+"\n")
		filetowrite.close()
    
	cmd = ["grep -c '>' del1_"+filenumber+".fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	nseq1=out.decode('ascii').strip()

	cmd = ["grep -c '>' del2_"+filenumber+".fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	nseq2=out.decode('ascii').strip()
	
	filetowrite = open("removedseq"+filenumber+".txt", "w")
	filetowrite.write(fname+"\t"+nseq1+"\t"+nseq2+"\t"+str(tccounter)+"\n")
	filetowrite.close()

	print()
	print (fname)
	print ("Old number of sequences: "+nseq1)
	print ("New number of sequences: "+nseq2)
	print ("Final number of sequences: "+str(tccounter))
	print()
	

	
from ruffus import *

def fnames():
	files=[]
	dirList=os.listdir(readspath)
	for fname in dirList:
		files.append(readspath+fname)

	parameters=[[0 for i in range(3)] for j in range(len(files))]
	for d1 in range(len(files)):
		parameters[d1][0]= files[d1]
		parameters[d1][1]= files[d1].replace(readspath, pathresults).replace(".fasta", ".newick")
		parameters[d1][2]= str(d1+1)

	for job_parameters in parameters:
		yield job_parameters

@files(fnames)
def parallel_task(input_file, output_file, i):
	print (i+"/"+nsequences+" "+input_file.replace(readspath,"").replace(".fasta",""))

	fname=input_file.replace(readspath,"")

	loop(input_file, output_file, i, aamodel, str(nproc), clusters, ffile)
	
pipeline_run([parallel_task], verbose=1, multiprocess=nmultiprocess)


dirList=os.listdir("./resultsseqtemp/")
ffilesm=[]
for ff in dirList:
	ffilesm.append(ff)

ffilesm.sort()

for j in range(len(ffilesm)):
	print ("cat ./resultsseqtemp/"+ffilesm[j]+" >> ./resultsseq/"+ffile)
	os.system("cat ./resultsseqtemp/"+ffilesm[j]+" >> ./resultsseq/"+ffile)

cmd = ["cat removedseq*.txt > removedseq"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm removedseq*.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mv removedseq removedseq.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["cat resultsclusters*.txt >> resultsclusters"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm resultsclusters*.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mv resultsclusters resultsclusters.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm temp*.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm del?_*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm *.job*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm query*.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm reference*.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm -r ./outdir*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm -r ./temp*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r ./resultsseqtemp"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()


print()
print ("I'm done!")
