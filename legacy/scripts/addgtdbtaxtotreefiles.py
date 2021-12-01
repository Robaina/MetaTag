#!/usr/bin/env python
# 

import re, csv
from Bio import SeqIO

regex = re.compile('[^a-zA-Z0-9]')
regexletters = re.compile('[^A-Z]')

handle = open("sequencesLongLabels.fasta", "r")
genome=[]; longlbl=[]; longerlbl=[]
for xx in SeqIO.parse(handle, "fasta"):
	for y in range(len(str(xx.name))):
		if ord(str(xx.name)[y])<48 or ord(str(xx.name)[y])>57:
			break
	y+=1
	b=str(xx.name)[y:]
	b=b[:b.find("__")]

	genome.append(b)
	longlbl.append(regex.sub("", str(xx.description)))
	
	longerlbl.append("")
	
handle.close()

for i in range(len(genome)):
	print (genome[i])
	
	act=0
	# MARDB
	if genome[i].startswith("complete"):
		countlines=0
		with open("/usr/gonzalez/cleangenomes/CurrentComplete.tsv") as file:
			tsv_file = csv.reader(file, delimiter="\t")
			for line in tsv_file:
				countlines+=1
				if countlines==1:
					pos=line.index("taxon_lineage_names")
				
				for j in range(len(line)):
					if regex.sub("", line[j])==regex.sub("", genome[i][8:]):
						longerlbl[i]=line[pos]
						print (line[pos])
						act=1
						break
				if act==1:
					break

	elif regexletters.sub("", genome[i].split("_")[-1])=="MMP":
		countlines=0
		with open("/usr/gonzalez/cleangenomes/CurrentPartial.tsv") as file:
			tsv_file = csv.reader(file, delimiter="\t")
			for line in tsv_file:
				countlines+=1
				if countlines==1:
					pos=line.index("mmp_ID")
					postax=line.index("taxon_lineage_names")
					
				if genome[i].split("_")[-1]==line[pos]:
					longerlbl[i]=line[postax]
					print (longerlbl[i])
					act=1
					break
					
	elif genome[i].split("_")[0]=="TARA" or genome[i].split("_")[0]=="MALA" or genome[i].split("_")[0]=="BGEO" or genome[i].split("_")[0]=="BATS" or genome[i].split("_")[0]=="GORG" or genome[i].split("_")[0]=="HOTS":
		countlines=0
		with open("/usr/gonzalez/cleangenomes/Browsethegenomecollection.csv") as file:
			tsv_file = csv.reader(file, delimiter=",")
			for line in tsv_file:
				countlines+=1
				if countlines==1:
					pos=line.index("GTDB Taxonomy")
					
				if genome[i]==line[0]:
					longerlbl[i]=line[pos]
					print (longerlbl[i])
					act=1
					break
		
	elif genome[i].startswith("OceanDNA"):
		countlines=0
		with open("/usr/gonzalez/cleangenomes/S3.tsv") as file:
			tsv_file = csv.reader(file, delimiter="\t")
			for line in tsv_file:
				countlines+=1
				if countlines==1:
					pos=line.index("gtdb_classification")
					
				if genome[i]==line[0].replace("-", "_"):
					longerlbl[i]=line[pos]
					print (longerlbl[i])
					act=1
					break

	elif act==0:
		countlines=0
		with open("/usr/gonzalez/cleangenomes/CurrentPartial.tsv") as file:
			tsv_file = csv.reader(file, delimiter="\t")
			for line in tsv_file:
				countlines+=1
				if countlines==1:
					pos=line.index("full_scientific_name")
					postax=line.index("taxon_lineage_names")
					
				if countlines>1:
					if regex.sub("", genome[i])==regex.sub("", line[pos]):
						longerlbl[i]=line[postax]
						print (longerlbl[i])
						act=1
						break		

	else:
		if act==0:
			print ("Not found!")
			#exit()

tree2 = open("0newtreelonger.txt", 'r').read()
print (tree2)

from Bio import Phylo

tree = Phylo.read("0newtreelonger.txt", "newick")
for leaf in tree.get_terminals(): 
	print (leaf.name)

	if regex.sub("", str(leaf.name)) in longlbl:
		a=regex.sub("", str(leaf.name))
		#print (a)
		#print (longerlbl[longlbl.index(a)])
		b=longerlbl[longlbl.index(a)].strip()
		b=regex.sub("_", b)
		while b.count("__")>0:
			b=b.replace("__", "_")
		#print (b)
		
		b=b.replace("_", " ")
		b=b.strip().replace(" ", "_")
		
		tree2=tree2.replace(str(leaf.name), str(leaf.name)+" "+b)
	
print (tree2)

filetowrite = open("0newtreelonger2.txt", "w")
filetowrite.write(tree2+"\n")
filetowrite.close()

print ("No problem!")
exit()