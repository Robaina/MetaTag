#!/usr/bin/env python
# 


"""
es sencillo

es solo complicado porque lo hice para que use pfam, tigrfam o un hmm que hago

pfam1 y pfam2 son los dos hmm que busca

gap es el número de genes entre ellos

lengthnucleotides es el fragmento del contig que guarda con esos dos genes

que busca.

cuando están a menos de 10 genes entre ellos

primero busca con un hmm y si lo encuentra entonces busca con el otro

si encuentra los dos entonces mira cuántos genes hay entre ellos

lo da por bueno cuando hay 10 o menos genes entre los hits para cada hmm

entiendes?

guarda los peptidos y la zona del contig donde están esos dos genes
"""

import os, re, subprocess
from Bio import SeqIO

regex = re.compile('[^a-zA-Z0-9.]')

nproc=1 # Goes slower if greater than 1.
gap=10

path="/home/gonzalez/hmmer/aAcidiferrobacteraceae/"
#path="/usr/gonzalez/prodigal/aMARgenomesCompleteSept2020/"
#path="/usr/gonzalez/prodigal/aMARgenomesPartialSept2020/"

lengthnucleotides=40000

pathtonucleotides="/home/gonzalez/hmmer/filesAcidiferrobacteraceae/"
#pathtonucleotides="/usr/gonzalez/prodigal/filesMARgenomesCompleteSept2020/"
#pathtonucleotides="/usr/gonzalez/prodigal/filesMARgenomesPartialSept2020/"


#pfam1="HmmbuildDmdA_Database_RBH.hmm" #1E-130
#pfam1="DmdB.hmm"
#pfam1="PF14064" # HmuY
#pfam1="PF00177" # Ribosomal_S7
#pfam1="PF00338" # Ribosomal_S10
pfam1="PF07980" # SusD_RagB
#pfam1="PF12741" # SusD-like
#pfam1="PF14322" # SusD-like_3
#pfam1="PF12771" # SusD-like_2
pfam1="pfams.txt"
pfam1="DmdA.hmm" #1E-130
pfam1="TIGR04486" # thiosulfohydrolase SoxB
#pfam1="TIGR03860"
#pfam1="sulfotransferases.hmm"

#pfam2="HmmbuildDddP.hmm" # -E 1E-90
#pfam2="PF00378" # Enoyl-CoA hydratase/isomerase
#pfam2="PF02770" # Acyl-CoA dehydrogenase, middle domain
#pfam2="PF00107" # Zinc-binding dehydrogenase
pfam2="DmdC.hmm" # All pfams in DmdC
#pfam2="PF00501" # DmdB
#pfam2="PF00441" # DmdC
#pfam2="TIGR01186" # glycine betaine/L-proline transport ATP binding subunit
#pfam2="PF13442" # SoxD and some other things
#pfam2="PF00576" # Hues/Transthyretin family
#pfam2="TIGR01326" # O-acetylhomoserine aminocarboxypropyltransferase/cysteine synthase
#pfam2="PF02737" # 3-hydroxyacyl-CoA dehydrogenase, NAD binding domain
#pfam2="TIGR01930" # acetyl-CoA C-acyltransferase
#pfam2="DmdB.hmm"
#pfam2="TIGR01326" # O-acetylhomoserine aminocarboxypropyltransferase/cysteine synthase
#pfam2="PF00593" # TonB_dep_Rec
pfam2="DsrA.hmm"
#pfam2="TIGR02823" # AcuI, putative quinone oxidoreductase, YhdH/YhfP family
#pfam2="PF04069" # OpuAC
#pfam2="PF02028" # BCCT
#pfam2="PF00732" # DddA

evalue1="-E 1E-130"
evalue1="--cut_ga"
evalue2="-E 1E-220"

cmd = ["rm *.h3?"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm -r resultsnucleotides"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir resultsnucleotides"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm pfam1"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm pfam2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

act=0
if pfam1.startswith("PF"):
	print ("cat /usr/gonzalez/hmmer/Pfam-A.hmm | grep '"+pfam1+"'")
	cmd = ["cat /usr/gonzalez/hmmer/Pfam-A.hmm | grep '"+pfam1+"'"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	out=out.decode('ascii')
	pfamtemp= regex.sub("", out[6:])
	print ("hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+pfamtemp+" > pfam1")
	cmd = ["hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+pfamtemp+" > pfam1"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	act=1

if pfam1.startswith("TIGR"):
	cmd = ["cp /usr/gonzalez/hmmer/TIGRFAMs_15.0_HMM.tar.gz ."]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["tar -zxvf TIGRFAMs_15.0_HMM.tar.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["mv "+ pfam1 +".HMM pfam1"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["rm *.HMM"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["rm TIGRFAMs_15.0_HMM.tar.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	act=1
	
if pfam1=="pfams.txt":
	filetowrite = open("pfam2", "w")
	filetowrite.close()
	
	filetoread = open("pfams.txt", "r")
	for line in filetoread:
		line=line.strip()
		if len(line)>0:
			line=line[:line.find(" ")]
			cmd = ["grep -m 1 '"+line+"' /usr/gonzalez/hmmer/Pfam-A.hmm"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
			print (out.decode('ascii'), err.decode('ascii'))
			result=out.decode('ascii')[3:].strip()
			print ("hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+result+" >> pfam1")
			cmd = ["hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+result+" >> pfam1"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
			print (out.decode('ascii'), err.decode('ascii'))
			
	cmd = ["grep 'DESC  ' pfams"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii').strip())
	filetoread.close()
	act=1

if act==0:
	print ("cp /usr/gonzalez/hmmer/"+pfam1+" pfam1")
	cmd = ["cp /usr/gonzalez/hmmer/"+pfam1+" pfam1"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

act=0
if pfam2.startswith("PF"):
	print ("cat /usr/gonzalez/hmmer/Pfam-A.hmm | grep '"+pfam2+"'")
	cmd = ["cat /usr/gonzalez/hmmer/Pfam-A.hmm | grep '"+pfam2+"'"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	out=out.decode('ascii')
	pfamtemp= regex.sub("", out[6:])
	print ("hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+pfamtemp+" > pfam2")
	cmd = ["hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+pfamtemp+" > pfam2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	act=1
	
if pfam2.startswith("TIGR"):
	cmd = ["cp /usr/gonzalez/hmmer/TIGRFAMs_15.0_HMM.tar.gz ."]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["tar -zxvf TIGRFAMs_15.0_HMM.tar.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["mv "+ pfam2 +".HMM pfam2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["rm *.HMM"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["rm TIGRFAMs_15.0_HMM.tar.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	act=1
	
if pfam2=="pfams.txt":
	filetowrite = open("pfam2", "w")
	filetowrite.close()
	
	filetoread = open("pfams.txt", "r")
	for line in filetoread:
		line=line.strip()
		if len(line)>0:
			line=line[:line.find(" ")]
			cmd = ["grep -m 1 '"+line+"' /usr/gonzalez/hmmer/Pfam-A.hmm"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
			print (out.decode('ascii'), err.decode('ascii'))
			result=out.decode('ascii')[3:].strip()
			print ("hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+result+" >> pfam2")
			cmd = ["hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+result+" >> pfam2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
			print (out.decode('ascii'), err.decode('ascii'))
			
	cmd = ["grep 'DESC  ' pfams"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii').strip())
	filetoread.close()
	act=1

if act==0:
	print ("cp /usr/gonzalez/hmmer/"+pfam2+" pfam2")
	cmd = ["cp /usr/gonzalez/hmmer/"+pfam2+" pfam2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

print ("hmmpress pfam1")
cmd = ["hmmpress pfam1"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))

print ("hmmpress pfam2")
cmd = ["hmmpress pfam2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))

filetowrite = open("DmdA"+pfam1.replace(".","")+pfam2.replace(".","")+".fasta", "w")
filetowrite2 = open("DmdAnext"+pfam1.replace(".","")+pfam2.replace(".","")+".fasta", "w")
filetowrite4 = open("DmdAonly"+pfam1.replace(".","")+pfam2.replace(".","")+".fasta", "w")

filetowrite5 = open("temp.txt", "w")
filetowrite5.close()

def loop(path, fname, evalue1, evalue2):
	print ("hmmscan --tblout results1.txt "+evalue1+" -o /dev/null --cpu "+str(nproc)+" pfam1 FASTA1.txt")
	cmd = ["hmmscan --tblout results1.txt "+evalue1+" -o /dev/null --cpu "+str(nproc)+" pfam1 FASTA1.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	
	nhits=0
	newhits1=0
	filetoread1 = open("results1.txt", "r")
	for line in filetoread1:
		line=line.strip()
		if len(line)>0 and not(line.startswith("#")):
			newhits1+=1
			print (line)
			evalue11=line[line.find(" "):].strip(); evalue11=evalue11[evalue11.find(" "):].strip(); evalue11=evalue11[evalue11.find(" "):].strip(); evalue11=evalue11[evalue11.find(" "):].strip()
			evalue11=evalue11[:evalue11.find(" ")].strip()
			print (evalue11)

			act=0
			pep1=line[line.find(" "):].strip()
			pep1= pep1[pep1.find(" "):].strip()
			pep1= pep1[:pep1.find(" ")].strip()
			
			contig1=pep1[2+pep1.find("__"):]
			pos1=contig1[1+contig1.find("_"):]
			pos1=int(pos1[:pos1.find("_")])

			contig1=str(int(contig1[:contig1.find("_")]))

			handle = open(path+fname, "r")
			filetowrite3 = open("FASTA2.txt", "w")
			for x in SeqIO.parse(handle, "fasta"):
				contig2=str(x.name)[2+str(x.name).find("__"):]
				contig2=str(int(contig2[:contig2.find("_")]))
				if contig2==contig1:
					filetowrite3.write(">"+str(x.name)+"\n")
					filetowrite3.write(str(x.seq)+"\n")
			filetowrite3.close()
			handle.close()
			
			print ("hmmscan --tblout results2.txt "+evalue2+" -o /dev/null --cpu "+str(nproc)+" pfam2 FASTA2.txt")
			cmd = ["hmmscan --tblout results2.txt "+evalue2+" -o /dev/null --cpu "+str(nproc)+" pfam2 FASTA2.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
			print (out.decode('ascii'), err.decode('ascii'))
			
			filetoread2 = open("results2.txt", "r")
			for line in filetoread2:
				line=line.strip()
				if len(line)>0 and not(line.startswith("#")):
					
					print (line)
					
					evalue21=line[line.find(" "):].strip(); evalue21=evalue21[evalue21.find(" "):].strip(); evalue21=evalue21[evalue21.find(" "):].strip(); evalue21=evalue21[evalue21.find(" "):].strip()
					evalue21=evalue21[:evalue21.find(" ")].strip()
					print (evalue21)
					
					pep2=line[line.find(" "):].strip()
					pep2= pep2[pep2.find(" "):].strip()
					pep2= pep2[:pep2.find(" ")].strip()
					print (pep2)
					
					contig2=pep2[2+pep2.find("__"):]
					pos2=contig2[1+contig2.find("_"):]
					pos2=int(pos2[:pos2.find("_")])
					print (pos2)
					
					if abs(pos1-pos2)<=gap:
						nhits+=1; act=1
						handle = open(path+fname, "r")
						for x in SeqIO.parse(handle, "fasta"):
							if str(x.name)==pep1:
								#filetowrite.write(">"+pep1+"\n")
								#filetowrite.write(">"+pep1+"_next_"+"_"+str(abs(pos1-pos2))+"_"+str(evalue11)+"\n")
								filetowrite.write(">"+pep1+"_next_"+str(abs(pos1-pos2))+"\n")
								filetowrite.write(str(x.seq)+"\n")
								break
						handle.close()
						
						handle = open(path+fname, "r")
						for x in SeqIO.parse(handle, "fasta"):
							if str(x.name)==pep2.strip():
								#filetowrite2.write(">"+pep2+"_next_"+"_"+str(abs(pos1-pos2))+"\n")
								#filetowrite2.write(">"+pep2+"\n")
								#filetowrite2.write(">"+pep2+"_next_"+"_"+str(evalue21)+"\n")
								#filetowrite2.write(">"+pep2+"_next_"+"_"+str(abs(pos1-pos2))+"_"+str(evalue21)+"\n")
								filetowrite2.write(">"+pep2+"_next_"+str(abs(pos1-pos2))+"\n")
								filetowrite2.write(str(x.seq)+"\n")
								break
						handle.close()

						print (pathtonucleotides+fname)
						handle = open(pathtonucleotides+fname, "r")
						
						a=pep2[pep2.find("__")+2:]
						a=a[1+a.find("_"):]
						a=a[1+a.find("_"):]
						posnuc=int(a[:a.find("_")])
						
						
						print (posnuc)
						posnuc1=posnuc-lengthnucleotides
						if posnuc1<0:
							posnuc1=1
						posnuc2=int(posnuc+lengthnucleotides)

						for x in SeqIO.parse(handle, "fasta"):
							if str(x.name)==str(contig1):
								filetowrite3 = open("./resultsnucleotides/"+fname[:fname.find(".")]+"_contig_"+str(contig1)+".fasta", "w")
								filetowrite3.write(">"+str(x.name)+"\n")
								filetowrite3.write(str(x.seq)[posnuc1:posnuc2]+"\n")
								filetowrite3.close()
								break
						handle.close()
			filetoread2.close()
			if act==0:
				print (pep1)
				handle = open(path+fname, "r")
				for x in SeqIO.parse(handle, "fasta"):
					if str(x.name)==pep1:
						filetowrite4.write(">"+pep1+"_alone\n")
						filetowrite4.write(str(x.seq)+"\n")
						break
				handle.close()
	filetoread1.close()
	
	filetowrite5 = open("temp.txt", "a")
	filetowrite5.write(fname+"\t"+str(nhits)+"\n")
	filetowrite5.close()
	
	cmd = ["rm results?.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

	return (newhits1)

dirList=os.listdir(path)
jj=0; ntotalhits1=0
for fname in dirList:
	jj+=1
	print(); print (str(jj)+"   "+fname)
	cmd = ["cp "+path+fname+" FASTA1.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	newhits1=loop(path, fname, evalue1, evalue2)
	ntotalhits1+=newhits1
	
print ("Total hits first model: "+str(ntotalhits1))

filetowrite4.close()
filetowrite2.close()
filetowrite.close()

cmd = ["rm *.h3?"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm FASTA?.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm pfam?"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["grep -c '>' DmdA"+pfam1.replace(".","")+pfam2.replace(".","")+".fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
out=out.decode('ascii')
print (out.strip()+" tandem hits.")

print() 
print() 	
print ("I'm done!")
	
exit()
