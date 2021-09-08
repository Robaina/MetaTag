#!/usr/bin/env python
#

import os, subprocess, sys
from Bio import SeqIO
from Bio.Seq import Seq

global ppath, nsequences


###########################################################################################################################

ppath="/usr/gonzalez/prodigal/aMARgenomesPartialSept2020/"
ppath="/usr/gonzalez/prodigal/aMARgenomesCompleteSept2020/"
#ppath="/usr/gonzalez/prodigal/aTARA/"
#ppath="/home/gonzalez/hmmer/aAcidiferrobacteraceae/"
#ppath="/usr/gonzalez/paoli2021/a/"
#ppath="/home/gonzalez/hmmer/files/"
#ppath="/usr/gonzalez/metagenomes/refseq/"
ppath="/home/gonzalez/hmmer/envision/"
#ppath="/usr/gonzalez/metagenomes/MarPeptides/"

pathtonucleotides="/usr/gonzalez/prodigal/filesMARgenomesPartialSept2020/"
#pathtonucleotides="/usr/gonzalez/simo/acroporaassemblies/nucl/"
pathtonucleotides="/usr/gonzalez/prodigal/filesMARgenomesCompleteSept2020/"
#pathtonucleotides="/usr/gonzalez/prodigal/filesAcidiferrobacteraceae/"
#pathtonucleotides="/usr/gonzalez/paoli2021/files/"
#pathtonucleotides="/home/gonzalez/hmmer/files/"
#pathtonucleotides="/home/gonzalez/hmmer/genomesnt/"

lengthnucleotides=20000 #1E7

pfammodel="TIGR01287"  # NifH
evalue="1E-90"

if len(sys.argv)>1:
	pfammodel=sys.argv[1].strip()
	
	question="y"
	question1="y"
	question2="n"
	question3="n"
	question4="n"
	
	x=20
	
	

###########################################################################################################################

cmd = ["ls -1 "+ppath+" | wc -l"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (err.decode('ascii'))
nsequences=out.decode('ascii').strip()
print ("Number of sequences: "+nsequences)
'''
pfammodel="TIGR01326" # O-acetylhomoserine aminocarboxypropyltransferase/cysteine synthase, MetY
pfammodel="HmmbuildDmdA_Database_RBH.hmm"
pfammodel="TPPTonBFlavos.hmm"
pfammodel="TPPTonBGammas.hmm"
pfammodel="RecAFlavos.hmm"
pfammodel="Pfam-A.hmm"
pfammodel="TIGRFAMs_15.0_HMM"
pfammodel="TIGR01782"
pfammodel="PF18761" # heliorhodopsin
pfammodel="PF01036" # Bac_rhodopsin
pfammodel="rhodopsinflavos.hmm" # 1E-130
pfammodel="PF14064" # HmuY
pfammodel="TonBHeme.hmm"
pfammodel="TIGR01785" # TonB-hemin 
pfammodel="PF00593" # TonB dependent receptor
pfammodel="TIGR04056" # TonB-linked outer membrane protein, SusC/RagA family
pfammodel="TonBTIGR01783.hmm"
pfammodel="TIGR02012" # RecA
pfammodel="PF00154" # RecA
pfammodel="TIGR02066" # DsrB 
pfammodel="TIGR01115" # PufM
pfammodel="TIGR01157" # PufL
pfammodel="pfams.txt"
pfammodel="TIGR04183" # Por_Secre_tail
pfammodel="TIGR04183" # MetH methionine synthase
pfammodel="PF05694" # MtoX  doi: 10.1038/ismej.2017.148
pfammodel="TIGR01779" # BtuB
pfammodel="TPPTonB.hmm"
pfammodel="B12TonB.hmm"
pfammodel="PF00016" # rubisco large subunit
pfammodel="PF00485" # phosphoribulokinase 
pfammodel="PF00101" # RuBisCO_small
pfammodel="PF00365" # 6-phosphofructokinase, EMP
pfammodel="TIGR01182" # 2-keto-3-deoxygluconate-6-phosphate aldolase, ED
pfammodel="4_hydroxybutyryl_CoA_dehydratase.hmm" # Thaumarchaeota evalue="1E-200"
pfammodel="PF00374" # NiFeSe_Hases


pfammodel="AcuI.hmm"
pfammodel="Alma1.hmm"
pfammodel="DddK.hmm"
pfammodel="DLL-Cupin.hmm"
pfammodel="PF04392"
pfammodel="DmdA.hmm"
pfammodel="DmdB.hmm"
pfammodel="DmdC.hmm"
pfammodel="gammas1.hmm"
pfammodel="DddD.hmm"
pfammodel="DddProseos.hmm"
pfammodel="DddQ.hmm"
pfammodel="DddY.hmm"
pfammodel="DddL.hmm"
fammodel="DddW.hmm"
pfammodel="DsyBLabrenzia.hmm"
pfammodel="DSYB.hmm"
pfammodel="PF02028" # BCCT
pfammodel="PF04069" # OpuAC
pfammodel="TIGR01326" # MetY
pfammodel="MtoX.hmm"
pfammodel="TdaF.hmm"
pfammodel="TIGR04486" # SoxB
pfammodel="DsrA.hmm"
pfammodel="TIGR02416" # CODH

# http://tigrfams.jcvi.org/cgi-bin/Listing.cgi
pfammodel="TIGR01776" # 0 hits
pfammodel="TIGR01778" # TonB-copper 47 hits
pfammodel="TIGR01779" # TonB-B12 46 hits, mostly vibrio
pfammodel="TIGR01782" # TonB-Xanth-Caul 861 hits but mechanism is not known
pfammodel="TIGR01783" # TonB-siderophor 1235 hits
pfammodel="TIGR01786" # TonB-hemlactrns 271 

pfammodel="PF17850" # CysA
pfammodel="TIGR02034" # CysN
pfammodel="PF01507" # CysD
pfammodel="TIGR00434" # CysH
pfammodel="PF01747"  # Sulfate adenylyltransferase
pfammodel="TIGR02061"  # AprA
pfammodel="TIGR02060"  # AprB
pfammodel="TIGR01287"  # NifH
pfammodel="TIGR01580"  # NarG

evalue="1E-110" # TonBHeme.hmm
evalue="1E-200" # TPPTonBFlavos.hmm, TPPTonBGammas.hmm, RecAFlavos.hmm
evalue="1E-150" # TonBTIGR01783.hmm

evalue="1E-130" # DmdA
evalue="1E-100" # DmdB
evalue="1E-200" # AcuI
evalue="1E-70" # DddK
evalue="1E-60" # Alma1.hmm 
evalue="1E-80" # DddL
evalue="1E-200" # DddY.hmm, DddProseos.hmm, DddD.hmm
evalue="1E-60" # DddW
evalue="1E-50" # DddQ
evalue="1E-230" # gammas1.hmm
evalue="1E-200" # MtoX.hmm
evalue="1E-220" # DsrA.hmm
evalue="1E-110" # TdaF.hmm
evalue="1E-100" # TPPTonB.hmm
evalue="1E-60" # B12TonB.hmm
evalue="1E-100" # DsyBLabrenzia.hmm
evalue="1E-90" # DSYB.hmm
evalue="1E-200" # Tmm.hmm
'''
pfammodel=pfammodel.strip(); print (pfammodel); print ()

if len(sys.argv)==1:

	question = input("Should I start over and remove results folders? ")
	#question="yes"
	print (question)
	
	question=question.lower()
	question=question[0]
	
	question1 =input("Gathering score? ")
	#question1="no"
	question1=question1.lower()
	question1=question1[0]
	
	question2 =input("Retrieve contigs with hits? ")
	#question2="no"
	question2=question2.lower()
	question2=question2[0]
	
	question3 =input("Retrieve gene hits? ")
	#question3="yes"
	question3=question3.lower()
	question3=question3[0]
	
	question4 =input("Retrieve genome if there are hits? ")
	question4=question4.lower()
	question4=question4[0]

if question1=="y":
	gatheringscore=1
else:
	gatheringscore=0
	print ("Highest evalue: "+evalue)

print ("Wait!")

if question=="y":
	cmd = ["rm -r /usr/gonzalez/hmmer/results"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["mkdir /usr/gonzalez/hmmer/results"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm -r /usr/gonzalez/hmmer/results2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["mkdir /usr/gonzalez/hmmer/results2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm -r /usr/gonzalez/hmmer/resultspep"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["mkdir /usr/gonzalez/hmmer/resultspep"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r /usr/gonzalez/hmmer/resultsgenes"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
if question3=="y":
	cmd = ["mkdir /usr/gonzalez/hmmer/resultsgenes"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r ./resultsgenomes"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
if question4=="y":
	cmd = ["mkdir ./resultsgenomes"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r /usr/gonzalez/hmmer/resultsnucleotides"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
if question2=="y":
	cmd = ["mkdir /usr/gonzalez/hmmer/resultsnucleotides"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	#cmd = ["rm -r ./resultspos"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	#cmd = ["mkdir ./resultspos"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
else:
	dirList=os.listdir("/usr/gonzalez/hmmer/results/")
	for fname in dirList:
		if os.path.getsize("/usr/gonzalez/hmmer/results/"+fname)==0:
			cmd = ["rm /usr/gonzalez/hmmer/results/"+fname]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
#'''
cmd = ["rm *.hmm"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()	
cmd = ["rm *.h3?"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
#'''
if pfammodel=="TIGRFAMs_15.0_HMM":
	cmd = ["cp /usr/gonzalez/hmmer/TIGRFAMs_15.0_HMM.tar.gz ."]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["tar -zxvf TIGRFAMs_15.0_HMM.tar.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["cat *.HMM > TIGRFAMs_15.0_HMM"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["rm *.HMM"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))

elif pfammodel=="Pfam-A.hmm":
	print ("cp /usr/gonzalez/hmmer/Pfam-A.hmm .")
	cmd = ["cp /usr/gonzalez/hmmer/Pfam-A.hmm ."]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	
elif pfammodel.startswith("TIGR"):
	cmd = ["cp /usr/gonzalez/hmmer/TIGRFAMs_15.0_HMM.tar.gz ."]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["tar -zxvf TIGRFAMs_15.0_HMM.tar.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["mv "+pfammodel+".HMM "+pfammodel]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["rm *.HMM"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["grep 'DESC  ' "+pfammodel]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))       
	cmd = ["rm TIGRFAMs_15.0_HMM.tar.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

elif pfammodel.startswith("PF"):
	cmd = ["grep -m 1 '"+pfammodel+"' /usr/gonzalez/hmmer/Pfam-A.hmm"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	result=out.decode('ascii')[3:].strip()
	print ("hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+result+" > "+pfammodel)
	cmd = ["hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+result+" > "+pfammodel]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	cmd = ["grep -m 1 'DESC  ' "+pfammodel]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))

elif pfammodel=="pfams.txt":
	#'''
	filetowrite = open("pfams", "w")
	filetowrite.close()
	
	filetoread = open("pfams.txt", "r")
	for line in filetoread:
		print (line)
		line=line[:7].strip()
		if len(line)>0:
			if not(line.startswith("#")):
				cmd = ["grep -m 1 '"+line+"' /usr/gonzalez/hmmer/Pfam-A.hmm"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
				print (out.decode('ascii'), err.decode('ascii'))
				result=out.decode('ascii')[3:].strip()
				print ("hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+result+" >> pfams")
				cmd = ["hmmfetch /usr/gonzalez/hmmer/Pfam-A.hmm "+result+" >> pfams"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
				print (out.decode('ascii'), err.decode('ascii'))
				
		else:
			break

	filetoread.close()

	cmd = ["grep 'DESC  ' pfams"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii').strip())
	#'''
	pfammodel="pfams"
	
else:
	cmd = ["cp /usr/gonzalez/hmmer/"+pfammodel+" ."]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
#'''
if len(sys.argv)==1:
	ans=input("Do I keep going? ")
else:
	ans="yes"

if ans.lower()[0]=="n":
	cmd = ["rm "+pfammodel+".h3*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	if pfammodel[-4:]!=".hmm":
		cmd = ["rm "+pfammodel]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm TIGRFAMs_15.0_HMM"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	exit()

print ("hmmpress "+pfammodel)
cmd = ["hmmpress "+pfammodel]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print (out.decode('ascii'), err.decode('ascii'))
#'''
print("")
if len(sys.argv)==1:
	x = input("How many processors? ")
	#x=1
	x=int(x)
	if x>20:
		x=20
print("")

def retrievedesc(output_file):
	towrite=""
	filetoread = open("/usr/gonzalez/hmmer/results/"+output_file, "r")
	lblpeptides=[]
	mhmmacc=[]
	for line in filetoread:
		line=line.strip()
		if not(line.startswith("#")) and len(line)>0:
			seqid=line[line.find(" "):].strip(); seqid=seqid[seqid.find(" "):].strip(); seqid=seqid[:seqid.find(" ")].strip()
			lblpeptides.append(seqid)
			hmmacc=line[line.find(" "):].strip(); hmmacc=hmmacc[:hmmacc.find(" ")].strip()
			mhmmacc.append(hmmacc)
			hmmname=line[:line.find(" ")].strip()
			for i in range(18):
				line=line[line.find(" "):].strip()
			towrite+=seqid+chr(9)+hmmacc+chr(9)+hmmname+chr(9)+line+"\n"
	filetoread.close()
	
	nhits=0
	if len(towrite)>0:
		filetowrite = open("/usr/gonzalez/hmmer/results2/"+output_file, "w")
		filetowrite.write(towrite)
		filetowrite.close()
		
		filetowrite = open("/usr/gonzalez/hmmer/resultspep/"+output_file[:output_file.find(".")]+".fasta", "w")
		handle = open(ppath+output_file[:output_file.find(".")]+".fasta", "r")		
		for xx in SeqIO.parse(handle, "fasta"):
			if str(xx.name) in lblpeptides:
				nhits+=1
				filetowrite.write(">"+str(xx.description))
                
				towrite=""
				for i in range(len(mhmmacc)):
					towrite+=mhmmacc[i]+" "
                       
				filetowrite.write(" "+towrite.strip()+"\n")
                
				filetowrite.write(str(xx.seq)+"\n")
		handle.close()
		filetowrite.close()
		
	return nhits

if question=="y":		
	from ruffus import *

	def fnames():
		files=[]
		dirList=os.listdir(ppath)
		for fname in dirList:
			files.append(ppath+fname)

		parameters=[[0 for i in range(3)] for j in range(len(files))]
		for d1 in range(len(files)):
			parameters[d1][0]= files[d1]
			resultfile=files[d1].replace(ppath,""); resultfile=resultfile[:resultfile.find(".")]+".txt"
			parameters[d1][1]= resultfile
			parameters[d1][2]= str(d1+1)

		for job_parameters in parameters:
			yield job_parameters

	@files(fnames)
	def parallel_task(input_file, output_file, i):
		print (i+"/"+nsequences+" "+output_file[:output_file.find(".")])

		if os.path.getsize(input_file)>0:
			
			if gatheringscore==1:
				if int(nsequences)==1:
					print ("hmmscan --tblout /usr/gonzalez/hmmer/results/"+output_file+" -o /dev/null --cut_ga "+pfammodel+" "+input_file)
					cmd = ["hmmscan --tblout /usr/gonzalez/hmmer/results/"+output_file+" -o /dev/null --cut_ga "+pfammodel+" "+input_file]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
					print (out.decode('ascii'), err.decode('ascii'))
				else:
					print ("hmmscan --tblout /usr/gonzalez/hmmer/results/"+output_file+" -o /dev/null --cut_ga --cpu 1 "+pfammodel+" "+input_file)
					cmd = ["hmmscan --tblout /usr/gonzalez/hmmer/results/"+output_file+" -o /dev/null --cut_ga --cpu 1 "+pfammodel+" "+input_file]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
					print (out.decode('ascii'), err.decode('ascii'))
			else:
				if int(nsequences)==1:
					print ("hmmscan --tblout /usr/gonzalez/hmmer/results/"+output_file+" -o /dev/null -E "+evalue+" "+pfammodel+" "+input_file)
					cmd = ["hmmscan --tblout /usr/gonzalez/hmmer/results/"+output_file+" -o /dev/null -E "+evalue+" "+pfammodel+" "+input_file]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
					print (out.decode('ascii'), err.decode('ascii'))			
				else:
					print ("hmmscan --tblout /usr/gonzalez/hmmer/results/"+output_file+" -o /dev/null -E "+evalue+" --cpu 1 "+pfammodel+" "+input_file)
					cmd = ["hmmscan --tblout /usr/gonzalez/hmmer/results/"+output_file+" -o /dev/null -E "+evalue+" --cpu 1 "+pfammodel+" "+input_file]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
					print (out.decode('ascii'), err.decode('ascii'))

			nhits=retrievedesc(output_file)
			print (i+"/"+nsequences+" "+output_file[:output_file.find(".")])
			print ("No. hits: "+str(nhits)); print()

	pipeline_run([parallel_task], verbose=1, multiprocess=x)

filetowrite = open("/usr/gonzalez/hmmer/results.fasta", "w")
filetowriteother = open("/usr/gonzalez/hmmer/resultsother.fasta", "w")

filetowritetemp = open("temp.txt", "w")
filetowritetemp.write("Genome\tNo. peptides\t"+pfammodel+"\n")

dirList=os.listdir("/usr/gonzalez/hmmer/results/")
ffiles=[]
for fname in dirList:
	ffiles.append(fname[:fname.find(".")])
		
ffiles.sort()

genomehits=0
for x in range(len(ffiles)):
	print (str(x+1)+"/"+str(len(ffiles))+"     "+ffiles[x])
	filetowritetemp.write(ffiles[x]+"\t")
	filetoread = open("/usr/gonzalez/hmmer/results/"+ffiles[x]+".txt", "r")
	hits=[]
	hitstigrfam=[]
	evalue=[]
	for line in filetoread:
		if not(line.startswith("#")):
			hitstigrfam.append(line[:line.find(" ")].strip())
			line=line[line.find(" "):].strip()
			line=line[line.find(" "):].strip()
			hits.append(line[:line.find(" ")].strip())
			line=line[line.find(" "):].strip()
			line=line[line.find(" "):].strip()
			evalue.append(float(line[:line.find(" ")].strip()))
	filetoread.close()
	
	if len(hits)>0:
		genomehits+=1

	handle = open(ppath+ffiles[x]+".fasta", "r")
	j=0; npeptides=0
	for xx in SeqIO.parse(handle, "fasta"):
		act=0
		if len(hits)>0:
			if str(xx.name) in hits:
				filetowrite.write(">"+str(xx.description)+"\n")
				#filetowrite.write(">"+str(xx.description)+"_"+pfammodel+"\n")
				#filetowrite.write(">"+str(xx.description)+"_"+str(evalue[j])+"\n")
				filetowrite.write(str(xx.seq)+"\n")
				j+=1
				act=1
			#else:    
			#	filetowriteother.write(">"+str(xx.description)+"\n")
			#	filetowriteother.write(str(xx.seq)+"\n")                
		#else:
		if act==0 and len(sys.argv)==1:
			filetowriteother.write(">"+str(xx.description)+"\n")
			filetowriteother.write(str(xx.seq)+"\n")
		npeptides+=1
	handle.close()
			
	for i in range(len(hits)):
		if question3=="n":
			break
		filetowrite2 = open("/usr/gonzalez/hmmer/resultsgenes/"+ffiles[x]+".fasta", "a")
		a=hits[i][2+hits[i].find("__"):]
		contig=int(a.split("_")[0])

		pos1=int(a.split("_")[2])-1
		pos2=int(a.split("_")[3])
		direction=a.split("_")[4]
		
		handle = open(pathtonucleotides+ffiles[x]+".fasta", "r")
		act=0
		for xx in SeqIO.parse(handle, "fasta"):
			if int(str(xx.name))==contig:
				filetowrite2.write(">"+hits[i]+"\n")
				towrite=str(xx.seq)[pos1:pos2]
				if direction=="neg":
					towrite=str(Seq(str(xx.seq)[pos1:pos2]).reverse_complement())

				filetowrite2.write(towrite+"\n")
				act=1
				break
		handle.close()
		filetowrite2.close()
	
		if act==0:
			print ("I couldn't find the contig!")
			print (hits[i])
			print (ffiles[x])
			print (pathtonucleotides+ffiles[x]+".fasta")
			exit()
	
		filetowrite2.close()
		
	if question4=="y" and len(hits)>0:
		print ("cp "+pathtonucleotides+ffiles[x]+".fasta ./resultsgenomes")
		cmd = ["cp "+pathtonucleotides+ffiles[x]+".fasta ./resultsgenomes"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
		print (out.decode('ascii'), err.decode('ascii'))
	
	for i in range(len(hits)):
		if question2=="n":
			break
		filetowrite2 = open("/usr/gonzalez/hmmer/resultsnucleotides/"+ffiles[x]+"_"+str(i+1)+".fasta", "w")
		a=hits[i][hits[i].find("__")+2:]
		
		#print ("a: "+a)

		m=a.split("_")
		
		#print (m)
		
		#print (len(a)-len("_".join(m[:-4])))
		
		#print (len("_".join(m[-4:])))
		
		contig=int(a[:len(a)-len("_".join(m[-4:]))-1])
		
		#contig+="_"+"_".join(m[-4:])
		
		#print (contig)
        
		
		#exit()

		
		#contig=str(int(a[:a.find("_")]))
		

		

		
		#a=a[1+a.find("_"):]; a=a[1+a.find("_"):]; pos=int(a[:a.find("_")])
		
		pos=int(m[-2])
		
		pos1=pos-lengthnucleotides
		if pos1<0:
			pos1=1
		pos2=int(pos+lengthnucleotides)

		handle = open(pathtonucleotides+ffiles[x]+".fasta", "r")
		act=0
		for xx in SeqIO.parse(handle, "fasta"):
			if str(xx.name)==str(contig):
				filetowrite2.write(">"+str(xx.description)+"_"+str(i+1)+"\n")
				filetowrite2.write(str(xx.seq)[pos1:pos2]+"\n")
				act=1
				break
		handle.close()
		filetowrite2.close()
		if act==0:
			print (contig)
			print ("Problem!")
			exit()
	
	filetowritetemp.write(str(npeptides)+"\t"+str(len(hits))+"\n")
	print ("No. peptides: "+str(npeptides))
	print ("No. hits: "+str(len(hits)))

filetowritetemp.close()
filetowriteother.close()
filetowrite.close()

dirList=os.listdir("/usr/gonzalez/hmmer/results2/")
ffiles=[]
for fname in dirList:
	ffiles.append(fname)
		
ffiles.sort()

filetowrite = open("results.txt", "w")

for x in range(len(ffiles)):
	s = open("/usr/gonzalez/hmmer/results2/"+ffiles[x], 'r').read()
	filetowrite.write(s+20*"-"+"\n")
filetowrite.close()

cmd = ["grep -c '^>' /usr/gonzalez/hmmer/results.fasta"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
print(); print ("HMM model: "+pfammodel); print (str(int(out.decode('ascii')))+" hits.")
print (str(genomehits)+" sequences with at least one hit.")
'''
filetowrite = open("resultspfams.txt", "w")
filetoread = open("results.txt", "r")
pfams=[]
for line in filetoread:
	m=line.split()
	if len(m)>1:
		a=m[1][:m[1].find(".")]
		if not(a in pfams):
			pfams.append(a)
			filetowrite.write(a+"\n")
filetoread.close()
filetowrite.close()
'''
#'''
cmd = ["rm *.h3?"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm pfams"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm "+pfammodel]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm *.tar.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
#'''
cmd = ["cp /usr/gonzalez/hmmer/results.fasta ."]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

os.system("python /home/gonzalez/send_email.py")
print ()
print ("I'm done!")
exit()
