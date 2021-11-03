from Bio.SeqIO.FastaIO import SimpleFastaParser

newseqs = {}
equivalent = {}
# parse fasta file with the low-level SimpleFastaParser, reads it as a tuple
with open("marref_sample10pc.fna") as sequences:
	for k, seq in enumerate(SimpleFastaParser(sequences)):
		newseqs[k]=seq[1]
		equivalent[k]=seq[0]
		

ofile = open("marsample.fna", "w")
for i in newseqs.keys():
	ofile.write(">{}\n{}\n".format(i, newseqs[i]))
ofile.close()


ofile = open("marsample_names_equivalence.txt", "w")
for i in equivalent.keys():
	ofile.write("{}\t{}\n".format(i, equivalent[i]))
ofile.close()

