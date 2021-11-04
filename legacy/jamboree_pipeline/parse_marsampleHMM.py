from Bio import SearchIO
import csv

# Read the output table of database search with TIGRFAMs HMMs for nosZ
hmm_qresult = SearchIO.parse('marsample_nosZ_tigrfam.hmm',  'hmmsearch3-domtab')

# Filter table to make a list of hits IDs
hit_ids = []
for qresult in hmm_qresult:
    for i in range(len(qresult)):
        hit_ids.append(qresult[i].id)

# Write list to a text file
file = open('hits_nosZ_marsample.txt', 'w')
for index in range(len(hit_ids)):
    file.write(str(hit_ids[index]) + "\n")
file.close()

