## Obtaining reference tree

1. Preprocess database. Only run once to remove duplicates and sequences containing illegal symbols:
```
python3 src/preprocess.py --help
```

2. Make reference (protein-specific) peptide database. Run to filter origial database by given tigrfam or pfam. Also reduces redundancy of peptide database and allows filtering database sequences by minimum and/or maximum sequence length
```
python3 src/makedatabase.py --help
```

3. Perform multiple sequence alignment (MSA) on reference database and infer reference tree
```
python3 src/buildtree.py --help
```
## Curating reference tree

* Remove tree and reference msa outliers
```
python3 src/removetreeoutliers.py --help
```
* Manual curation and classification of clusters

## Placement of query short reads on reference tree

1. Preprocess query sequences. Remove illegal symbols from nucleotide or peptide sequence (if already translated). Can also translate nucleotide sequence with prodigal if not already translated.
```
python3 src/preprocess.py --help
```
2. Place sequences onto tree. Run either papara or hmmalign to align query sequences to reference MSA. Call epa-ng to place sequences onto tree. Call gappa to make new tree newick file with placements.
```
python3 src/placesequences.py --help
```
3. Assign taxonomy and function to placed sequences. Assign taxonomy and function to placed query sequences.
```
python3 src/labelplacements.py --help
```
4. Quantify placed queries. Filter placed sequences according to cluster/placement quality score and quantify remaining placed sequences according to function and taxonomy
```
python3 src/countplacements.py --help
```