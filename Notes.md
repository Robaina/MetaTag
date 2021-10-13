# Tree reconstruction

1. iqtree took 3 hours to ake a tree while fasttree took 4 minutes with the same input
   I can handle output files within python, moving them to specified location.

2. Both iqtree and fasttree produce unrooted trees.

3. May be needed to clean sequence labels of punctuation marks, since 
   software to do placement seems to complain (https://www.polarmicrobes.org/phylogenetic-placement-revisited/)

4. Trees are naturally messy because we are using many sequences. Ways to improve this:
   4.1 Remove uninformative (very similar) sequences: cd-hit
   4.2 Run trimal
   4.3 Remove outlier branches: treeshrink or similar

5. Reference peptide database should contain less than 1000 sequences (600 - 800) so manual tree analysis is manageable.
   CD-hit may be used to reduce database. Is  there a way to reduce it to given maximum number of sequences already implemented in CD-hit?
   Running Cd-hit with default params reduced nxr database from 1268 seqs to 510 seqs. That's a lot.


# TODO

1. Think about automatizing short read labelling based on (already labelled) tree placement

2. Python's pathlib much better than os.path to handle paths and files. See this [post](https://medium.com/@ageitgey/python-3-quick-tip-the-easy-way-to-deal-with-file-paths-on-windows-mac-and-linux-11a072b58d5f)

3. Contree used by Jose, tree best model (sometimes yeast or insect models appear, seems wrong).This is because iqtree is testing different substitution models suring optimization, among which those appear.

# Code meeting notes:

1. Preprocessing:
   1.1. Jose clean.py made for DNA seqs, and removes illegal symbols but keeps sequence. This may introduce artifacts (i.e., deletions) in the sequences. I remove these sequences from dataset, also made function to check peptide sequences.
   1.2. Parallezation of tasks over input files no longer depends on ruffus. Implemented with 
   multiprocessing Pool instead. Addtional arguments may be passed (**kwargs) 
   1.3 Working directly from peptides or translating with prodigal?
   