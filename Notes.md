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


# TODO

1. Check/implement sequence preprocessing: e.g., removal of short sequences, trimming...

2. Implement placement

3. Think about automatizing short read labelling based on (already labelled) tree placement