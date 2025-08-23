# sequence-tools
A repository for custom tools for sequence and tree data manipulation.

# python
- stripSeqs.py

parse the output from [meshclust](https://github.com/BioinformaticsToolsmith/MeShClust) and produce a fasta file of reads for each cluster

- grampaRename.py

change the tip labels of gene trees so that they are suitable for use in [GRAMPA](https://github.com/gwct/grampa)

- flyeLogParser.py

Parse the output of several [flye](https://github.com/mikolmogorov/Flye) assembly runs and select single contigs

- resolveMultiContigs.py

Parse the output of flye and flyeLogParser.py, and attempt to pick a single contig for assembilies with multiple. Selection is based on the quality of mapping against some reference

- plotVCF.py

Take in the output of summariseVCF.sh and plot a histogram of SNP frequencies. Filtering based on read depths

- ete_multi_root.py

An [ete4](https://github.com/etetoolkit/ete4) based script that re-roots all newick trees in a file onto a single taxa

# shell scripts
- summariseVCF.sh

Take in a VCF file and report the read depth and allele ratios. Plot frequencies using plotVCF.py using various filters.

- splitIPRLoci.sh

Take in the .loci output from [ipyrad](https://github.com/dereneaton/ipyrad) and convert this into a bunch of .fasta files.
