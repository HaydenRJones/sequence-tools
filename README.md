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
