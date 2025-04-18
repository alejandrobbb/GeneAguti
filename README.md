This script performs a comparative analysis of the ASIP gene (Agouti Signaling Protein) across three species, using data retrieved directly from the NCBI API.

Using Biopython, the script:

Downloads ASIP gene sequences in .fasta format for each species.

Performs pairwise local alignments (localxx) to evaluate genetic similarity.

Calculates the percent identity between sequences to assess genetic conservation.

This analysis enables exploration of how a gene associated with the expression of a shared phenotype such as brown coat coloration is evolutionarily conserved among mammals.
