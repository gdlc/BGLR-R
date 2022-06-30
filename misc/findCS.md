## A forward approach to identify credible sets in Bayesian Variable Selection Models

The function `findCS()`, available [here](https://github.com/gdlc/BGLR-R/blob/master/misc/findCS.R) can be used to detect sets of SNPs jointly associated with a phenotype.

The function take as inputs the samples from the posterior distribution of a variable selection model (e.g., BayesC or BayesB priors), map information (chromosome and physical position of the SNPs) plus a few algorithm parameters to identify credible sets.

Credible set inferences can be useful when one has groups of SNPs in LD with elevated inclusion probability with no-individual SNP reaching a sufficiently high inclusion probability.

The following [handout](https://www.dropbox.com/s/hwv6e83kl2d8aq9/HANDOUT_CredSetInference.pdf?dl=0) describes the types of problems where this function can be useful, and provides a toy-simulation example that illustrates how to use the function.

