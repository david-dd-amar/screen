# screen
Empirical Bayes algorithms for controlling the false discovery rate

We present new algorithms and tools for Empirical Bayes replicability analysis of large scale data measured in many experiments.

The tools take a matrix of p-values in whichs rows are objects (e.g., genes or snps) and columns are studies. We then report the fdr of rows for a desired level of replication (e.g., at least 20% of the studies), taking into consideration dependencies among the studies. These are inferred without the need to specify prior knowledge.

For more details, see our [paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005700#sec029).
