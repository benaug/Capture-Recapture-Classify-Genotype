Capture-recapture-classify using genotypes for individual identity

Nonspatial application of the approach used in the genotype SPIM repo, with some differences. No spatial information is used (obviously), observation model is different (Poisson more realistic with space, less so without), and therefore the ID update is different. The observation model here is Bernoulli for detection then zero-truncated Poisson for number of samples given detection. I use a categorical sampler here instead of Metropolis-Hastings I used in genotype SPIM. I like this observation model better and should add it to genotype SPIM.

Currently M0 only with no variation in genotyping error probabilities as a function of covariates. I can add Mt, Mh, Mb, etc, and 
sample covariates on genotyping error probabilities.

Disclaimer: I haven't tested this extensively, but it appears the samplers are sufficient, at least when you generate consensus genotypes (only doing so crudely, could be smarter) for each sample and only initialize samples to the same ID if their consensus genotypes match. A joint N/y.true update (analogous to z-ID update, but distributions are on N and y.true) might be needed to improve convergence/sampling. Will investigate.

Second, if the observed genotype data is very sparse and the false positive probabilities are initialized near 1, particularly for heterozygotes, it looks like it can get stuck there. The posterior is multimodal, but this mode is obviously the wrong one. In this case, initialize the genotyping error rates to realistic values or change the prior to reflect that false alleles are very rare.

Found and fixed a bug in the categorical sampler on 4/9. Appears to work well with sparse genotyping information now.
