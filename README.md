Capture-recapture-classify using genotypes for individual identity

Nonspatial application of the approach used in the genotype SPIM repo, with some differences. No spatial information is used (obviously!), observation model is different, and therefore the ID update is different. The observation model here is Bernoulli for detection then zero-truncated Poisson for number of samples given detection. I use a categorical sampler here instead of Metropolis-Hastings I used in genotype SPIM. I like this observation model better and should add it to genotype SPIM.

Currently M0 only with no variation in genotyping error probabilities as a function of covariates. I can add Mt, Mh, Mb, etc, and 
sample covariates on genotyping error probabilities.

Disclaimer: I haven't tested this extensively, but it appears the samplers are sufficient, at least when you generate consensus genotypes for each sample and only initialize samples to the same ID if their consensus genotypes match. Basically, you want to favor initializing nonmatching samples to different individuals and combine samples during convergence, rather than combining a bunch that do not match, and needing to split clusters during convergence. Second, if the observed genotype data is very sparse and the false positive probabilities are initialized near 1, particularly for heterozygotes, it looks like it can get stuck there. So, in this case, initialize the genotyping error rates to more realistic values or change the prior.

I may need to add a joint z-ID update to improve mixing, or remove bottlenecks preventing convergence, if they occur. Will investigate. Did not see any complete bottlenecks for genotype SPIM, but that uses spatial information not available here. None of the issues above were seen for genotype SPIM, spatial info is useful!
