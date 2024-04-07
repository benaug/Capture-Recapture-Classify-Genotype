Classical capture-recapture-classify using genotypes for individual identity

Analogous approach as used in genotype SPIM repo, with some differences. No spatial information used here, observation model is different, and the ID update is different. 
The observation model here is Bernoulli for detection then zero-truncated Poisson for number of samples given detection. 
I use a categorical sampler here instead of Metropolis-Hastings I used in genotype SPIM. 
I like this observation model better and should add it to genotype SPIM.

Currently M0 only with no variation in genotyping error probabilities as a function of covariates. 

I haven't tested this extensively, but it appears the samplers are sufficient. I may need to add a joint z-ID update to improve mixing. Will investigate.
