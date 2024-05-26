NimModel <- nimbleCode({
  ###priors##
  #process model
  lambda.N ~ dunif(0,1000) #expected abundance
  #genotype frequency priors
  for(m in 1:n.loci){
    for(k in 1:n.levels[m]){
      #dirichlet prior parameters, 1 is "all equal", could be smarter, can put distributions on log(alpha)
      alpha[m,k] <- 1
    }
    gammaMat[m,1:n.levels[m]] ~ ddirch(alpha[m,1:n.levels[m]]) #ragged matrix of genotype frequencies
  }
  
  #observation model - Bernoulli detection, zero-truncated poisson samples|detection
  logit(p.y) ~ dlogis(0,1) #detection probability
  lambda.y ~ dunif(0,5) #parameter for number of samples|detection
  
  #genotype observation process
  #genotyping error priors for heterozygote and homozygote loci-level genotypes for each sample type
  for(sl in 1:samp.levels){
    for(i in 1:3){
      #dirichlet prior parameters, 1 is "all equal", seems to work fine, could be smarter
      alpha.het[sl,i] <- 1
    }
    for(i in 1:2){
      #dirichlet prior parameters, 1 is "all equal", seems to work fine, could be smarter
      alpha.hom[sl,i] <- 1
    }
    p.geno.het[sl,1:3] ~ ddirch(alpha.het[sl,1:3])
    p.geno.hom[sl,1:2] ~ ddirch(alpha.hom[sl,1:2])
  }
  
  ##likelihoods##
  N ~ dpois(lambda.N) #realized abundance
  #data augmentation "under the hood", jointly update N/z, 
  #no distribution induced on z, just turns obsmod on/off, used in y.true/ID update
  for(i in 1:M){
    for(m in 1:n.loci){
      #Individual genotypes, all augmented individuals
      G.true[i,m] ~ dcat(gammaMat[m,1:n.levels[m]])
    }
    #Observation model is zero-truncated Poisson hurdle
    #vectorized over occasions, only evaluated when z[i]=1
    #y.true/ID are latent, updated "under the hood"
    y.true[i,1:K] ~ dHurdleVector(p.y=p.y,lambda.y=lambda.y,K=K,z=z[i])
  }
  #genotype classification probability array
  #getTheta() divides the probability of each error type by the number of ways
  #the error can be made and structures into classification array
  for(sl in 1:samp.levels){
    theta[1:n.loci,sl,1:max.levels,1:max.levels] <- getTheta(ptype = ptype[1:n.loci,1:max.levels,1:max.levels],
                                                          p.geno.het = p.geno.het[sl,1:3],
                                                          p.geno.hom = p.geno.hom[sl,1:2],
                                                          n.levels=n.levels[1:n.loci])
  }
  
  #genotype observation process, vectorized over reps
  for(l in 1:n.samples){
    for(m in 1:n.loci){
      #custom distribution to skip missing values.
      #ID updates won't always work if you sample unobserved data.
      G.obs[l,m,1:n.rep] ~ dcat2(theta=theta[m,samp.type[l],G.true[ID[l],m],1:n.levels[m]],
                                 n.rep = n.rep,na.ind=na.ind[l,m,1:n.rep]) 
    }
  }
  #calculate number of inds captured, intermediate object to compute n below
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:K])
  #must use ID and G.latent somewhere to make nimble happy. Sticking them here, not used in function.
  #G.latent used in custom G.true update
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples],G.latent=G.latent[1:M,1:n.loci])
})# end model
