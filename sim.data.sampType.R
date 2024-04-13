sim.data.sampType <-
  function(N=NA,p.y=NA,lambda.y=NA,K=NA,p.geno.het=NA,
           p.geno.hom=NA,n.loci=NA,n.rep=NA,pi.samp.type=NA,
           p.amp=NA,gamma=NA,IDcovs=NA,ptype=NA,seed=NA){
    #error checks
    if(length(gamma)!=n.loci)stop("gamma must be of length n.loci")
    for(l in 1:n.loci){
      if(length(gamma[[l]])!=length(IDcovs[[l]]))stop("gamma[[l]] must have one element per element of IDcovs[[l]]")
      if(sum(gamma[[l]])!=1)stop("gamma[[l]] must sum to 1")
    }
    samp.levels <- length(pi.samp.type)
    for(i in 1:samp.levels){
      if(sum(p.geno.hom[[i]])!=1)stop("p.geno.hom must sum to 1 for each sampType")
      if(sum(p.geno.het[[i]])!=1)stop("p.geno.het must sum to 1 for each sampType")
      if(length(p.geno.hom[[i]])!=2)stop("p.geno.hom must be of length 2 for each sampType")
      if(length(p.geno.het[[i]])!=3)warning("p.geno.het must be of length 3 unless simulating SNPs")
    }
    if(length(p.amp)!=samp.levels)stop("p.amp must be of length samp.levels")
    
    if(!is.na(seed)){
      set.seed(seed)
    }
    
    #simulate IDcovs
    n.levels <- unlist(lapply(gamma,length))
    G.true <- matrix(NA,nrow=N,ncol=n.loci) #all IDcovs in population.
    for(i in 1:N){
      for(j in 1:n.loci){
        G.true[i,j] <- sample(IDcovs[[j]],1,prob=gamma[[j]])
      }
    }
    library(VGAM)
    #Capture individuals
    y.det <- y.count <- matrix(0,N,K)
    for(i in 1:N){
      for(k in 1:K){
        y.det[i,k] <- rbinom(1,1,p.y)
        if(y.det[i,k]==1){
          y.count[i,k] <- rzapois(1,lambda.y,pobs0=0)
        }
      }
    }
    
    #discard uncaptured inds and aggregate true IDcovs for all samples, keeping track of where they came from with A matrix (used to put observed data back together)
    caught <- which(apply(y.count,c(1),sum)>0)
    y.true <- y.count
    y <- y.count[caught,]
    n <- length(caught)
    n.samples <- sum(y)
    G.true <- G.true[caught,]
    if(n.loci==1){
      G.true <- matrix(G.true,ncol=1)
    }
    G.cap <- matrix(NA,nrow=n.samples,ncol=n.loci)

    #disaggregate samples
    ID <- this.k <- rep(NA,n.samples)
    idx <- 1
    for(i in 1:n){ #loop through inds (uncaptured already removed)
      for(k in 1:K){ #then occasions
        if(y[i,k]>0){ #is there at least one sample here?
          for(l in 1:y[i,k]){ #then samples
            ID[idx] <- i #ID numbers don't count uncaptured guys
            this.k[idx] <- k
            G.cap[idx,]=G.true[i,]
            idx <- idx+1
          }
        }
      }
    }

    #double check that observed data can reconstruct true data
    y.obs <- matrix(0,max(ID),K)
    for(l in 1:n.samples){
      y.obs[ID[l],this.k[l]] <- y.obs[ID[l],this.k[l]] + 1
    }
    if(!all(y.obs==y)){
      stop("Error in data simulator!")
    }
    
    #simulate sample type covariates
    samp.type <- sample(samp.levels,nrow(G.cap),replace=TRUE,prob=pi.samp.type)
    #build genotype classification array
    theta <- vector("list",n.loci)
    for(m in 1:n.loci){
      theta[[m]] <- array(0,dim=c(samp.levels,n.levels[m],n.levels[m]))
      for(sl in 1:samp.levels){
        for(l in 1:n.levels[m]){
          if(!any(ptype[[m]][l,]==2)){#homozygote
            theta[[m]][sl,l,which(ptype[[m]][l,]==1)] <- (p.geno.hom[[sl]][1])
            theta[[m]][sl,l,which(ptype[[m]][l,]==3)] <- (p.geno.hom[[sl]][2])*(1/sum(ptype[[m]][l,]==3))
          }else{
            theta[[m]][sl,l,which(ptype[[m]][l,]==1)] <- (p.geno.het[[sl]][1])
            theta[[m]][sl,l,which(ptype[[m]][l,]==2)] <- (p.geno.het[[sl]][2])*(1/sum(ptype[[m]][l,]==2))
            theta[[m]][sl,l,which(ptype[[m]][l,]==3)] <- (p.geno.het[[sl]][3])*(1/sum(ptype[[m]][l,]==3))
          }
        }
      }
    }
    
    #loci-level amplification and genotyping observation process
    G.error <- array(NA,dim=c(n.samples,n.loci,n.rep))
    for(k in 1:n.rep){
      G.error[,,k] <- G.cap
      for(l in 1:n.loci){
        for(i in 1:n.samples){
          if(rbinom(1,1,p.amp[samp.type[i]])==1){ #scores missing completely at random
            for(j in 1:n.levels[l]){
              if(G.cap[i,l]==j){
                G.error[i,l,k] <- sample(IDcovs[[l]],1,prob=theta[[l]][samp.type[i],j,])
              }
            }
          }else{
            G.error[i,l,k] <- NA
          }
        }
      }
    }
    #find errors that occurred
    G.obstype <- array(0,dim=c(n.samples,n.loci,n.rep))
    for(k in 1:n.rep){
      for(l in 1:n.loci){
        for(i in 1:n.samples){
          if(is.na(G.error[i,l,k]))next#is missing
          if(G.error[i,l,k]==G.cap[i,l]){#is correct
            G.obstype[i,l,k] <- 1
          }else{#is an error
            G.obstype[i,l,k] <- ptype[[l]][G.cap[i,l],G.error[i,l,k]]
          }
        }
      }
    }
   
    out <- list(y=y,this.k=this.k,G.true=G.true,G.obs=G.error,n.loci=n.loci,n.levels=n.levels,
             n.samples=length(this.k),IDlist=list(n.loci=n.loci,IDcovs=IDcovs,ptype=ptype),
             ID=ID,K=K,n=nrow(y),samp.type=samp.type,G.obstype=G.obstype,seed=seed)
    
    return(out)
  }