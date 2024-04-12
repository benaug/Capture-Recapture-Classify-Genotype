sim.data.Mh <-
  function(N=NA,beta0.p.y=NA,sd.p.y=NA,lambda.y=NA,K=NA,p.geno.het=NA,
           p.geno.hom=NA,n.loci=NA,n.rep=NA,
           pID=NA,gamma=NA,IDcovs=NA,ptype=NA,seed=NA){
    #error checks
    if(length(gamma)!=n.loci)stop("gamma must be of length n.loci")
    for(l in 1:n.loci){
      if(length(gamma[[l]])!=length(IDcovs[[l]]))stop("gamma[[l]] must have one element per element of IDcovs[[l]]")
      if(sum(gamma[[l]])!=1)stop("gamma[[l]] must sum to 1")
    }
    if(sum(p.geno.hom)!=1)stop("p.geno.hom must sum to 1")
    if(sum(p.geno.het)!=1)stop("p.geno.het must sum to 1")
    if(length(p.geno.hom)!=2)stop("p.geno.hom must be of length 2")
    if(length(p.geno.het)!=3)warning("p.geno.het must be of length 3 unless simulating SNPs")
    
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
    p.y <- plogis(rnorm(N,beta0.p.y,sd.p.y))
    y.det <- y.count <- matrix(0,N,K)
    for(i in 1:N){
      for(k in 1:K){
        y.det[i,k] <- rbinom(1,1,p.y[i])
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
    # ID <- rep(1:nrow(y),times=rowSums(y))
    
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
    theta=vector("list",n.loci)
    for(m in 1:n.loci){
      theta[[m]] <- matrix(0,nrow=n.levels[m],ncol=n.levels[m])
      for(l in 1:n.levels[m]){
        if(!any(ptype[[m]][l,]==2)){#homozygote
          theta[[m]][l,which(ptype[[m]][l,]==1)] <- (p.geno.hom[1])
          theta[[m]][l,which(ptype[[m]][l,]==3)] <- (p.geno.hom[2])*(1/sum(ptype[[m]][l,]==3))
        }else{
          theta[[m]][l,which(ptype[[m]][l,]==1)] <- (p.geno.het[1])
          theta[[m]][l,which(ptype[[m]][l,]==2)] <- (p.geno.het[2])*(1/sum(ptype[[m]][l,]==2))
          theta[[m]][l,which(ptype[[m]][l,]==3)] <- (p.geno.het[3])*(1/sum(ptype[[m]][l,]==3))
        }
      }
    }
    
    #observation error
    G.error <- array(NA,dim=c(n.samples,n.loci,n.rep))
    for(k in 1:n.rep){
      G.error[,,k] <- G.cap
      for(l in 1:n.loci){
        for(i in 1:n.samples){
          if(rbinom(1,1,pID[l])==1){
            for(j in 1:n.levels[l]){
              if(G.cap[i,l]==j){
                G.error[i,l,k] <- sample(IDcovs[[l]],1,prob=theta[[l]][j,])
              }
            }
          }else{
            G.error[i,l,k] <- NA
          }
        }
      }
    }
    #find errors that occurred
    G.Obstype <- array(0,dim=c(n.samples,n.loci,n.rep))
    for(k in 1:n.rep){
      for(l in 1:n.loci){
        for(i in 1:n.samples){
          if(is.na(G.error[i,l,k]))next#is missing
          if(G.error[i,l,k]==G.cap[i,l]){#is correct
            G.Obstype[i,l,k] <- 1
          }else{#is an error
            G.Obstype[i,l,k] <- ptype[[l]][G.error[i,l,k],G.cap[i,l]]
          }
        }
      }
    }
    getmode  <-  function(v) {
      uniqv  <-  unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    #generate a crude consensus genotype
    G.consensus <- apply(G.error,c(1,2),function(x){getmode(x[x!=0])})
    #how many of the crude consensus genotypes are corrupted?
    corrupted <- sum(apply(G.consensus==G.cap,1,function(x){any(x==FALSE)}),na.rm=TRUE)
    
    out <- list(y=y,this.k=this.k,G.true=G.true,G.obs=G.error,n.loci=n.loci,n.levels=n.levels,
             n.samples=length(this.k),IDlist=list(n.loci=n.loci,IDcovs=IDcovs,ptype=ptype),
             ID=ID,K=K,n=nrow(y),corrupted=corrupted,G.Obstype=G.Obstype,p.y=p.y,seed=seed)
    
    return(out)
  }