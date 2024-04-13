#Zero-truncated Poisson Hurdle distribution for observation model
dHurdleVector <- nimbleFunction(
  run = function(x = double(1),p.y=double(0),lambda.y=double(0),z=double(0),K=double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z==1){
      for(k in 1:K){
        if(x[k]==0){ #if not captured, Bernoulli
          logProb <- logProb + dbinom(0,size=1,prob=p.y,log=TRUE)
        }else{ #if captured, Bernoulli and ZT Poisson
          logProb <- logProb + dbinom(1,size=1,prob=p.y,log=TRUE) + log(dpois(x[k],lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
        }
      }
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

#RNG dummy to make nimble happy
rHurdleVector <- nimbleFunction(
  run = function(n = integer(0),p.y=double(0),lambda.y=double(0),z=double(0),K=double(0)) {
    returnType(double(1))
    return(rep(0,K))
  }
)

#sum up number of captures of each individual
#intermediate objects that can be useful for
#joint z-ID updates, but not used. Used to compute
#number of guys captured, n
Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    K <- nimDim(y.true)[2]
    capcounts <- numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i] <- sum(y.true[i,1:K])
    }
    return(capcounts)
  }
)
#count up the number of captured individuals
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1),G.latent=double(2)){ #don't need ID or G.latent, but nimble requires is it used in a function 
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

#build theta. Not as efficient as possible because we rebuild both
#homozygotes and heterozygotes every time p.geno.het or p.geno.hom elements are updated.
#could block those updates.
getTheta <- nimbleFunction(
  run = function(ptype = double(3), p.geno.het = double(1), p.geno.hom = double(1), n.levels = double(1)) {
    returnType(double(3))
    n.loci <- nimDim(ptype)[1]
    max.levels <- nimDim(ptype)[2]
    thetaArray <- array(NA,dim=c(n.loci,max.levels,max.levels))
    for(m in 1:n.loci){
      for(l in 1:n.levels[m]){
        if(!any(ptype[m,l,1:n.levels[m]]==2)){#homozygote bc no possible allelic dropout
          tmp3 <- (1/sum(ptype[m,l,1:n.levels[m]]==3))
          for(l2 in 1:n.levels[m]){
            if(ptype[m,l,l2]==3){#false allele
              thetaArray[m,l,l2] <- p.geno.hom[2]*tmp3
            }else{#correct
              thetaArray[m,l,l2] <- p.geno.hom[1]
            }
          }
        }else{ #heterozygote
          tmp2 <- (1/sum(ptype[m,l,1:n.levels[m]]==2))
          tmp3 <- (1/sum(ptype[m,l,1:n.levels[m]]==3))
          for(l2 in 1:n.levels[m]){
            if(ptype[m,l,l2]==3){#false allele
              thetaArray[m,l,l2] <- p.geno.het[3]*tmp3
            }else if(ptype[m,l,l2]==2){#allelic dropout
              thetaArray[m,l,l2] <- p.geno.het[2]*tmp2
            }else{#correct
              thetaArray[m,l,l2] <- p.geno.het[1]
            }
          }
        }
      }
    }
    return(thetaArray)
  }
)

#Categorical distribution that skips missing values
dcat2 <- nimbleFunction(
  run = function(x = double(1), theta = double(1), n.rep = integer(0),
                 na.ind=double(1), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    for(rep in 1:n.rep){
      if(na.ind[rep]==FALSE){
        logProb <- logProb + log(theta[x[rep]])
      }
    }
    return(logProb)
  }
)
#dummy rng to make nimble happy, not used
rcat2 <- nimbleFunction(
  run = function(n = integer(0), theta = double(1), n.rep = integer(0),
                 na.ind=double(1)) {
    returnType(double(1))
    out <- numeric(n.rep,value=999)
    return(out)
  }
)

#Custom sampler to update G.true, using information in G.latent to determine proposal distribution
#Metropolis-Hastings update here allows other parameters to be a function of G.true
#Unsure if this is the most efficient way to calculate the likelihood for the proposal.
#I assume nimble is attempting to update the likelihood of every G.obs[l,m,rep] and skipping
#the samples not involved. But we know which ones are involved. "these.samps" below.
#Anyways, perhaps an inefficiency here. This is resolved in GSampler2 below.
GSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    n.levels <- control$n.levels
    n.rep <- control$n.rep
    na.ind <- control$na.ind
    calcNodes <- model$getDependencies(target)
    i <- control$i
    m <- control$m
    samp.type <- control$samp.type
  },
  run = function() {
    G.probs <- model$gammaMat[m,1:n.levels[m]] #pull out genotype frequencies
    if(model$G.latent[i,m]==FALSE){ #these individual-loci have samples allocated to them currently
      #build proposal distribution using genotype frequency likelihood and classification likelihood.
      these.samps <- which(model$ID==i)
      G.obs <- model$G.obs
      error.probs <- rep(1,n.levels[m]) #error probs|category level
      for(i2 in 1:length(these.samps)){
        for(obs in 1:n.rep){
          if(na.ind[these.samps[i2],obs]==FALSE){ #if observed
            error.probs <- error.probs*model$theta[m,samp.type[these.samps[i2]],1:n.levels[m],G.obs[these.samps[i2],m,obs]]
          }
        }
      }
      G.probs <- G.probs*error.probs
      G.probs <- G.probs/sum(G.probs)
    }else{ #these individual-loci do not have samples allocated to them currently
      #build proposal distribution using only genotype frequency likelihood. This is the
      #full conditional if no other parameters depend on G.true
      G.probs <- G.probs/sum(G.probs)
    }
    
    #MH step
    G.true.lp.initial <- model$getLogProb(calcNodes) #initial logProb
    prop.back <- G.probs[model$G.true[i,m]] #backwards proposal prob
    G.prop <- rcat(1,G.probs[1:n.levels[m]])
    prop.for <- G.probs[G.prop] #forwards proposal prob
    model$G.true[i,m] <<- G.prop #store in model
    G.true.lp.proposed <- model$calculate(calcNodes)#proposed logProb
    log_MH_ratio <- (G.true.lp.proposed+log(prop.back)) - (G.true.lp.initial+log(prop.for))
    # log_MH_ratio
    accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

#Alternative G.true sampler that uses much less RAM and is faster.
#However! No parameters other than G.obs can depend on G.true with this
#version unless you add their likelihoods to the MH update.
GSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    M <- control$M
    n.loci <- control$n.loci
    n.levels <- control$n.levels
    n.samples <- control$n.samples
    n.rep <- control$n.rep
    na.ind <- control$na.ind
    G.true.nodes <- control$G.true.nodes
    G.obs.nodes <- control$G.obs.nodes
    calcNodes <- control$calcNodes
    samp.type <- control$samp.type
  },
  run = function() {
    node.idx <- 1 #this is used to keep up with the correct nodes in G.true.nodes which MUST BE 1D, for each i and m
    #must loop over loci first, then individuals for this to work correctly.
    G.obs <- model$G.obs
    for(m in 1:n.loci){
      G.obs.m.idx <- seq((m-1)*n.samples+1,(m-1)*n.samples+n.samples,1) #get index for cov m for all samples
      for(i in 1:M){
        G.probs <- model$gammaMat[m,1:n.levels[m]] #pull out genotype frequencies
        if(model$G.latent[i,m]==FALSE){ #these individual-loci have samples allocated to them currently
          #must use MH
          #build proposal distribution using genotype frequency likelihood and classification likelihood.
          these.samps <- which(model$ID==i)
          G.obs.idx <- G.obs.m.idx[these.samps] #pull out nodes for these samples at cov m
          n.these.samps <- length(these.samps)
          error.probs <- rep(1,n.levels[m]) #error probs|category level
          for(i2 in 1:n.these.samps){
            for(obs in 1:n.rep){
              if(na.ind[these.samps[i2],m,obs]==FALSE){ #if observed
                error.probs <- error.probs*model$theta[m,samp.type[these.samps[i2]],1:n.levels[m],G.obs[these.samps[i2],m,obs]]
              }
            }
          }
          G.probs <- G.probs*error.probs
          G.probs <- G.probs/sum(G.probs)

          #MH step
          G.true.lp.initial <- model$getLogProb(G.true.nodes[node.idx]) #initial logProb for G.true
          G.obs.lp.initial <- model$getLogProb(G.obs.nodes[G.obs.idx]) #initial logProb for G.obs
          prop.back <- G.probs[model$G.true[i,m]] #backwards proposal prob
          G.prop <- rcat(1,G.probs[1:n.levels[m]]) #proposal
          if(G.prop!=model$G.true[i,m]){ #we can skip this if we propose the current value
            prop.for <- G.probs[G.prop] #forwards proposal prob
            model$G.true[i,m] <<- G.prop #store in model
            G.true.lp.proposed <- model$calculate(G.true.nodes[node.idx])#proposed logProb for G.true
            G.obs.lp.proposed <- model$calculate(G.obs.nodes[G.obs.idx])#proposed logProb for G.true
            log_MH_ratio <- (G.true.lp.proposed+G.obs.lp.proposed+log(prop.back)) -
              (G.true.lp.initial+G.obs.lp.initial+log(prop.for))
            accept <- decide(log_MH_ratio)
            if(accept) {
              mvSaved["G.true",1][i,m] <<- model[["G.true"]][i,m] #move to mvSaved
            } else {
              model[["G.true"]][i,m] <<- mvSaved["G.true",1][i,m] #set back to init
              model$calculate(G.true.nodes[node.idx]) #set log prob back to init
              model$calculate(G.obs.nodes[G.obs.idx]) #set log prob back to init
            }
          }
        }else{ #these individual-loci do not have samples allocated to them currently
          #build proposal distribution using only genotype frequency likelihood. This is the
          #full conditional if no other parameters depend on G.true. So always accept
          G.probs <- G.probs/sum(G.probs)
          G.prop <- rcat(1,G.probs[1:n.levels[m]])
          model$G.true[i,m] <<- G.prop #store in model
          model$calculate(G.true.nodes[node.idx]) #update G.true logprob. No G.obs logprob.
          mvSaved["G.true",1][i,m] <<- model[["G.true"]][i,m]
        }
        node.idx=node.idx+1 #increment
      }
    }
    #copy back to mySaved to update logProbs. should be done above, already though.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)


#y.true/ID update
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    K <- control$K
    G.obs <- control$G.obs
    n.loci <- control$n.loci
    n.rep <- control$n.rep
    n.levels <- control$n.levels
    n.samples <- control$n.samples
    samp.type <- control$samp.type
    this.k <- control$this.k
    na.ind <- control$na.ind
    calcNodes <- model$getDependencies(c("y.true","G.obs"))
  },
  run = function() {
    G.true <- model$G.true
    y.true <- model$y.true
    ID <- model$ID
    z <- model$z
    theta <- model$theta
    p.y <- model$p.y[1]
    lambda.y <- model$lambda.y[1]
    
    #Precalculate y likelihood in 2D (1D in model)
    lp.y <- matrix(0,nrow=M,ncol=K)
    for(i in 1:M){
      if(z[i]==1){
        for(k in 1:K){
          if(y.true[i,k]==0){ #if not captured
            lp.y[i,k] <- dbinom(0,size=1,prob=p.y,log=TRUE)
          }else{ #if captured
            lp.y[i,k] <- dbinom(1,size=1,prob=p.y,log=TRUE) + log(dpois(y.true[i,k],lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
          }
        }
      }
    }
    
    #categorical update
    for(l in 1:n.samples){
      # y.cand <- y.true #necessary for every sample for correct proposal probs
      #G.probs proportional to genotyping error likelihood over all loci and reps
      lp.G <- rep(0,M)
      lp.y.total <- rep(sum(lp.y[,this.k[l]]),M)*z
      for(i in 1:M){
        if(z[i]==1){
          for(m in 1:n.loci){
            for(rep in 1:n.rep){
              if(na.ind[l,m,rep]==FALSE){ #if observed
                lp.G[i] <- lp.G[i] + log(theta[m,samp.type[l],G.true[i,m],G.obs[l,m,rep]])
              }
            }
          }
          if(i!=ID[l]){#otherwise lp.y.total is already correct (current total)
            #move sample from old guy to new guy
            #old guy
            y.tmp1 <- y.true[ID[l],this.k[l]] - 1
            lp.tmp1 <- 0
            if(y.tmp1==0){ #if not captured
              lp.tmp1 <- dbinom(0,size=1,prob=p.y,log=TRUE)
            }else{ #if captured
              lp.tmp1 <- dbinom(1,size=1,prob=p.y,log=TRUE) + log(dpois(y.tmp1,lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
            }
            #new guy
            y.tmp2 <- y.true[i,this.k[l]] + 1
            lp.tmp2 <- 0
            if(y.tmp2==0){ #if not captured, will never happen in giving a guy this sample
              lp.tmp2 <- dbinom(0,size=1,prob=p.y,log=TRUE)
            }else{ #if captured
              lp.tmp2 <- dbinom(1,size=1,prob=p.y,log=TRUE) + log(dpois(y.tmp2,lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
            }
            #add likelihood difference for each guy from current.
            lp.y.total[i] <- lp.y.total[i] + (lp.y[ID[l],this.k[l]] - lp.tmp1) + (lp.y[i,this.k[l]] - lp.tmp2)
          }
        }else{ #can't propose this guy if z==0
          lp.y.total[i] <- lp.G[i] <- -Inf
        }
      }

      total.probs <- exp(lp.y.total + lp.G)
      total.probs <- total.probs/sum(total.probs)
      prop.i <- rcat(1,prob=total.probs)
      focal.i <- ID[l]

      if(ID[l]!=prop.i){ #skip if propose same ID
        ID[l] <- prop.i
        swapped <- c(focal.i,prop.i)
        #update y.true
        y.true[swapped[1],this.k[l]] <- y.true[swapped[1],this.k[l]]-1
        y.true[swapped[2],this.k[l]] <- y.true[swapped[2],this.k[l]]+1
        ##update lp.y
        for(i in 1:2){
          if(y.true[swapped[i],this.k[l]]==0){ #if not captured
            lp.y[swapped[i],this.k[l]] <- dbinom(0,size=1,prob=p.y,log=TRUE)
          }else{ #if captured
            lp.y[swapped[i],this.k[l]] <- dbinom(1,size=1,prob=p.y,log=TRUE) + log(dpois(y.true[swapped[i],this.k[l]],lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
          }
        }
      }
    }
    
    #update G.latent after ID changes
    G.latent <- matrix(TRUE,nrow=M,ncol=n.loci)
    for(l in 1:n.samples){
      for(m in 1:n.loci){
        for(rep in 1:n.rep){
          if(na.ind[l,m,rep]==FALSE){ #if this sample is not latent
            G.latent[ID[l],m] <- FALSE #then its ID is not latent
          }
        }
      }
    }
    #put everything back into the model$stuff after updating y.sight.true, y.sight.true.event
    model$y.true <<- y.true
    model$ID <<- ID
    model$G.latent <<- G.latent
    model$calculate(calcNodes) #update dependencies, likelihoods
    #Need to update ID in mvSaved manually, can't calculate. 
    #only required if using again in custom update, which we aren't
    mvSaved["ID",1] <<- model[["ID"]]
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#Required custom update for N/z (unless you modify this back to regular data augmentation)
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        # find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(model$capcounts[pick]>0){#is this an individual with samples?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)