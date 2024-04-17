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
  run = function(ptype = double(2), p.geno.het = double(1), p.geno.hom = double(1)) {
    returnType(double(2))
    thetaMatrix <- matrix(NA,3,3)
    for(l in 1:3){
      if(!any(ptype[l,]==2)){#homozygote bc no possible allelic dropout
        tmp3 <- 1/sum(ptype[l,]==3)
        for(l2 in 1:3){
          if(ptype[l,l2]==3){#false allele
            thetaMatrix[l,l2] <- p.geno.hom[2]*tmp3
          }else{#correct
            thetaMatrix[l,l2] <- p.geno.hom[1]
          }
        }
      }else{ #heterozygote
        tmp2 <- 1/sum(ptype[l,]==2)
        tmp3 <- 1/sum(ptype[l,]==3)
        for(l2 in 1:3){
          if(ptype[l,l2]==3){#false allele
            print("error, SNPs shouldn't have false allele events for true heterozygotes")
          }else if(ptype[l,l2]==2){#allelic dropout
            thetaMatrix[l,l2] <- p.geno.het[2]*tmp2
          }else{#correct
            thetaMatrix[l,l2] <- p.geno.het[1]
          }
        }
      }
    }
    return(thetaMatrix)
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
    n.rep <- control$n.rep
    na.ind <- control$na.ind
    calcNodes <- model$getDependencies(target)
    i <- control$i
    m <- control$m
  },
  run = function() {
    G.probs <- model$gammaMat[m,1:3] #pull out genotype frequencies
    if(model$G.latent[i,m]==FALSE){ #these individual-loci have samples allocated to them currently
      #build proposal distribution using genotype frequency likelihood and classification likelihood.
      these.samps <- which(model$ID==i)
      G.obs <- model$G.obs
      error.probs <- rep(1,3) #error probs|category level
      for(i2 in 1:length(these.samps)){
        for(obs in 1:n.rep){
          if(na.ind[these.samps[i2],obs]==FALSE){ #if observed
            error.probs <- error.probs*model$theta[1:3,G.obs[these.samps[i2],m,obs]]
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
    G.prop <- rcat(1,G.probs[1:3])
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
    n.samples <- control$n.samples
    n.rep <- control$n.rep
    na.ind <- control$na.ind
    G.true.nodes <- control$G.true.nodes
    G.obs.nodes <- control$G.obs.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    node.idx <- 1 #this is used to keep up with the correct nodes in G.true.nodes which MUST BE 1D, for each i and m
    #must loop over loci first, then individuals for this to work correctly.
    G.obs <- model$G.obs
    for(m in 1:n.loci){
      G.obs.m.idx <- seq((m-1)*n.samples+1,(m-1)*n.samples+n.samples,1) #get index for cov m for all samples
      for(i in 1:M){
        G.probs <- model$gammaMat[m,1:3] #pull out genotype frequencies
        if(model$G.latent[i,m]==FALSE){ #these individual-loci have samples allocated to them currently
          #must use MH
          #build proposal distribution using genotype frequency likelihood and classification likelihood.
          these.samps <- which(model$ID==i)
          G.obs.idx <- G.obs.m.idx[these.samps] #pull out nodes for these samples at cov m
          n.these.samps <- length(these.samps)
          error.probs <- rep(1,3) #error probs|category level
          for(i2 in 1:n.these.samps){
            for(obs in 1:n.rep){
              if(na.ind[these.samps[i2],m,obs]==FALSE){ #if observed
                error.probs <- error.probs*model$theta[1:3,G.obs[these.samps[i2],m,obs]]
              }
            }
          }
          G.probs <- G.probs*error.probs
          G.probs <- G.probs/sum(G.probs)

          #MH step
          G.true.lp.initial <- model$getLogProb(G.true.nodes[node.idx]) #initial logProb for G.true
          G.obs.lp.initial <- model$getLogProb(G.obs.nodes[G.obs.idx]) #initial logProb for G.obs
          prop.back <- G.probs[model$G.true[i,m]] #backwards proposal prob
          G.prop <- rcat(1,G.probs[1:3]) #proposal
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
          G.prop <- rcat(1,G.probs[1:3])
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
    n.samples <- control$n.samples
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
    #could pull this out of model instead, vectorized over reps...
    lp.G <- array(0,dim=c(n.samples,n.loci,n.rep))
    for(l in 1:n.samples){
      for(m in 1:n.loci){
        for(rep in 1:n.rep){
          if(na.ind[l,m,rep]==FALSE){
            lp.G[l,m,rep] <- dcat(G.obs[l,m,rep],theta[G.true[ID[l],m],1:3],log=TRUE)
          }
        }
      }
    }
    lp.y.cand <- lp.y
    lp.G.cand <- lp.G
    ID.cand <- ID
    y.cand <- y.true
    
    for(l in 1:n.samples){
      #get proposal distribution for sample l
      lp.G.prop <- rep(0,M) #considers how well sample l matches all G.true in population (symmetric)
      lp.y.prop <- rep(0,M) #considers how well sample l fits observation model (not symmetric)
      for(i in 1:M){
        if(z[i]==1){
          for(m in 1:n.loci){
            for(rep in 1:n.rep){
              if(na.ind[l,m,rep]==FALSE){ #if observed
                lp.G.prop[i] <- lp.G.prop[i] + log(theta[G.true[i,m],G.obs[l,m,rep]])
              }
            }
          }
          if(i!=ID[l]){#new state
            y.tmp1 <- y.true[i,this.k[l]] + 1 #if we add sample here
            y.tmp2 <- y.true[ID[l],this.k[l]] - 1 #if we add sample here
          }else{
            y.tmp1 <- y.true[i,this.k[l]] #current state
          }
          if(y.tmp1==0){ #if not captured
            lp.y.prop[i] <- dbinom(0,size=1,prob=p.y,log=TRUE)
          }else{ #if captured
            lp.y.prop[i] <- dbinom(1,size=1,prob=p.y,log=TRUE) + 
              log(dpois(y.tmp1,lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
          }
          if(i!=ID[l]){
            if(y.tmp2==0){ #if not captured
              lp.y.prop[i] <- lp.y.prop[i] + dbinom(0,size=1,prob=p.y,log=TRUE)
            }else{ #if captured
              lp.y.prop[i] <- lp.y.prop[i] + dbinom(1,size=1,prob=p.y,log=TRUE) + 
                log(dpois(y.tmp2,lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
            }
          }
        }else{ #can't propose this guy if z==0
          lp.y.prop[i] <- lp.G.prop[i] <- -Inf
        }
      }
      prop.probs <- exp(lp.G.prop + lp.y.prop)
      prop.probs <- prop.probs/sum(prop.probs)
      prop.i <- rcat(1,prob=prop.probs)
      focal.i <- ID[l]
      
      if(ID[l]!=prop.i){ #skip if propose same ID
        #update this ID
        ID.cand[l] <- prop.i
        swapped <- c(focal.i,prop.i)
        #update y.true
        y.cand[swapped[1],this.k[l]] <- y.true[swapped[1],this.k[l]] - 1
        y.cand[swapped[2],this.k[l]] <- y.true[swapped[2],this.k[l]] + 1
        ##update lp.y
        for(i in 1:2){
          if(y.cand[swapped[i],this.k[l]]==0){ #if not captured
            lp.y.cand[swapped[i],this.k[l]] <- dbinom(0,size=1,prob=p.y,log=TRUE)
          }else{ #if captured
            lp.y.cand[swapped[i],this.k[l]] <- dbinom(1,size=1,prob=p.y,log=TRUE) + log(dpois(y.cand[swapped[i],this.k[l]],lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
          }
        }
        #update lp.G
        for(m in 1:n.loci){
          for(rep in 1:n.rep){
            if(na.ind[l,m,rep]==FALSE){
              lp.G.cand[l,m,rep] <- dcat(G.obs[l,m,rep],theta[G.true[ID.cand[l],m],1:3],log=TRUE)
            }else{
              lp.G.cand[l,m,rep] <- 0
            }
          }
        }
        
        #get backwards proposal distribution for sample l
        #genotype model probs are symmetric, but observation model probs are not
        lp.y.prop.back <- rep(0,M) #considers how well sample l fits observation model (not symmetric)
        for(i in 1:M){
          if(z[i]==1){
            if(i!=ID.cand[l]){#new state
              y.tmp1 <- y.cand[i,this.k[l]] + 1 #if we add sample here
              y.tmp2 <- y.cand[ID.cand[l],this.k[l]] - 1 #if we add sample here
            }else{
              y.tmp1 <- y.cand[i,this.k[l]] #current state
            }
            if(y.tmp1==0){ #if not captured
              lp.y.prop.back[i] <- dbinom(0,size=1,prob=p.y,log=TRUE)
            }else{ #if captured
              lp.y.prop.back[i] <- dbinom(1,size=1,prob=p.y,log=TRUE) + 
                log(dpois(y.tmp1,lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
            }
            if(i!=ID.cand[l]){
              if(y.tmp2==0){ #if not captured
                lp.y.prop.back[i] <- lp.y.prop.back[i] + dbinom(0,size=1,prob=p.y,log=TRUE)
              }else{ #if captured
                lp.y.prop.back[i] <- lp.y.prop.back[i] + dbinom(1,size=1,prob=p.y,log=TRUE) + 
                  log(dpois(y.tmp2,lambda=lambda.y)/(1-dpois(0,lambda=lambda.y)))
              }
            }
          }else{ #can't propose this guy if z==0
            lp.y.prop.back[i] <- -Inf
          }
        }
        
        prop.probs.back <- exp(lp.G.prop + lp.y.prop.back)
        prop.probs.back <- prop.probs.back/sum(prop.probs.back)
        
        prop.prob.for <- prop.probs[swapped[2]]
        prop.prob.back <- prop.probs.back[swapped[1]]
        
        #probability we select this y[i,k] to update by selecting an ID at random
        select.prob.for <- sum(ID==ID[l]&this.k==this.k[l])/n.samples
        select.prob.back <- sum(ID.cand==ID.cand[l]&this.k==this.k[l])/n.samples
        
        if(runif(1)<exp((sum(lp.y.cand[swapped,this.k[l]]) + sum(lp.G.cand[l,,]))-
                        (sum(lp.y[swapped,this.k[l]]) + sum(lp.G[l,,])))*
           (prop.prob.back/prop.prob.for)*(select.prob.back/select.prob.for)){
          y.true[swapped,this.k[l]] <- y.cand[swapped,this.k[l]]
          lp.y[swapped,this.k[l]] <- lp.y.cand[swapped,this.k[l]]
          lp.G[l,,] <- lp.G.cand[l,,]
          ID[l] <- ID.cand[l]
        }else{#set these back
          y.cand[swapped,this.k[l]] <- y.true[swapped,this.k[l]]
          lp.y.cand[swapped,this.k[l]] <- lp.y[swapped,this.k[l]]
          lp.G.cand[l,,] <- lp.G[l,,]
          ID.cand[l] <- ID[l]
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