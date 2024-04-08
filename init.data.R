init.data <- function(data=NA,M=NA,inits=inits,n.cluster.init=30,initTrue=FALSE){
  this.k <- data$this.k
  K <- data$K
  n.cov <- data$n.cov
  n.levels <- data$n.levels
  IDcovs <- data$IDlist$IDcovs
  
  n.samples <- length(this.k)
  G.obs <- data$G.obs
  if(!is.array(G.obs))stop("G.obs must be an array")
  n.rep <- dim(G.obs)[3]
  n.cov <- dim(G.obs)[2]
  
  ##pull out initial values for gamma
  gammaMat <- inits$gammaMat
  
  #initialize G.obs.true, the true sample-level full categorical identities
  #initializing to the most commonly observed sample by category values across observers
  G.obs.true <- matrix(NA,nrow=n.samples,ncol=n.cov)
  for(i in 1:n.samples){
    for(l in 1:n.cov){
      vec <- G.obs[i,l,]
      if(any(is.na(vec))){
        vec <- vec[-which(is.na(vec))]
      }
      if(length(vec)>0){
        choose <- as.integer(names(which(table(vec)==max(table(vec)))))
        if(length(choose)>1){
          pick <- sample(1:length(choose),1)
          choose <- choose[pick]
        }
        G.obs.true[i,l] <- choose
      }else{
        G.obs.true[i,l] <- sample(IDcovs[[l]],1,replace=FALSE,prob=gammaMat[l,1:n.levels[l]])
      }
    }
  }
  if(initTrue){
    ID <- data$ID
  }else{
    #find unique consensus genotypes
    G.obs.true2 <- apply(G.obs.true,1,paste,collapse="")
    unique.genos <- unique(G.obs.true2)
    n.unique.genos <- length(unique.genos)
    if(n.unique.genos>M)stop("More unique consensus genotypes them M, raise M to initialize or write your own algorithm")
    ID <- match(G.obs.true2,unique.genos)
    # G.clusters <- kmeans(G.obs.true,round(n.samples/2)) #use kmeans to cluster crude consensus genotypes
    # ID <- G.clusters$cluster
  }

  y.true <- matrix(0,M,K)
  for(l in 1:n.samples){
    y.true[ID[l],this.k[l]] <- y.true[ID[l],this.k[l]] + 1
  }
  
  #Initialize z
  z <- 1*(apply(y.true,1,sum)>0)
  N <- sum(z)
  
  #Initialize G.true
  G.true <- matrix(0, nrow=M,ncol=n.cov)
  for(i in 1:max(ID)){
    idx <- which(ID==i)
    if(length(idx)==1){
      G.true[i,] <- G.obs.true[idx,]
    }else if(length(idx)>1){
      if(ncol(G.obs.true)>1){
        G.true[i,] <- apply(G.obs.true[idx,],2, max) #consensus
      }else{
        G.true[i,] <- max(G.obs.true[idx,])
      }
    }#else need to fill below
  }
  #Fill in missing values, create indicator for them
  G.latent <- G.true==0
  for(j in 1:n.cov){
    fix <- G.true[,j]==0
    G.true[fix,j] <- sample(IDcovs[[j]],sum(fix),replace=TRUE,prob=gammaMat[j,1:n.levels[j]])
  }
  
  #Converting NAs in G.obs to some other arbitrary value that is not a genotype number to avoid NA warnings from nimble
  G.obs <- data$G.obs
  G.obs.NA.indicator <- is.na(G.obs)
  G.obs[is.na(G.obs)] <- 9999999
  if(n.rep==1){
    stop("n.rep==1 not currently handled, won't work without space.")
  }
  if(n.cov==1){
    stop("n.cov==1 not currently handled")
  }
  return(list(y.true=y.true,z=z,N=N,G.true=G.true,ID=ID,n.samples=n.samples,this.k=this.k,
              G.latent=G.latent,G.obs=G.obs,G.obs.NA.indicator=G.obs.NA.indicator))
}
