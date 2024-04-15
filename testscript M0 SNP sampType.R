#This version is for SNPs. I've never run it by a geneticist, so tell me if it needs modification
#I'm assuming the only heterozygote, 1-2, can be observed as 1-1 and 2-2 (allelic dropout), and
#then each homozygote 1-1, and 2-2, could be observed as a false allele, 1-2, 1-1 (for true 2-2), or 2-2 (for true 1-1)
#I'm assuming the false allele rate is the same for 1-1 and 2-2, and all error rates are the same over all SNP loci

#This version allows genotyping error rates to vary by categorical sample covariates
#continuous covariates can be done in practice, but very slow-need a theta matrix for
#every sample, rather than one for each category level

#I have the data simulator set to only 15 SNPs. If SNPs not informative enough
#(number of SNPs, number of reps, amplification rate, genotyping error rates), you may get positive bias

#Don't try to add time or individual detection effects in M0 scripts
#unless you make corresponding changes to custom updates

library(coda)
library(nimble)
source("sim.data.sampType.R")
source("init.data.R")
source("build.genos.R")
source("map.genos.R")
source("NimbleModel M0 SNP sampType.R")
source("Nimble Functions M0 SNP sampType.R")

#must run this line for nimble to do what we want
nimbleOptions(determinePredictiveNodesInModel = FALSE)

#First, let's structure some SNPs.
n.loci <- 50 #number of SNP loci, 25 will leave uncertainty in ID with settings below
unique.genos <- vector("list")
for(m in 1:n.loci){
  unique.genos[[m]] <- matrix(c(1,1,2,2,1,2),nrow=3,byrow=TRUE)
}
unique.genos[[1]] #11, 22, 12


#This function creates objects to determine which classifications are 1) correct, 2) false allele, and 3) allelic dropout
built.genos <- build.genos(unique.genos)
ptype <- built.genos$ptype #list of length n.loci with each element being an 3 x 3 matrix containing 
#indicator for error type for true genotypes along rows and classified genotype along columns

#Normal capture-recapture stuff
N <- 75 #realized abundance
p.y <- 0.2 #capture probability
lambda.y <- 1 #expected number of samples given capture (ZT Poisson)
K <- 5 #number of capture occasions
n.rep <- 3 #number of PCR reps per sample. This repo assumes at least 2 (1 allowed in genoSPIM, but generally need replication)

IDcovs <- vector("list",n.loci) #enumerating genotypes here for simulation and data initialization
for(i in 1:n.loci){
  IDcovs[[i]] <- 1:nrow(unique.genos[[i]])
}
gamma <- vector("list",n.loci)
for(i in 1:n.loci){
  #This simulates equal genotype frequencies. I don't have SNP data to use, probably not realistic
  gamma[[i]] <- rep(1/3,3)
}

#Genotype observation process parameters. Can have failed amplification (missing completely at random) and genotyping error
#Difference from microsat code here: true heterozygotes don't have false allele events (there is only 1 heterozygote!)
p.amp <- c(0.999,0.25) #sample by replication amplification probabilities (controls level of missing scores in G.obs)
samp.levels <- 2 #number of sample type covariates. Each type has it's own genotyping error rates.
p.geno.het <- vector("list",samp.levels)
p.geno.hom <- vector("list",samp.levels)
#P(correct, allelic dropout) for heterozygotes
p.geno.het[[1]] <- c(0.99,0.01) #high quality
p.geno.het[[2]] <- c(0.65,0.35) #low quality
#P(correct,false allele) for homozygotes
p.geno.hom[[1]] <- c(0.99,0.01) #high quality
p.geno.hom[[2]] <- c(0.65,0.35) #low quality

pi.samp.type <- c(0.52,0.48) #frequencies of each sample type



#can use same data simulator for MSATs for SNPs, will give you a warning that p.geno.het is only of length 2.
data <- sim.data.sampType(N=N,p.y=p.y,lambda.y=lambda.y,K=K,#cap-recap parameters/constants
                  n.loci=n.loci,p.amp=p.amp,n.rep=n.rep,
                  p.geno.hom=p.geno.hom,p.geno.het=p.geno.het,
                  gamma=gamma,IDcovs=IDcovs,ptype=ptype,pi.samp.type=pi.samp.type)

#The observed data are 
#1) the occasion of every "count member". E.g., a count of 3 
#(3 genetic samples at 1 occasion) has 3 count members
head(data$this.k)
#2) observed genotype replicates for every sample. Can have missing data indicated with NA
t(data$G.obs[1,,]) #observed genotypes for 1st count member

#Can use the map function to see which genotypes the enumerated genotypes correspond to
ind <- 1 #change ind number ot look at different individuals
data$G.true[ind,] #True genotype of individual 1, enumerated
map.genos(data$G.true[ind,],unique.genos) #converted back to actual genotypes
#can compare observed genotypes to true genotypes above
#allelic dropout events are when heterozygotes are observed as homozygotes. false allele is any other error.
these.samps <- which(data$ID==ind)
for(i in 1:length(these.samps)){
  print(map.genos(t(data$G.obs[these.samps[i],,]),unique.genos))
}

##OK, let's fit the model

#Data augmentation level - must be larger than N. 
#If N ever hits M during MCMC after convergence, raise M and start over
M <- 200
if(M<N)stop("M must be larger than simulate N")

#set some gamma inits. Using equal across locus-level genotypes here
#note, gamma is a ragged matrix for use in nimble.
gammaMat <- matrix(0,nrow=n.loci,ncol=3)
for(l in 1:n.loci){
  gammaMat[l,1:3] <- rep(1/3,3)
}

#provide some ballpark inits for gamma to help initialize data
inits <- list(gammaMat=gammaMat) #using gammaMat inits to initialize G.true not determined by observed data
nimbuild <- init.data(data=data,M=M,inits=inits,initTrue=FALSE) #can initialize simulated data sets at true IDs for testing

n.samples <- data$n.samples

#ptype is a ragged array that tells use which genotype classifications are correct, allelic dropout, or false alleles
#for SNPs, these are the same for all loci, so we'll convert to a matrix
ptypeMatrix <- built.genos$ptypeArray[1,,]

#OK, what is the data we use to fit the model?
#1) data$G.obs, the observed genotypes
#2) y.true/ID. Note, you can build y.true with ID, the individual of each sample, and this.k, the occasion of each sample.
#we specify a distribution on y.true, not ID. If we initialize ID, y.true is initialized.
#If we update ID, y.true is updated. So, we just initialize these, provide to nimble
#as inits, and update them on each MCMC iteration

#supply data to nimble
Nimdata <- list(G.obs=nimbuild$G.obs,samp.type=data$samp.type) #adding sample type covariates here

#inits for nimble
#Without spatial information, the possibility of false alleles makes convergence difficult with sparse G.obs
#Helps to provide ballpark inits for p.geno.het especially, providing them for p.geno.hom here, too.
#Comments above are likely less relevant to SNPs where heterozygotes can't have a false positive,
#but setting some ballpark inits here.
p.geno.het.init <- matrix(NA,nrow=2,ncol=2)
p.geno.hom.init <- matrix(NA,nrow=2,ncol=2)
p.geno.het.init[1,] <- c(0.9,0.1)
p.geno.het.init[2,] <- c(0.9,0.1)
p.geno.hom.init[1,] <- c(0.9,0.1)
p.geno.hom.init[2,] <- c(0.9,0.1)

Niminits <- list(z=nimbuild$z,N=nimbuild$N, #must initialize N to be sum(z) for this data augmentation approach
                 G.true=nimbuild$G.true,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                 y.true=nimbuild$y.true,G.latent=nimbuild$G.latent,
                 p.geno.het=p.geno.het.init,p.geno.hom=p.geno.hom.init,
                 gammaMat=gammaMat)

#constants for Nimble
constants <- list(M=M,K=K,n.samples=n.samples,n.loci=n.loci,n.rep=n.rep,
                na.ind=nimbuild$G.obs.NA.indicator, #tells nimble which observed genotype scores are missing
               ptype=ptypeMatrix,samp.levels=samp.levels)

# set parameters to monitor
parameters <- c('lambda.N','p.y','lambda.y','N','n','p.geno.het','p.geno.hom','gammaMat')

#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c('ID',"G.true") #monitoring the IDs and the true genotypes so we can explore them below
nt <- 1 #thinning rate
nt2 <- 50 #thin more

# Build the model, configure the mcmc, and compile
# can ignore warnings about 1) ID in constants 2) possible size mismatch for G.obs.
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#tell nimble which nodes to configure so we don't waste time for samplers we will replace below
#if you add parameters to the model file, need to add them here.
config.nodes <- c('log_lambda.N','logit_p.y','lambda.y','p.geno.het','p.geno.hom','gammaMat')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = FALSE,
                      monitors2=parameters2,thin2=nt2,
                      nodes=config.nodes) 

#Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
#Remove these if assigned above
#conf$removeSampler("G.obs")
#conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",K,"]"),
                type = 'IDSampler',control = list(M=M,K=K,n.loci=n.loci,n.samples=n.samples,samp.type=data$samp.type,
                                                  n.rep=n.rep,this.k=nimbuild$this.k,G.obs=data$G.obs,
                                                  na.ind=nimbuild$G.obs.NA.indicator),
                silent = TRUE)

#replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
# conf$removeSampler("G.true")
# # this is the "safe" version. It will use a lot of RAM, but will be correct if any parameters
# # depend on G.true besides G.obs, which seems unlikely for genotypes. But if, say, G.true[,1] is "sex", and
# # you specify that sigma varies by sex, this update is correct and the more efficient one below will not be.
# for(i in 1:M){
#   for(m in 1:n.loci){
#     conf$addSampler(target = paste("G.true[",i,",",m,"]", sep=""),
#                     type = 'GSampler',
#                     control = list(i = i,m=m,n.rep=n.rep,samp.type=data$samp.type,
#                                    na.ind=nimbuild$G.obs.NA.indicator[,m,]), silent = TRUE)
#   }
# }
# this is the low RAM version. No parameters can depend on G.true except G.obs
# identify G.true nodes here. Must be in matrix with individuals down rows and loci across columns.
# This update only works with "reps" vectorized in bugs code. Must modify this sampler if you unvectorize those.
G.true.nodes <- Rmodel$expandNodeNames(paste0("G.true[1:",M,",1:",n.loci,"]"))
G.obs.nodes <- Rmodel$expandNodeNames(paste0("G.obs[1:",n.samples,",1:",n.loci,",1:",n.rep,"]"))
calcNodes <- c(G.true.nodes,G.obs.nodes)
conf$addSampler(target = paste0("G.true[1:",M,",1:",n.loci,"]"),
                type = 'GSampler2',
                control = list(M=M,n.loci=n.loci,n.rep=n.rep,samp.type=data$samp.type,
                               na.ind=nimbuild$G.obs.NA.indicator,n.samples=nimbuild$n.samples,
                               G.true.nodes=G.true.nodes,G.obs.nodes=G.obs.nodes,
                               calcNodes=calcNodes), silent = TRUE)

###required sampler replacement for "alternative data augmentation" z/ID update (using distributions on N/y.true)
# how many N/z proposals per iteration? Not sure what is optimal, setting to 50% of M here.
# optimal should be higher when uncertainty in N is higher.
z.ups <- round(M*0.5) 
# conf$removeSampler("N")
#nodes used for update, calcNodes + z nodes
y.nodes <- Rmodel$expandNodeNames(paste("y.true[1:",M,",1:",K,"]"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,y.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,
                                                 y.nodes=y.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),
                silent = TRUE)

#can block these if highly correlated.
#often works best if you use both the independent and block updates, so maybe don't remove the independent updates.
# conf$removeSampler(c("log_lambda.N","logit_p.y"))
conf$addSampler(target = c("log_lambda.N","logit_p.y"),
type = 'RW_block',control=list(adaptive=TRUE,tries=1),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
#Can ignore nimble warnings about NA or NaN in ptype and theta
#Can ignore nimble warnings about G.obs value NA or NaN, due to padding to keep dimensions constant for nimble
start.time2 <- Sys.time()
#starting with short run
Cmcmc$run(2000,reset=FALSE) #can extend run by rerunning this line, e.g. run longer if not converged, or want more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

#remove gammaMat posteriors (not that interesting and tons of them) and plot
idx <- grep("gammaMat",colnames(mvSamples))
plot(mcmc(mvSamples[500:nrow(mvSamples),-idx]))

data$n #number of individuals captured to compare to posterior for n. No uncertainty with enough genotype info.

##Explore ID posteriors - not really interesting if no uncertainty
#Assuming ID posterior was monitored in mvSamples2
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
nrow(mvSamples2) #Need enough posterior iterations for reliable inference. If not, reduce thinning and/or run longer
idx <- grep("ID",colnames(mvSamples2))
plot(mcmc(mvSamples2[2:nrow(mvSamples2),idx]))

library(MCMCglmm)
burnin <- 50
IDpost <- round(posterior.mode(mcmc(mvSamples2[burnin:nrow(mvSamples2),idx])))
#For simulated data sets, comparing posterior mode ID to truth.
#Numbers will not be the same, but all samples with same true ID will have
#same ID in posterior mode when posterior mode is exactly correct. Numbers just don't match up.
cbind(data$ID,round(IDpost))

#calculate posterior probability of pairwise sample matches
#P(sample x belongs to same individual as sample y)
burnin <- 50 #where to start. Don't start at 1, is NA.
n.iter <- nrow(mvSamples2)-burnin+1
pair.probs <- matrix(NA,n.samples,n.samples)
for(i in 1:n.samples){
  for(j in 1:n.samples){
    count <- 0
    for(iter in burnin:n.iter){
      count <- count+1*(mvSamples2[iter,idx[j]]==mvSamples2[iter,idx[i]])
    }
    pair.probs[i,j] <- count/(n.iter-burnin+1)
  }
}

this.samp <- 1 #sample number to look at
# this.samp <- this.samp + 1
pair.probs[this.samp,] #probability this sample is from same individual as all other samples
pair.probs[this.samp,data$ID==data$ID[this.samp]] #for simulated data, these are the other samples truly from same individual

#inspect G.true (true genotype) posteriors
idx <- grep("G.true",colnames(mvSamples2))

#posterior mode of true genotypes. Note, this is the posterior mode of each loci individually
#I expect this should usually be the same at the posterior mode complete genotype, but this
#can also be calculated from this posterior
library(MCMCglmm)
burnin <- 50
G.mode <- round(posterior.mode(mcmc(mvSamples2[burnin:nrow(mvSamples2),idx])))
G.mode <- matrix(G.mode,nrow=M)
#rearrange all G.true samples to look at range of values instead of just the mode
G.samps <- mvSamples2[burnin:nrow(mvSamples2),idx]
G.samps <- array(t(G.samps),dim=c(M,n.loci,nrow(G.samps)))

#look at posterior mode genotype of each individual (not numbered the same as in true data,
#but numbers are the same as in IDpost)
ind <- 1 #change ind number to look at different individuals
G.mode[ind,] #True genotype of focal individual, enumerated
map.genos(G.mode[ind,],unique.genos) #converted back to actual genotypes

#which samples were most commonly assigned to this individual? (assumes you calculated IDpost above)
these.samps <- which(IDpost==ind)
if(length(these.samps>0)){
  for(i in 1:length(these.samps)){
    print(map.genos(t(data$G.obs[these.samps[i],,]),unique.genos))
  }
}else{
  "No sample's posterior mode was this individual"
}

#here we can look at the entire posterior of true genotypes for this individual
#Note, individuals with samples strongly linked to them will have precisely
#estimated true genotypes while individuals without samples strongly linked
#will have very imprecisely estimated true genotypes. If no samples ever allocated,
#you are just drawing true genotypes from the estimated population-level genotype frequencies
out <- t(apply(G.samps[ind,,],2,FUN=map.genos,unique.genos))
head(out,10)

loci <- 1
loci <- loci+1
table(out[,loci])