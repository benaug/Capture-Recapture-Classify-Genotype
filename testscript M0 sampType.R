#This version allows genotyping error rates to vary by categorical sample covariates
#continuous covariates can be done in practice, but very slow-need a theta matrix for
#every sample, rather than one for each category level

#In genoSPIM, the genotyping error settings below work great, and they appear to work well
#here, too, but with more uncertainty in ID. However, I expect you will have a bad time
#if you don't have a good number of higher quality samples when using low quality samples
#with low amplification rates and high genotyping error rates. Using spatial information
#helps disproportionately more in these scenarios.

#Don't try to add time or individual detection effects in M0 scripts
#unless you make corresponding changes to custom updates

library(coda)
library(nimble)
source("sim.data.sampType.R")
source("init.data.R")
source("build.genos.R")
source("map.genos.R")
source("NimbleModel M0 sampType.R")
source("Nimble Functions M0 sampType.R")

#must run this line for nimble to do what we want
nimbleOptions(determinePredictiveNodesInModel = FALSE)

#We will use the fisher data set parameter estimates from Augustine et al. (2020)
#to specify the genetic portion of the simulation settings

#This is a list of locus-level genotypes. Each list element corresponds to a locus
#and contains a matrix of all possible combinations of locus-level genotypes, one per row.
#the order of alleles in a genotype does not matter
load("unique.genos.RData")
str(unique.genos)
#A real data set should enumerate the observed loci-level genotypes in the order listed here,
#e.g. 152-152 is 1, 152-154 is 2, 152-156 is 3, etc.
unique.genos[[1]] #1st loci

n.levels <- unlist(lapply(unique.genos,nrow)) #how many loci-level genotypes per locus?

#This function creates objects to determine which classifications are 1) correct, 2) false allele, and 3) allelic dropout
built.genos <- build.genos(unique.genos)
ptype <- built.genos$ptype #list of length n.loci with each element being an n.levels[m] x n.levels[m] matrix containing 
#indicator for error type for true genotypes along rows and classified genotype along columns
#note, "ptype" in list form is used in the data simulator, but ptypeArray is a ragged array used in nimble

#####example 1: true homozygote:
#the first genotype at locus 1 is 152-152. Because it is homozygous, there cannot be allelic dropout.
unique.genos[[1]][1,]
#Now, looking at ptype for true locus-level genotype 1 at locus 1:
#we see only a correct classification or a false allele are possible.
ptype[[1]][1,]

#there is only 1 correct classification event, indicated with a 1:
unique.genos[[1]][ptype[[1]][1,]==1,]

#here are the classifications that are the result of false alleles (everything that is not correct)
unique.genos[[1]][ptype[[1]][1,]==3,]

#####example 2: true heterozygote:
#the second genotype at locus 1 is 152-154. Because it is heterozygous, there can be allelic dropout or false alleles
unique.genos[[1]][2,]
#Now, looking at ptype for true locus-level genotype 2 at locus 1:
#we see correct classification, allelic dropout, or false allele are possible.
ptype[[1]][2,]

#there is only 1 correct classification event, indicated with a 1:
unique.genos[[1]][ptype[[1]][2,]==1,]

#here are the classifications that are the result of allelic dropout (152 or 154 can drop out)
unique.genos[[1]][ptype[[1]][2,]==2,]

#here are the classifications that are the result of false alleles
unique.genos[[1]][ptype[[1]][2,]==3,]

#OK, moving on to the genotype frequencies.
load("gammameans.RData")#loci-level genotype frequency estimates for fisher data set
#These are the frequencies of each locus-level genotype at each locus
str(gammameans)

#Now we have the information required to simulate a data set similar to the fisher data set

#First, let's decide how many loci to use. This repo assumes you have at least 2 (I didn't put in work to allow 1, can be done in theory)
#with 9 loci, there is rarely uncertainty in ID, unless genotyping error is high. Can use fewer loci
#here to get more uncertainty in sample IDs
n.loci <- 9

#discard unused information if you don't use them all
if(n.loci!=9){
  for(i in 9:(n.loci+1)){
    gammameans[[i]] <- NULL
    unique.genos[[i]] <- NULL
    ptype[[i]] <- NULL
  }
}
n.levels <- unlist(lapply(unique.genos,nrow)) #update n.levels in case some loci discarded
#now all lists of length "n.loci"
str(gammameans)
str(unique.genos)
str(ptype)
n.levels

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
  # gamma[[i]] <- rep(1/n.levels[i],n.levels[i]) #This simulates equal genotype frequencies
  gamma[[i]] <- gammameans[[i]] #This uses the frequencies estimated from fisher data set
}

#Genotype observation process parameters. Can have failed amplification (missing completely at random) and genotyping error
#using estimates from fisher data set below
samp.levels <- 2 #number of sample type covariates. Each type has it's own genotyping error rates.
#p.amp below is for each sample type in this script instead of each loci as in others
p.amp <- c(0.999,0.25) #sample by replication amplification probabilities (controls level of missing scores in G.obs)
p.geno.het <- vector("list",samp.levels)
p.geno.hom <- vector("list",samp.levels)
#P(correct, allelic dropout,false allele) for heterozygotes (using fisher ests here)
p.geno.het[[1]] <- c(0.806,0.185,0.009) 
p.geno.het[[2]] <- c(0.489,0.496,0.015)
#P(correct,false allele) for homozygotes
p.geno.hom[[1]] <- c(0.994,0.006)
p.geno.hom[[2]] <- c(0.999,0.001)
pi.samp.type <- c(0.52,0.48) #frequencies of each sample type

data <- sim.data.sampType(N=N,p.y=p.y,lambda.y=lambda.y,K=K,#cap-recap parameters/constants
                  n.loci=n.loci,p.amp=p.amp,n.rep=n.rep,
                  p.geno.hom=p.geno.hom,p.geno.het=p.geno.het,
                  gamma=gamma,IDcovs=IDcovs,ptype=ptype,
                  pi.samp.type=pi.samp.type)

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
M <- 150
if(M<N)stop("M must be larger than simulate N")

#set some gamma inits. Using equal across locus-level genotypes here
#note, gamma is a ragged matrix for use in nimble.
gammaMat <- matrix(0,nrow=n.loci,ncol=max(n.levels))
for(l in 1:n.loci){
  gammaMat[l,1:n.levels[l]] <- rep(1/n.levels[l],n.levels[l])
}

#provide some ballpark inits for gamma to help initialize data
inits <- list(gammaMat=gammaMat) #using gammaMat inits to initialize G.true not determined by observed data
nimbuild <- init.data(data=data,M=M,inits=inits,initTrue=FALSE) #can initialize simulated data sets at true IDs for testing

n.samples <- data$n.samples

#ptype is a ragged array that tells use which genotype classifications are correct, allelic dropout, or false alleles
ptype <- built.genos$ptypeArray

#OK, what is the data we use to fit the model?
#1) data$G.obs, the observed genotypes
#2) y.true/ID. Note, you can build y.true with ID, the individual of each sample, and this.k, the occasion of each sample.
#we specify a distribution on y.true, not ID. If we initialize ID, y.true is initialized.
#If we update ID, y.true is updated. So, we just initialize these, provide to nimble
#as inits, and update them on each MCMC iteration

#supply data to nimble
Nimdata <- list(G.obs=nimbuild$G.obs,samp.type=data$samp.type) #adding sample covariates here

#inits for nimble
#Without spatial information, the possiblity of false alleles makes convergence difficult with sparse G.obs
#Helps to provide ballpark inits for p.geno.het especially, providing them for p.geno.hom here, too.
#essentially, if you initialize the false allele probability near 1, the samples will be allocated very poorly
#on the first iterations and the false allele probability can get stuck near 1.
p.geno.het.init <- matrix(NA,nrow=2,ncol=3)
p.geno.hom.init <- matrix(NA,nrow=2,ncol=2)
p.geno.het.init[1,] <- c(0.9,0.05,0.05)
p.geno.het.init[2,] <- c(0.9,0.05,0.05)
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
                n.levels=n.levels,max.levels=max(n.levels),ptype=ptype,samp.levels=samp.levels)

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
                                                  na.ind=nimbuild$G.obs.NA.indicator,n.levels=n.levels),
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
#                     control = list(i = i,m=m,n.levels=n.levels,n.rep=n.rep,samp.type=data$samp.type,
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
                control = list(M=M,n.loci=n.loci,n.levels=n.levels,n.rep=n.rep,samp.type=data$samp.type,
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
Cmcmc$run(5000,reset=FALSE) #can extend run by rerunning this line, e.g. run longer if not converged, or want more samples
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