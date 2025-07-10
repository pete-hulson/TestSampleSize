#########################
# Testing effective sample size computation for bootstrapping
#########################

install.packages("dirmult")
library(dirmult)
Norm = function(Vec) Vec / sum(Vec)

Rmultinom = function(n, size, prob){
  Return = NULL
  if(length(size)==1 & is.vector(prob)) Return = rmultinom(n=n, size=size, prob=prob)
  if(length(size)>=2 & is.vector(prob)) for(i in 1:n) Return = cbind(Return, rmultinom(n=1, size=size[i], prob=prob))
  if(length(size)==1 & is.matrix(prob)) for(i in 1:n) Return = cbind(Return, rmultinom(n=1, size=size, prob=prob[,i]))
  if(length(size)>=2 & is.matrix(prob)) for(i in 1:n) Return = cbind(Return, rmultinom(n=1, size=size[i], prob=prob[,i]))
  return(Return)
}

SampFn = function(Vec){
  Long = rep.int(1:length(Vec), times=Vec)
  Samp = sample(Long, replace=TRUE)
  New = tabulate(Samp, nbins=length(Vec))
  return(New)
}

Prop = rdirichlet(1, alpha=rep(1,10))
#Prop = Norm(matrix(N_atsm[,Nyears,StrataI,Mode_i][,1],nrow=1))

#####################
# Check effective sample size theory
#####################
# effective sample size of each raw data row
Size = 100
N = 1e2
Rand = t( Rmultinom(N, size=Size, prob=Prop[1,]) )
#Rand = CompTrue
Rand = Rand / Size
Neff = apply(Rand, MARGIN=1, FUN=function(Boot){sum( Boot * (1-Boot)) / sum( (Boot-Prop)^2 )})
mean(Neff)
1 / mean(1/Neff)    # Should be equal to Size

#####################
# Check boostrap for:
# 1. rows (tows) and samples within rows (individuals)          -- BIASED (double counts variance)
# 2. Aggregating across tows, and then samples within aggregate -- FINE
# 3. Just rows (tows)                                           -- FINE (MY RECOMMENDATION)
# 4. Just samples within rows (individuals)                     -- FINE
#####################
# effective sample size for mean across rows
Size = 10
N = 100
Nsims = 100
Save = matrix(NA, nrow=Nsims, ncol=8)
for(SimI in 1:nrow(Save)){
  # two-stage bootstrap
  Prop = rdirichlet(1, alpha=rep(1,10))
  Rand = t( Rmultinom(N, size=Size, prob=Prop[1,]) )
  aOrig = Norm(colSums(Rand))
  Neff = aBoot = NULL
  for(BootI in 1:Nsims){
    WhichRow = sample(1:nrow(Rand), replace=TRUE)
    CompBoot = Rand[WhichRow,]
    CompBoot = t(apply(CompBoot, MARGIN=1, FUN=SampFn))
    aBoot = rbind(aBoot, Norm(colSums(CompBoot)) )
    Neff = c(Neff, sum( aBoot[BootI,] * (1-aBoot[BootI,])) / sum( (aBoot[BootI,]-aOrig)^2 ) )
  }
  Save[SimI,1] = mean(Neff)
  Save[SimI,2] = 1 / mean( 1/Neff )
  
  # one-stage "aggregated" bootstrap
  Neff = aBoot = NULL
  for(BootI in 1:Nsims){
    CompBoot = SampFn(colSums(Rand))
    aBoot = rbind(aBoot, Norm(CompBoot) )
    Neff = c(Neff, sum( aBoot[BootI,] * (1-aBoot[BootI,])) / sum( (aBoot[BootI,]-aOrig)^2 ) )
  }
  Save[SimI,3] = mean(Neff)
  Save[SimI,4] = 1 / mean( 1/Neff )
  
  # one-stage "data-level" bootstrap
  Prop = rdirichlet(1, alpha=rep(1,10))
  Rand = t( Rmultinom(N, size=Size, prob=Prop[1,]) )
  aOrig = Norm(colSums(Rand))
  Neff = aBoot = NULL
  for(BootI in 1:Nsims){
    CompBoot = t(apply(Rand, MARGIN=1, FUN=SampFn))
    aBoot = rbind(aBoot, Norm(colSums(CompBoot)) )
    Neff = c(Neff, sum( aBoot[BootI,] * (1-aBoot[BootI,])) / sum( (aBoot[BootI,]-aOrig)^2 ) )
  }
  Save[SimI,5] = mean(Neff)
  Save[SimI,6] = 1 / mean( 1/Neff )
  
  # one-stage "tow-level" bootstrap
  Prop = rdirichlet(1, alpha=rep(1,10))
  Rand = t( Rmultinom(N, size=Size, prob=Prop[1,]) )
  aOrig = Norm(colSums(Rand))
  Neff = aBoot = NULL
  for(BootI in 1:Nsims){
    WhichRow = sample(1:nrow(Rand), replace=TRUE)
    CompBoot = Rand[WhichRow,]
    aBoot = rbind(aBoot, Norm(colSums(CompBoot)) )
    Neff = c(Neff, sum( aBoot[BootI,] * (1-aBoot[BootI,])) / sum( (aBoot[BootI,]-aOrig)^2 ) )
  }
  Save[SimI,7] = mean(Neff)
  Save[SimI,8] = 1 / mean( 1/Neff )
}
colMeans(Save)



# pete's test ----

# effective sample size for single (tow) and double (tow-subsample) bootstrapping

# set up sampling variables
H = 100 # number of 'hauls'
Nsims = 100 # number of sim replicates
R = 25 # number of replicates of sim
test = NULL
for(r in 1:R){
  Save = matrix(NA, nrow=Nsims, ncol=6)
  
  for(SimI in 1:nrow(Save)){

    # set up pop'n comp
    Prop = rdirichlet(1, alpha=rep(1,10))
    
    # number of fish in haul
    N = round(runif(H, 10, 1000))
    # number of fish in subsample
    n = round(runif(H, 1, 10))
    
    # stage 1 sample
    s1 = NULL
    for(h in 1:H) s1 = rbind(s1, t(Rmultinom(1, size = N[h], prob = Prop[1,])))

    # stage 2 sample
    s2 = NULL
    for(h in 1:H) s2 = rbind(s2, t(Rmultinom(1, size = n[h], prob = Norm(s1[h,]))))
    Prop_surv = Norm(colSums(s2))
    
    # 'survey' effective sample size
    Neff_surv = sum(Prop_surv * (1 - Prop_surv)) / sum((Prop_surv - Prop)^2)
    Save[SimI, 1] = Neff_surv
    
    # test bootstrap
    Neff1 = aBoot1 = NULL
    Neff2 = aBoo21 = NULL
    for(BootI in 1:Nsims){
      # sample hauls
      WhichRow = sample(1:nrow(s2), replace = TRUE)
      # one-stage "tow-level" bootstrap
      CompBoot1 = s2[WhichRow,]
      aBoot1 = Norm(colSums(CompBoot1))
      Neff1 = c(Neff1, sum(aBoot1 * (1 - aBoot1)) / sum((aBoot1 - Prop_surv)^2))
      # two-stage bootstrap
      CompBoot2 = t(apply(CompBoot1, MARGIN = 1, FUN = SampFn))
      aBoot2 = Norm(colSums(CompBoot2))
      Neff2 = c(Neff2, sum(aBoot2 * (1 - aBoot2)) / sum((aBoot2 - Prop_surv)^2))
    }
    Save[SimI,2] = mean(Neff1)
    Save[SimI,3] = 1 / mean(1/Neff1)
    Save[SimI,4] = mean(Neff2)
    Save[SimI,5] = 1 / mean(1/Neff2)
    
    Save[SimI,6] = sum(n)

    
  }
  test = rbind(test, colMeans(Save))
}
colnames(test) <- c("Neff_survey", "Neff_tow", "Neff_tow_1", "Neff_tow_subsample", "Neff_tow_subsample_1", "n_total")
test
colMeans(test)
# observations:
# 1. tow-level bootstrap Neff is slightly larger than 'true' survey Neff (Neff_tow > Neff_survey)
# 2. tow-subsample bootstrap Neff is smaller than 'true' survey Neff (Neff_tow_subsample < Neff_survey)
