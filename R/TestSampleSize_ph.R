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








##### testing

# effective sample size for mean across rows

# set up sampling variables
H = 100 # number of hauls=
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
test
colMeans(test)





## to do:
# do a comparison in surveyISS between full bootstrap and haul only bootstrap
# check out simulation to see if resampling haul causes issues






























# effective sample size for mean across rows

# set up sampling variables
N = 500 # number of pop'n units/schools
H = 100 # number of hauls
SubSize = 10 # size of subsample
Nsims = 100 # number of sim replicates

R = 25
test = NULL
for(r in 1:R){
  Save = matrix(NA, nrow=Nsims, ncol=4)
  
  for(SimI in 1:nrow(Save)){
    
    # 'true' prop in pop'n units/schools
    Prop = rdirichlet(N, alpha=rep(1,10))
    # 'true' pop'n comp
    aTrue = colMeans(Prop)
    
    # pop'n units sampled in hauls
    SampProp = Prop[sample(1:N, H, replace = TRUE),]
    colMeans(SampProp)
    # get haul sample
    Rand = NULL
    for(i in 1:H) Rand = rbind(Rand, t(rmultinom(1, SubSize, SampProp[i,])))
    # 'survey' comps
    aSurv = Norm(colSums(Rand))
    
    # effective sample size of 'survey'
    Neff_surv = sum(aSurv * (1 - aSurv)) / sum((aSurv - aTrue)^2)
    Save[SimI, 1] = Neff_surv
    
    # resample hauls (to be consistent across tests)
    WhichRow = sample(1:nrow(Rand), replace = TRUE)
    
    # two-stage bootstrap
    Neff = aBoot = NULL
    for(BootI in 1:Nsims){
      CompBoot = Rand[WhichRow,]
      CompBoot = t(apply(CompBoot, MARGIN=1, FUN=SampFn))
      aBoot = rbind(aBoot, Norm(colSums(CompBoot)))
      Neff = c(Neff, sum( aBoot[BootI,] * (1-aBoot[BootI,])) / sum( (aBoot[BootI,]-aSurv)^2 ) )
    }
    Save[SimI,2] = mean(Neff)
    
    # one-stage "tow-level" bootstrap
    Neff = aBoot = NULL
    for(BootI in 1:Nsims){
      CompBoot = Rand[WhichRow,]
      aBoot = rbind(aBoot, Norm(colSums(CompBoot)))
      Neff = c(Neff, sum( aBoot[BootI,] * (1-aBoot[BootI,])) / sum( (aBoot[BootI,]-aSurv)^2 ) )
    }
    Save[SimI,3] = mean(Neff)
    
    # two-stage: (1) resample tows, (2) multinomial draw based on school proportions
    Neff = aBoot = NULL
    for(BootI in 1:Nsims){
      SampProp2 = SampProp[WhichRow,]
      Rand2 = NULL
      for(i in 1:H) Rand2 = rbind(Rand2, t(rmultinom(1, SubSize, SampProp2[i,])))
      aBoot = rbind(aBoot, Norm(colSums(Rand2)))
      Neff = c(Neff, sum( aBoot[BootI,] * (1-aBoot[BootI,])) / sum( (aBoot[BootI,]-aSurv)^2 ) )
    }
    Save[SimI,4] = mean(Neff)
    
  }
  test = rbind(test, colMeans(Save))
}
test

colMeans(test)








# effective sample size for mean across rows
Size = 10
N = 100
Nsims = 100
Save = matrix(NA, nrow=Nsims, ncol=8)
for(SimI in 1:nrow(Save)){
  
  # set up pop'n
  Prop = rdirichlet(1, alpha=rep(1,10))
  
  # set up sample
  Rand = t(Rmultinom(N, size = Size, prob = Prop[1,]))
  aOrig = Norm(colSums(Rand))
  
  # two-stage bootstrap
  Neff = aBoot = NULL
  for(BootI in 1:Nsims){
    WhichRow = sample(1:nrow(Rand), replace=TRUE)
    CompBoot = Rand[WhichRow,]
    CompBoot = t(apply(CompBoot, MARGIN=1, FUN=SampFn))
    aBoot = rbind(aBoot, Norm(colSums(CompBoot)) )
    Neff = c(Neff, sum(aBoot[BootI,] * (1 - aBoot[BootI,])) / sum((aBoot[BootI,] - aOrig)^2))
  }
  Save[SimI,1] = mean(Neff)
  Save[SimI,2] = 1 / mean(1 / Neff)
  
  # two-stage: 1st sample hauls, 2nd multi draw of Prop
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





























































# another test

# set up comp
Prop = rdirichlet(1, alpha=rep(1,10))
# set up pop'n unit
N = 1000
PopUnit = t(rmultinom(1, N, Prop[1,]))
# we know that, on average, the Neff_PopUnit will equal N
Neff_PopUnit = sum(Norm(PopUnit) * (1 - Norm(PopUnit))) / sum((Norm(PopUnit) - Prop[1,])^2 )

# now subsample pop'n unit
n = 100 # sample size
H = 100 # number of sample replicates (i.e., H for hauls)

Sub_PopUnit = t(rmultinom(H, n, Norm(PopUnit)))
Neff_Sub = NULL
for(i in 1:H) Neff_Sub = c(Neff_Sub, sum(Norm(Sub_PopUnit[i,]) * (1 - Norm(Sub_PopUnit[i,]))) / sum((Norm(Sub_PopUnit[i,]) - Norm(PopUnit))^2 ))
# we also know that, on average, Neff_Sub will equal n
mean(Neff_Sub)
1 / mean(1 / Neff_Sub)

# now, bootstrap the sample
BootSub_PopUnit = t(apply(Sub_PopUnit, MARGIN=1, FUN=SampFn))
Neff_BootSub = NULL
for(i in 1:H) Neff_BootSub = c(Neff_BootSub, sum(Norm(BootSub_PopUnit[i,]) * (1 - Norm(BootSub_PopUnit[i,]))) / sum((Norm(BootSub_PopUnit[i,]) - Norm(PopUnit))^2 ))

samp = colSums(BootSub_PopUnit)
Neff_BootSub_t = NULL
for(i in 1:H) Neff_BootSub_t = c(Neff_BootSub_t, sum(Norm(BootSub_PopUnit[i,]) * (1 - Norm(BootSub_PopUnit[i,]))) / sum((Norm(BootSub_PopUnit[i,]) - Norm(samp))^2 ))

mean(Neff_BootSub)
1 / mean(1 / Neff_BootSub)

mean(Neff_BootSub_t)
1 / mean(1 / Neff_BootSub_t)

sum(BootSub_PopUnit)





SampFn_2 = function(Vec, scalar){
  Long = rep.int(1:length(Vec), times=Vec)
  Samp = sample(Long, round(scalar * length(Long)), replace=TRUE)
  New = tabulate(Samp, nbins=length(Vec))
  return(New)
}



BStest = function(H, n, p, rep, scalar){
  
  Sub = t(rmultinom(H, n, p))
  
  Neff = NULL
  for(r in 1:rep){
    # bootstrap without scalar
    BootSub = t(apply(Sub, MARGIN=1, FUN=SampFn))
    samp = colSums(BootSub)
    Neff_r = NULL
    for(i in 1:H) Neff_r = c(Neff_r, sum(Norm(BootSub[i,]) * (1 - Norm(BootSub[i,]))) / sum((Norm(BootSub[i,]) - Norm(samp))^2 ))
    
    
  BootSub = t(apply(Sub, MARGIN=1, FUN=SampFn_2))
  samp = colSums(BootSub)
  Neff_r = NULL
  for(i in 1:H) Neff_r = c(Neff_r, sum(Norm(BootSub[i,]) * (1 - Norm(BootSub[i,]))) / sum((Norm(BootSub[i,]) - Norm(samp))^2 ))
  Neff = c(Neff, mean(Neff_r))
  }

  out = c(n, mean(Neff))
  return(out)
}

H = 100 # number of sample replicates (i.e., H for hauls)
p = rdirichlet(1, alpha=rep(1,10)) [1,]
rep = 1000 # number of bootstrap replicates
n = seq(from = 10, to = 1000, by = 10) # sub-sample sizes
test = NULL

for(i in 1:length(n)) test = rbind(test, BStest(H, n[i], p, rep))

plot(test[,1], test[,1] / test[,2], ylim = c(0.5, 1.5))
abline(h=1)
abline(h=mean(test[,1] / test[,2]), lty = 2)








CompSub = t(apply(Rand, MARGIN=1, FUN=SampFn))



Neff_PopUnit

Norm(PopUnit)



