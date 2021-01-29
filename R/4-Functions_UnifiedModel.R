################################################################################
# NBLNPL_LLfxn.R
# Code Author: Mizue Hisano, with edits by Sean Connolly
# Most recent changes: 24 January 2016 (tidy up for DRYAD)
# ----------------------------------------------------------------------
# Functions called by FitUnifiedModel.R to fit the unified model to a sites 
# (columns) x species (rows) matrix of abundances (counts of individuals).
# 
# If you use this code, please cite the following:
# Connolly et al. (2017) A unified model explains commonness and
# rarity on coral reefs. Ecology Letters.
#===============================================================================


Fit.nbv.ln = function(pars, obs) {
  # parameters are:
  # a,b: Power-law parameters for how spatial variance in abundance scales
  # with the mean ()var = a*mu^b)
  # log(m), sig: mean and sd of the natural log of metacommunity species
  # abundance
  # obs is a matrix with as many columns as sites and as many rows as species
  a = exp(pars[1])
  b = exp(pars[2])
  m = pars[3]^2
  sig = 1/abs(pars[4])
  
  # number of sites:
  S = ncol(obs)

  rcom = colSums(obs)/mean(colSums(obs))  # Relative community size
  # Note rcom is scaled so that it has a mean value of 1 (this necessary
  # so that the calculated variance is based on site-level mean abundance --
  # specifically the species mean abundance at average community size.

  # first calculate p0, the probability of not observing any species at any sites
  p0 = getP0.nbv.ln(a,b,m,sig,rcom)
  LL = rep(NA,nrow(obs))
  for (iSpecies in 1:nrow(obs)) {
    # calculate the likelihood for one species at a time
    # vector of abundances observed for this species
    ab = obs[iSpecies,]
    loglikeli = LL.nbv.ln(a,b,m,sig,ab,rcom)
    # add normalization factor since we condition on the species being
    # observed in at least one site:
    LL[iSpecies] = loglikeli - log(1-p0)
  }

  # returning negative log-likelihood
  return(-sum(LL))
}

getP0.nbv.ln = function(a, b, m, sig, rcom) {
  # Function to calculate probability of observing a vector of all zero 
  # abundances, if the probability distribution of abundances follows the
  # unified model

  # a,b = power-law parameters
  # log(m), sig = parameters of the lognormal distribution
  # rcom = vector of relative community sizes (average is 1)

  # number of sites
  S = length(rcom)

  thresh = exp(log(a)/(1-b)) # obtained by solving a*mu^b = mu for mu. 
                             # This is where variance=mean, and is the cutoff 
                             # for the transition between poisson vs negative
                             # binomial spatial distribution of abundances
   # if b>1 use poisson below threshold & neg bin above threshold
  if(b>1) {
    int.up = integrate(fnbv,lower=thresh,upper=Inf,a=a,b=b,m=m,sig=sig,
              ab=rep(0,S),rcom=rcom)$value
    int.low = integrate(fpoiv,lower=0,upper=thresh,m=m,sig=sig,
              ab=rep(0,S),rcom=rcom,subdivisions=2000)$value
  } else {
  # if b<1 use nb below threshold & pois above threshold
    int.up = integrate(fnbv,lower=0,upper=thresh,a=a,b=b,m=m,sig=sig,
              ab=rep(0,S),rcom=rcom,subdivisions=2000)$value
    int.low = integrate(fpoiv,lower=thresh,upper=Inf,m=m,sig=sig,
              ab=rep(0,S),rcom=rcom,subdivisions=2000)$value
    if(int.up==0&int.low==0) {
      int.up = integrate(fnbv,lower=0,upper=Inf,a=a,b=b,m=m,sig=sig,
              ab=rep(0,S),rcom=rcom)$value
    }
  }
  return(int.up + int.low)
}

LL.nbv.ln = function(a,b,m,sig,ab,rcom) {
  # Function to calculate the log-likelihood of observing abundance vector "ab"
  # for a species, given parameters a,b,m,sig
  # NOTE that this function does not include the zero-vector truncation, so
  # is not called directly from an optimizer

  # a,b power-law scaling parameters
  # log(m),sig: parameters of the lognormal distribution
  # ab: vector of abundances observed for this species
  # rcom: vector of relative community sizes
  
  # Note that this function may be problematic for data sets for which the MLE
  # of b is <1 or very close to one, and thus the else{} statement is executed.

  thresh = exp(log(a)/(1-b)) # obtained by solving a*mu^b = mu for mu
  # if b>1, use pois below 'thresh' and nb above 'thresh'
  # if b<1, use nb below 'thresh' and pois above 'thresh'

  if(b>1) {
    # use NB above threshold of var = mean
    part2a = integrate(fnbv,lower=thresh,upper=Inf,
              a=a,b=b,m=m,sig=sig,ab=ab,rcom=rcom)$value
    # use Pois below threshold of var=mean
    part2b = integrate(fpoiv,lower=0,upper=thresh,
              m=m,sig=sig,ab=ab,rcom=rcom,subdivisions=2000)$value
  } else {
    # use NB below threshold of var = mean
    part2a = try(integrate(fnbv,lower=0,upper=thresh,
              a=a,b=b,m=m,sig=sig,ab=ab,rcom=rcom)$value,T)
    # use Pois above threshold of var=mean
    part2b = integrate(fpoiv,lower=thresh,upper=Inf,
              m=m,sig=sig,ab=ab,rcom=rcom,subdivisions=2000)$value
    if(class(part2a)=="try-error"||(part2a==0&part2b==0)) {
      part2a = try(integrate(fnbv,lower=0,upper=Inf,
              a=a,b=b,m=m,sig=sig,ab=ab,rcom=rcom)$value,T)
      if(class(part2a)!="try-error")
        part2a = part2a - part2b else {
        # the integral doesn't work - it is probably divergent.
        # pass an arbitrary very small value to loglikeli to avoid these pars
        part2a = .Machine$double.xmin
      }
    }
  }
  # calculate the total loglikelihood
  # NOTE that parts 2a,b must be added before logarithm can be taken
  loglikeli = log(part2a + part2b)
  return(loglikeli)
}

fnbv = function(mu, a, b, m, sig, ab, rcom) {
  # Note this function is integrated only for mean abundances mu for which
  # the variance (a*mu^b) exceeds the mean.

  # a,b = power-law parameters
  # mu = mean for species i on the arithmetic scale (average across sites)
  # log(m), sig = parameters of the lognormal distribution
  # ab = vector of abundances
  # rcom = vector of relative community sizes (average is 1)

  v = a*mu^b   # note variance is assumed constant across sites,
               # based on estimated species' mean abundance at
               # the average local community size
  k = mu^2/(v-mu)   # likewise, k constant across sites
  mus = sapply(rcom,'*',mu)   # Different means for each site:
    # number of rows equals length of mu vector passed by integrate
    # number of columns equals to the number of sites

  if(length(rcom)>1)
    nbpart = colSums(sapply(1:length(mu),function(i)
              dnbinom(ab,size=k[i],mu=mus[i,],log=T))) else
    nbpart = sapply(1:length(mu),function(i)
              dnbinom(ab,size=k[i],mu=mus[i,],log=T))

  # return the probability
  exp(nbpart + dlnorm(mu,log(m),sig,log=T))
}

fpoiv = function(mu, m, sig, ab,rcom) {
  # This is used whenever the power law would predict
  # variance less than the mean -- thus we implicitly
  # assume that variance=mean for any mean abundance
  # below the threshold where mu < exp( log(a)/(1-b) )

  # mu = mean for species i on the arithmetic scale (average across sites)
  # log(m), sig = parameters of the lognormal distribution
  # ab = vector of abundances
  # rcom = vector of relative community sizes (average is 1)

  # let the value of mu vary according to the relative sample size of each site
  # but keep the average at mu

  mus = sapply(rcom,'*',mu)

  if(length(rcom)>1)
    poispart = colSums(apply(mus,1,function(i)dpois(ab,i,log=T))) else
    poispart = apply(mus,1,function(i)dpois(ab,i,log=T))

  # return the probability
  exp(poispart + dlnorm(mu,log(m),sig,log=T))
}
