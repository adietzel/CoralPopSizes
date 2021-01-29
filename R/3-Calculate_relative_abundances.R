################################################################################
# FitUnifiedModel.R
# Code Author: Sean Connolly with modifications by Andreas Dietzel
# ----------------------------------------------------------------------
# Calculate mean abundances of observed species
#===============================================================================

library(tidyverse)
library(poilog)
library(cowplot)
library(quantreg)
library(stats)
library(data.table)

load("RData/Input_Files.RData")

### Compare numerical abundance and total area occupied
Interc %>%
  group_by(Species) %>%
  summarise(N = n(),
            Cover = sum(intercept),
            .groups = "keep") %>%
  ggplot(aes(x = N, y = Cover)) +
  geom_point() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_smooth(method = "lm") +
  labs(x = "total number of intercepts", y = "total intercept length (cm)")

ggsave("figures/FigSXX_CoverNumeric.pdf", width = 4.5, height = 4)

AbundHab.Site <- Interc %>%
  group_by(site, Species) %>%
  summarise(abund = n(), .groups = "keep") %>%
  ungroup() %>%
  complete(site, Species, fill = list(abund = 0)) %>%
  mutate(Site.Nr = ((as.numeric(substr(site, 2, 2)) - 1) * 4) +
           as.numeric(substr(site, 3, 3))) %>%
  left_join(., SiteCode_LookUp[, 2:3], by = c("site" = "site_new"))

abund.wide <- AbundHab.Site %>%
  dplyr::select(Species, Site.Nr, RegHab, abund) %>%
  group_by(RegHab) %>%
  spread(Site.Nr, abund) %>% ungroup() %>%
  mutate(rowSum = rowSums(.[3:14])) %>%
  filter(rowSum > 0)

abundmat.list <- lapply(seq_along(unique(abund.wide$RegHab)), function (i) {

  RegHabCode <- unique(abund.wide$RegHab)[i]

  abundmat <- abund.wide %>%
    filter(RegHab == RegHabCode) %>%
    dplyr::select(-Species, -rowSum, -RegHab) %>%
    as.matrix(header=T)
})

names(abundmat.list) <- unique(abund.wide$RegHab)

# read in the functions used to fit the NB-LN-PL model
source("R/4-Functions_UnifiedModel.R")

# Fit the unified model to each local community
# ---------------------------------------------------------------------------

modeloutput <- data.frame()
rcom.df <- data.frame()

for (i in 1:length(abundmat.list)) {

  abundmat <- abundmat.list[[i]]
  # Find a set of reasonable initial par values to be passed onto nlminb
  # ---------------------------------------------------------------------------
  # For the power-law parameters, regress log(mean) against log(variance)
  mus <- rowMeans(abundmat)
  vars <- apply(abundmat,1,var)
  mvfit <- lm(log(vars)~log(mus))
  plini <- c(coef(mvfit)[1],log(coef(mvfit)[2]))
  names(plini) = c("log(a)","log(b)")
  # intercept (coefficient 1) is log(a) and coefficient 2 is b
  # (var = a*mean^b)

  # For the metacommunity SAD, fit the poisson lognormal to the summed
  # metacommunity abundance.
  metafit = poilogMLE(rowSums(abundmat))
  sitem = exp(metafit$par[1]-log(ncol(abundmat)))
  # rescale mu to site scale, and backtransform to arithmetic scale

  lnini = c(sqrt(sitem),1/metafit$par[2])
  names(lnini) = c("sqrt(m)","1/sig")
  # These transformations have been chosen because they produce relatively
  # symmetric profile likelihoods across the different metacommunities analysed
  # in the paper.

  init <- c(plini,lnini)
  # In practice, it is important to re-fit the model using different initial
  # parameter estimates, to ensure that true MLEs are found

  # Fit the NB-LN-PL model
  # ---------------------------
  fit.nblnpl <- nlminb(init,Fit.nbv.ln,lower=c(-Inf,-Inf,0,0),obs=abundmat)
  conv = fit.nblnpl$conv
  while(conv!=0) {
  # if(fit.nblnpl$conv!=0) {
    # in case of convergence issue, re-start from similar ini pars
    fit.nblnpl <- nlminb(rnorm(length(init),mean=init,sd=abs(init/4)),Fit.nbv.ln,obs=abundmat)
    conv = fit.nblnpl$conv
  }
  print(c("convergence ",fit.nblnpl$conv))
   # Be sure to check here if the fit converged. If it didn't, rerun the previous
   # step.

  # Back-transform parameter estimates. Note that the third parameter is the
  # mean of log-abundance for the SAD, BACK-TRANSFORMED to an arithmetic scale
  mle.nblnpl <- c(exp(fit.nblnpl$par[1:2]),fit.nblnpl$par[3]^2,1/abs(fit.nblnpl$par[4]))
  names(mle.nblnpl) = c("a","b","m","sig")
  mll.nblnpl <- -fit.nblnpl$obj
  rcom <- colSums(abundmat)/mean(colSums(abundmat))
    # a vector of relative community sizes
  rcom.df <- rbind(rcom.df, rcom)
  p0.nblnpl <- getP0.nbv.ln(mle.nblnpl[1],mle.nblnpl[2],mle.nblnpl[3],mle.nblnpl[4],rcom)
    # probability of observing 0 abundance at every site
  Spool.nblnpl <- nrow(abundmat)/(1-p0.nblnpl)
    # estimated community species-pool size
  output <- data.frame(a = mle.nblnpl[1],
                       b = mle.nblnpl[2],
                       m = mle.nblnpl[3],
                       sig = mle.nblnpl[4],
                       SpeciesPool = Spool.nblnpl)
  modeloutput <- rbind(modeloutput, output)
}
modeloutput$RegHab <- names(abundmat.list)

################################################################################
################################################################################

# Simulate new species by site abundance matrices for each community using
# the parameters m, sig, a and b
# ---------------------------------------------------------------------------

Bootstrap.Abund <- data.frame()
CommSize.all <- data.frame()

for (i in 1:length(modeloutput$RegHab)) {
  numsim <- 100000

  # Read in the data and fit the unified model to it
  a.mle <- modeloutput$a[i] # power law parameter a
  b.mle <- modeloutput$b[i] # power law parameter b
  m.mle <- modeloutput$m[i] # mean of lognormal on arithmetic scale
  sig.mle <- modeloutput$sig[i] # variance of lognormal distribution
  Spool.nblnpl <- modeloutput$SpeciesPool[i] # size of species pool
  rcom <- as.numeric(rcom.df[i,]) # vector of relative community sizes
  Sobs <- nrow(abundmat.list[[i]]) # number of observed species

  # number of sites in the sample
  nsites <- ncol(abundmat.list[[i]]) # number of sites

  # simulate data
  set.seed(1)
  simdata <- array(dim = c(round(Spool.nblnpl), 2, numsim))
  CommSize.df <- data.frame(RegHab = rep(modeloutput$RegHab[i], times = numsim),
                            CommSize = NA)

  for(isim in 1 : numsim) {
    # randomly draw a mean abundance using the fitted m and s
    mu.i <- rlnorm(round(Spool.nblnpl), meanlog = log(m.mle), sdlog = sig.mle)
    CommSize.df[isim,2] <- sum(mu.i)
    # using the observed relative sample sizes across sites of this group,
    # calculate species-and-site-specific mean abundance values
    mus.ij <- sapply(mu.i,'*',rcom)
    for(isp in 1:length(mu.i)) {
      # for each species i,
      # if var > mean (i.e. a*mu_i^b > mu_i --> mu_i > exp(log(a)/(1-b))
      # randomly draw an abundance value at each site j from a NB with
      # the sampled, adjusted mean mus.ij and the fitted aggregation parameter
      # k = mu_i^2/(a*mu_i^b-mu_i)
      if(mu.i[isp] > exp(log(a.mle) / (1 - b.mle))) {
        k.i <- mu.i[isp] ^ 2 / (a.mle * mu.i[isp] ^ b.mle - mu.i[isp])
        simdata[isp,,isim] <- c(mean(rnbinom(nsites, size = k.i, mu = mus.ij[,isp])), mu.i[isp])
      } else
        # if var <= mean (i.e. a*mu_i^b <= mu_i --> mu_i <= exp(log(a)/(1-b))
        # randomly draw an abundance value at each site j from a Poisson with
        # the sampled adjusted mean mus.ij
        simdata[isp,,isim] <- c(mean(rpois(nsites,lambda = mus.ij[, isp])), mu.i[isp])
    }
  }

  Abund <- rowSums(as.data.frame(abundmat.list[[i]])) / 12
  MaxAbund <- max(Abund) + 5

  CommSize.all <- bind_rows(CommSize.all, CommSize.df)

  simdata.long <- data.table(as.data.frame.table(simdata))
  simdata.long <- data.table::dcast(simdata.long, Var1 + Var3 ~ Var2, value.var = "Freq")
  colnames(simdata.long) <- c("Species", "Sim", "Abund_Obs", "Abund_True")

  simdata.long <- as.data.frame(simdata.long)
  simdata.long1 <- simdata.long[simdata.long$Abund_Obs < MaxAbund,]
  simdata.long2 <- simdata.long[(simdata.long$Abund_Obs %in% Abund),]

  simdata.long3 <- simdata.long2 %>%
    group_by(Abund_Obs) %>%
    top_n(500, Sim)

  simdata.long3$RegHab <- modeloutput$RegHab[i]

  Bootstrap.Abund <- bind_rows(Bootstrap.Abund, simdata.long3)
}

CommSize.all$Sim <- unique(simdata.long$Sim) # assign simulation run ID

### SAVING OUTPUT FILES
Bootstrap.RelAbund <- Bootstrap.Abund %>%
  left_join(.,CommSize.all, by = c("RegHab" = "RegHab", "Sim" = "Sim")) %>%
  mutate(RelAbund_True = Abund_True/CommSize)

Freq <- Bootstrap.Abund %>%
  group_by(RegHab, Abund_Obs) %>%
  summarise(count = n(), .groups = "keep")

MeanCommSize.all <- CommSize.all %>%
  group_by(RegHab) %>%
  summarise(MeanCommSize = mean(CommSize), .groups = "keep")

UnifiedModel.params <- modeloutput %>%
  left_join(., MeanCommSize.all, by = "RegHab")

save(Bootstrap.RelAbund, UnifiedModel.params,
     file = "RData/OUT3_Relative_Species_Abundances.RData")

