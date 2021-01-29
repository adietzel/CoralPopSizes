# Calculate combined proportion of unobserved species

# Load packages and data
library(tidyverse)
library(poilog)
library(data.table)

load("RData/Input_Files.RData")
load("RData/OUT3_Relative_Species_Abundances.RData")

# read in the functions used to fit the NB-LN-PL model
source("R/4-Functions_UnifiedModel.R")

modeloutput <- UnifiedModel.params

AbundHab.Site <- Interc %>%
  group_by(site, Species) %>%
  summarise(abund = n()) %>%
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

Bootstrap.Abund <- data.frame()
CommSize.all <- data.frame()

numsim <- 1000

for (i in 1:length(modeloutput$RegHab)) {

  # Read in the data and fit the unified model to it
  a.mle <- modeloutput$a[i] # power law parameter a
  b.mle <- modeloutput$b[i] # power law parameter b
  m.mle <- modeloutput$m[i] # mean of lognormal on arithmetic scale
  sig.mle <- modeloutput$sig[i] # variance of lognormal distribution
  Spool.nblnpl <- modeloutput$SpeciesPool[i] # size of species pool
  abundmat <- abundmat.list[[i]]
  rcom <- colSums(abundmat)/mean(colSums(abundmat))
  Sobs <- nrow(abundmat.list[[i]]) # number of observed species

  # number of sites in the sample
  nsites <- ncol(abundmat.list[[i]]) # number of sites

  # simulate data
  set.seed(1)
  simdata <- array(dim = c(round(Spool.nblnpl), 2, numsim))
  CommSize.df <- data.frame(RegHab = rep(modeloutput$RegHab[i],times = numsim),
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

  CommSize.all <- bind_rows(CommSize.all, CommSize.df)

  simdata.long <- data.table(as.data.frame.table(simdata))
  simdata.long <- data.table::dcast(simdata.long, Var1 + Var3 ~ Var2, value.var = "Freq")
  colnames(simdata.long) <- c("Species", "Sim", "Abund_Obs", "Abund_True")
  simdata.long <- as.data.frame(simdata.long)
  simdata.long$RegHab <- modeloutput$RegHab[i]

  Bootstrap.Abund <- bind_rows(Bootstrap.Abund, simdata.long)
}

CommSize.all$Sim <- unique(simdata.long$Sim)

Bootstrap.RelAbund <- Bootstrap.Abund %>%
  left_join(.,CommSize.all, by = c("RegHab"="RegHab", "Sim"="Sim")) %>%
  mutate(RelAbund_True = Abund_True/CommSize) %>%
  filter(Abund_Obs == 0) %>%
  group_by(RegHab, Sim) %>%
  summarise(Prop_Unobs = sum(RelAbund_True)*100) %>%
  group_by(RegHab) %>%
  summarise(Mean_Prop_Unobs = round(mean(Prop_Unobs), digits = 1),
            SD_Prop_Unobs = round(sd(Prop_Unobs), digits = 1))
