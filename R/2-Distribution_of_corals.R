# Calculate total number of corals in study domain
# Code author: Andreas Dietzel
# Last modified: February 2, 2020

rm(list = ls())

# Load libraries
library(tidyverse)
library(sp)
library(raster)
library(rstan)
library(rstanarm)
library(gstat)

### Load files
load("RData/Input_Files.RData")
load("RData/OUT1_Habitat_proportions.RData")

### Calculate coral cover and colony densities at survey sites
# Colony densities were calculated using method by Marsh et al. 1984
CoverDens <- Interc %>%
  group_by(SiteCode) %>%
  dplyr:::summarise(mean.diam = n()*pi / (2*sum(1/intercept)),
                    PCcover = sum(intercept)/30000, # total transect length = 30,000cm or 300m
                    nrInt = n(),
                    density.m2 = (nrInt/300) / ((mean.diam/(100)))) %>%
  dplyr::select(-mean.diam)

### Fit generalised additive model: colony density ~ coral cover
CoverDens.GAM.Bayes <- stan_gamm4(density.m2 ~ s(PCcover, k = 3),
                                  data = CoverDens, family = gaussian,
                                  prior_intercept = normal(0, 10), prior = normal(0, 10),
                                  prior_aux = cauchy(0, 5),
                                  chains = 4, iter = 17500, thin = 5, warmup = 5000,
                                  adapt_delta = .9995)

# Plotting gam fit
newdata <- with(CoverDens, data.frame(PCcover = seq(0, 1, len = 1000)))
fit <- posterior_linpred(CoverDens.GAM.Bayes, newdata = newdata)
newdata <- cbind(newdata, tidyMCMC(fit, conf.int = T, conf.level = .95,
                                   conf.method = "HPDinterval"))

GAM.pooled <- ggplot(newdata, aes(y = estimate, x = PCcover * 100)) +
  geom_line() +
  geom_point(data = CoverDens, aes(x = PCcover * 100, y = density.m2)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "blue", alpha = .3) +
  xlab("Coral cover (%)") + ylab("Colonies per m2") +
  xlim(0, 60) + ylim(0, 60)

### Interpolate coral cover across study domain
# Exclude coral cover values outside study domain (bbx)
CoralCover@data$CEIP <- sapply(over(CoralCover, geometry(bbx),
                                    returnList = T),length)

CoralCover <- CoralCover[CoralCover@data$CEIP == 1,]

# Interpolate coral cover using inverse distance weighted model
idw.CC <- idw(formula = CoralCover@data[,7] ~ 1,
              locations = CoralCover,
              newdata = grd.IP, idp = 2, nmax = 10)

# Rasterize output, exclude grid cells with reef area = 0 and
# convert to spatial pixels data frame
idw.r.CC <- raster(idw.CC)
idw.r.CC.Reef <- overlay(idw.r.CC, reefs.rast.CEIP,
                         fun = function(x,y) {x * (y/y)})
idw.r.CC.Reef.SpPix <- as(idw.r.CC.Reef, "SpatialPixelsDataFrame")

# Predict colony density from coral cover data in each raster cell using GAM
newdata.raw <- data.frame(
  PCcover = round(idw.r.CC.Reef@data@values/100, digits = 3),
  CellNo = 1:length(idw.r.CC.Reef@data@values))

newdata <- data.frame(PCcover = unique(newdata.raw$PCcover[
  !is.na(newdata.raw$PCcover)]))

fit <- posterior_linpred(CoverDens.GAM.Bayes, newdata = newdata, transform = T)

newdata <- cbind(newdata,
                 tidyMCMC(fit, conf.int = T,
                          conf.level = .95, conf.method = "HPDinterval"))

newdata.all <- newdata.raw %>%
  left_join(., newdata, by = "PCcover")

# Predict number of corals in each grid cell (mean, lower bound, upper bound)

# Mean
RasterCells.Dens.Mean <- idw.r.CC
RasterCells.Dens.Mean@data@values <- newdata.all$estimate
RasterCells.NC.Mean <- overlay(RasterCells.Dens.Mean,
                               reefs.rast.CEIP,
                               fun = function(x,y) {x * (y/100) * (res^2)})

NC.Scl.Mean <- round(cellStats(RasterCells.NC.Mean, "sum"), digits = 1) / (10^9)

# Lower bound
RasterCells.Dens.Lower <- idw.r.CC
RasterCells.Dens.Lower@data@values <- newdata.all$conf.low
RasterCells.NC.Lower <- overlay(RasterCells.Dens.Lower,
                                reefs.rast.CEIP,
                                fun = function(x,y) {x * (y/100) * (res^2)})

NC.Scl.Lower <- round(cellStats(RasterCells.NC.Lower, "sum"), digits = 1) / (10^9)

# Upper bound
RasterCells.Dens.Upper <- idw.r.CC
RasterCells.Dens.Upper@data@values <- newdata.all$conf.high
RasterCells.NC.Upper <- overlay(RasterCells.Dens.Upper,
                                reefs.rast.CEIP,
                                fun=function(x,y) {x * (y/100) * (res^2)})

NC.Scl.Upper <- round(cellStats(RasterCells.NC.Upper, "sum"), digits = 1) / (10^9)

# Summarize
RasterCells.NC <- data.frame(CellNo = 1:length(RasterCells.NC.Lower@data@values),
                             Lower = RasterCells.NC.Lower@data@values,
                             Mean = RasterCells.NC.Mean@data@values,
                             Upper = RasterCells.NC.Upper@data@values)

# Calculate reef area in each grid cell
ReefArea <- sum((reefs.rast.CEIP@data@values/100) * (res/1000)^2, na.rm = T)

NumberCorals <- data.frame(NC.Scl.Mean.uncorr = NC.Scl.Mean,
                           NC.Scl.Lower.uncorr = NC.Scl.Lower,
                           NC.Scl.Upper.uncorr = NC.Scl.Upper,
                           ReefArea.CEIP = ReefArea)

# Calculate total number of corals in study domain
# Use proportion of reef dominated by live coral to correct for reef area overestimate
NumberCorals.CEIP = sample_n(NumberCorals, 10000, replace = T) %>%
  mutate(HPLC.all2 = HPLC.posterior_rescale$HPLC[
    HPLC.posterior_rescale$HABITAT == "Total"]) %>%
  rowwise() %>%
  mutate(NC.Scl =
           rnorm(n = 1, mean = NC.Scl.Mean.uncorr,
                 sd = (NC.Scl.Upper.uncorr - NC.Scl.Lower.uncorr) / (2*1.96)) * HPLC.all2) %>%
  ungroup() %>%
  summarise(ReefArea = mean(ReefArea),
              NC.Scl.Mean = mean(NC.Scl),
              NC.Scl.Lower = quantile(NC.Scl,probs = .025),
              NC.Scl.Upper = quantile(NC.Scl,probs = .975)) %>%
  mutate(MeanCC = mean(CoralCover@data$HARD_CORAL/100),
         SdCC = sd(CoralCover@data$HARD_CORAL/100),
         type = "CEIP", Percent = 100, Location = "CEIP")

### Plotting
# Histogram of coral cover data from Bruno et al. 2016
Histo.Bruno <- ggplot(CoralCover@data, aes(HARD_CORAL)) +
  geom_histogram(breaks = seq(0,100,7), colour = "black", fill = "white") +
  geom_vline(xintercept = median(CoralCover@data$HARD_CORAL)) +
  geom_text(aes(x = 45, y = 135,
                label = round(median(CoralCover@data$HARD_CORAL, na.rm = T),
                              digits = 1))) +
  xlab("Coral cover (%)") + ylab("Frequency") +
  scale_x_continuous(limits = c(0, 90), breaks = c(0, 25, 50, 75))

# Histogram of coral cover from species abundance survey sites
Histo.Hughes <- ggplot(CoverDens, aes(PCcover*100)) +
  geom_histogram(breaks = seq(0, 100, 7), colour = "black", fill = "white") +
  geom_vline(xintercept = median(CoverDens$PCcover*100)) +
  xlab("Coral cover (%)") + ylab("Frequency") +
  geom_vline(xintercept = median(CoverDens$PCcover*100, na.rm = T)) +
  geom_text(aes(x = 29, y = 15,
                label = round(median(CoverDens$PCcover*100, na.rm = T),
                              digits = 1))) +
  scale_x_continuous(limits = c(0, 90), breaks = c(0, 25, 50, 75))

# Histogram of colony density estimates at survey sites
Histo.DensHughes <- ggplot(CoverDens, aes(density.m2)) +
  geom_histogram(bins = 10, breaks = seq(0, 50, 4), colour = "black", fill = "white") +
  xlab("Colonies per m2") + ylab("Frequency") +
  geom_vline(xintercept = median(CoverDens$density.m2, na.rm = T)) +
  geom_text(aes(x = 29, y = 12,
                label = round(median(CoverDens$density.m2, na.rm = T),
                              digits = 1))) +
  scale_x_continuous(breaks = c(0, 20, 40))

# Histogram of predicted colony densities in grid cells
Histo.DensPred <- newdata.all %>%
  filter(!is.na(estimate)) %>%
  ggplot(aes(estimate)) +
  geom_histogram(breaks = seq(0, 50, 3),
                 colour = "black", fill = "white") +
  geom_vline(xintercept = median(newdata.all$estimate, na.rm = T)) +
  geom_text(aes(x = 35, y = 490,
                label = round(median(newdata.all$estimate, na.rm = T),
                              digits = 1))) +
  scale_x_continuous(breaks = c(0, 20, 40)) +
  xlab("Colonies per m2") + ylab("Frequency")

# Histogram of interpolated coral cover data in each grid cell
Histo.Cover <- newdata.all %>%
  filter(!is.na(PCcover)) %>%
  ggplot(aes(PCcover * 100)) +
  geom_histogram(breaks = seq(0, 100, 7),
                 colour = "black", fill = "white") +
  geom_vline(xintercept = median(newdata.all$PCcover*100, na.rm = T)) +
  geom_text(aes(x = 45, y = 420,
                label = median(newdata.all$PCcover * 100, na.rm = T))) +
  xlab("Coral cover (%)") + ylab("Frequency") +
  scale_x_continuous(limits = c(0, 90), breaks = c(0, 25, 50, 75))

# Extended Data Figure 4: Model inputs and outputs for predicting coral abundances
cowplot::plot_grid(Histo.Hughes, Histo.DensHughes,
                   Histo.Bruno, GAM.pooled,
                   Histo.Cover, Histo.DensPred,
                   ncol = 2, align = "hv",
                   labels = c("a","b","c","d","e","f"))

ggsave("figures/ED_Figure_4_GAM.pdf", width = 6, height = 7)

# SAVE OUTPUT FILES
save(NumberCorals, RasterCells.NC,
     file = "RData/OUT2_Number_Corals.RData")

