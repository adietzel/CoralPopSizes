# Calculating proportions of crest, slope and flat habitat dominated by live coral
# Code author: Andreas Dietzel
# Last modified: February 2, 2020

rm(list = ls())

### Load packages
library(tidyverse)
library(rstanarm)
library(rstan)
library(brms)
library(ggridges)
library(tidybayes)
library(coda)

# Load files
load("RData/Input_Files.RData")

# Calculate proportion live coral dominated habitats for each map
HPLC <- HabitatMaps %>%
  bind_rows(., LOF) %>%
  filter(!HABITAT == "Exclude") %>% # categories such as "mangrove" or "open ocean water" excluded
  group_by(LOCATION, HABITAT, M_COVER) %>%
  summarise(Agg.Area = sum(AREA), .groups = "keep") %>%
  group_by(LOCATION) %>%
  mutate(Prop = Agg.Area / sum(Agg.Area)) %>%
  filter(M_COVER == "Coral") %>%
  filter(!(HABITAT == "Others" & Prop == 1)) %>%
  ungroup() %>%
  mutate(HABITAT = as.factor(HABITAT))

# Set priors and fit model
prior = c(set_prior("normal(0, 100)", class = "b"),
          set_prior("normal(0, 100)", class = "Intercept"),
          set_prior("cauchy(0, 5)", class = "sigma"))

brm.logit <- brm(logit_scaled(sqrt(sqrt(Prop))) ~ HABITAT, data = HPLC,
                 prior = prior,
                 chains = 4, iter = 17500, warmup = 5000, thin = 5)

# Extract model fit
newdata = data.frame(HABITAT = levels(HPLC$HABITAT))
Xmat <- model.matrix(~ HABITAT, data = newdata)
coefs <- as.matrix(as.data.frame(rstan:::extract(brm.logit$fit)))[, 1:4]
inv.logit <- binomial()$linkinv
fit = inv.logit(coefs %*% t(Xmat)) ^ 4 # back-transform
newdata = cbind(newdata, plyr:::adply(fit, 2, function(x) {
  data.frame(Mean = mean(x), Median = median(x), HPDinterval(as.mcmc(x)))
}))

# Plot model fit
ED_Fig_5a <- ggplot(newdata, aes(y = Median, x = HABITAT), colour = "blue") +
  geom_boxplot(data = HPLC, aes(x = HABITAT,y = Prop)) +
  geom_pointrange(aes(ymin = lower, ymax = upper, x = as.numeric(HABITAT)),
                  colour = "blue") +
  theme_classic() + coord_flip() +
  labs(y = "Proportion dominated by live coral", x = "")

# Extract posterior disttributions
Post <- as.data.frame(fit)
colnames(Post) <- levels(HPLC$HABITAT)
Post$Total <- rowSums(Post)
HPLC.posterior <- Post %>% gather(HABITAT, HPLC)

# Calculate colony densities by habitat at each site
DensHab <- Interc %>%
  group_by(SiteCode, Depth) %>%
  dplyr:::summarise(mean.diam = n() * pi / (2 * sum(1 / intercept)),
                    PCcover = sum(intercept) / 100,
                    density.m2 = (n() / .1) /
                      (mean.diam / (100 * 1000)) / (1000 ^ 2),
                    .groups = "keep") %>%
  ungroup() %>%
  mutate(Depth2 = ifelse(Depth == "flat", "Flat/Back Reef",
                        ifelse(Depth == "slope", "Slope/Fore Reef", "Crest")))

# Boxplot colony density by habitat
ED_Fig_5b <- ggplot(DensHab, aes(x = Depth2, y = density.m2)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "", y = "Colonies per m2") +
  ylim(0, 80) + coord_flip() + theme_classic()

# Relative colony densities in each habitat type
MeanDensHab <- DensHab %>%
  group_by(Depth) %>%
  summarise(MeanDens = mean(density.m2),
            .groups = "keep") %>%
  ungroup() %>%
  bind_rows(., data.frame(Depth = "others",
                          MeanDens = mean(MeanDensHab$MeanDens))) %>%
  mutate(RelMeanDens = MeanDens/mean(MeanDens))

### Rescale habitat proportions by mean relative colony density
Post_rescale <- Post
Post_rescale[,"Crest"] <- Post_rescale[, "Crest"] *
  MeanDensHab$RelMeanDens[MeanDensHab$Depth == "crest"]
Post_rescale[,"Flat/Back Reef"] <- Post_rescale[, "Flat/Back Reef"] *
  MeanDensHab$RelMeanDens[MeanDensHab$Depth == "flat"]
Post_rescale[,"Slope/Fore Reef"] <- Post_rescale[, "Slope/Fore Reef"] *
  MeanDensHab$RelMeanDens[MeanDensHab$Depth == "slope"]
Post_rescale[,"Total"] <- rowSums(Post_rescale[, 1:4])

# Combine original and rescaled posterior distributions in one file
HPLC.posterior_rescale <- Post_rescale %>%
  gather(HABITAT, HPLC) %>%
  mutate(type = "rescaled")

HPLC.posterior$type <- "not rescaled"

HPLC.comp <- bind_rows(HPLC.posterior_rescale, HPLC.posterior)

# Boxplot showing original and rescaled posterior distributions
ED_Fig_5c <- ggplot(HPLC.comp, aes(x = HPLC, y = HABITAT, fill = type)) +
    geom_density_ridges(alpha = .5) +
    labs(x = "Proportion of total reef with\nlive coral as dominant benthic cover",
         y = "Posterior density", fill = "") +
    theme_classic() + theme(legend.position = "top")

# Multi-panel figure
ED_Fig_5 <- cowplot::plot_grid(ED_Fig_5a, ED_Fig_5b, ED_Fig_5c,
                               rel_heights = c(.5,.4, .8),
                   ncol = 1, labels = c("a", "b", "c"), align = "v")

# ggsave("figures/ED_Figure_5_HPLC.pdf", width = 10, height = 15, units = "cm")

# Save output files
save(HPLC.posterior_rescale, HabMaps.sel,
     file = "RData/OUT1_Habitat_proportions.RData")

