#### ANOVA on IUCN status

rm(list = ls())

# Load libraries and files ####
library(tidyverse)
library(rstan)
library(rstanarm)
library(brms)
library(broom)
library(coda)

load("RData/Input_Files.RData")
load("RData/OUT6_PopSizes.RData")

PopSizes$IUCNstatus.bin <- ifelse(
  PopSizes$IUCNstatus %in% c("Endangered","Vulnerable"), "Elevated",
  ifelse(PopSizes$IUCNstatus %in% c("Least Concern", "Near Threatened"), "Low",
         PopSizes$IUCNstatus))

PopSizesIUCN <- PopSizes %>%
  mutate(Log10PopSize = log10(PopSizeMean * 10 ^ 9)) %>%
  filter(IUCNstatus.bin %in% c("Elevated","Low")) %>%
  mutate(IUCN.regroup = as.factor(ifelse(IUCNstatus %in% c("Vulnerable","Endangered"),
                                         "Elevated", as.character(IUCNstatus)))) %>%
  filter(IUCN.regroup %in% c("Near Threatened", "Least Concern", "Elevated"))

### Model fitting and comparison

IUCN.stanGaus <- stan_glm(Log10PopSize ~ IUCN.regroup, data = PopSizesIUCN,
                          family = gaussian,
                          prior_intercept = normal(0, 100),
                          prior = normal(0, 100),
                          prior_aux = cauchy(0, 5),
                          chains = 4, iter = 17500, thin = 5, warmup = 5000)

PopSizesIUCN <- PopSizesIUCN %>%
  dplyr::select(Species, Log10PopSize, IUCN.regroup)

# get the coefficients of the model runs
coefs <- as.data.frame(IUCN.stanGaus) %>% as.matrix()
coefs <- coefs[,1:3]

(cmat2 <- cbind("Elevated vs Low"= c(1,-.5,-.5),
                "Elevated vs Least Concern"= c(1,-1,0),
                "Elevated vs Near Threatened"= c(1,0,-1),
                "Near Threatened vs Least Concern"= c(0,-1,1)))

# Add intercept column which is 0 and make contrast matrix computer-readible

(cmat2 <- cbind(0, t(cmat2) %*% contr.treatment(3)))
fit <- coefs %*% t(cmat2)
fit <- 10 ^ fit

# Calculate highest posterior density intervals
newdata = data.frame(term = rownames(cmat2))
newdata = cbind(newdata, plyr:::adply(fit, 2, function(x) {
  data.frame(Mean = mean(x), Median = median(x), HPDinterval(as.mcmc(x), prob = .99))
}))

# Plot contrasts
(IUCN.effectsize.plot <- ggplot(newdata, aes(y = Mean, x = term)) +
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_discrete("") + ylab("effect size") +
    coord_flip() +
    theme(axis.text = element_text(size = 12)) +
    scale_y_continuous(breaks = c(0.25,0.5,0.75,1,1.25)))

### R2 value
rsq <- bayes_R2(IUCN.stanGaus)
hist(rsq)
R2.mcmc <- HPDinterval(as.mcmc(rsq))

### Summary plots
cols <- c("non-coral species" = "blue", "Least Concern" = "#60c659",
          "Near Threatened" = "#cce226", "Vulnerable" = "#f9e814",
          "Endangered" = "#FC7F3F", "Data Deficient" = "black",
          "Central and Eastern Indo-Pacific"="#22A884",
          "Species examined"="#7AD152")

SampleSize.IUCN <- PopSizes %>% group_by(IUCNstatus) %>%
  summarise(Freq = n(), .groups = "keep")

SampleSize.IUCN.bin <- PopSizes %>% group_by(IUCNstatus.bin) %>%
  summarise(Freq = n(), .groups = "keep")

### Plotting

Box.IUCN <- ggplot(PopSizes,
                   aes(x = reorder(IUCNstatus, desc(IUCNstatus)),
                       y = PopSizeMean, fill = IUCNstatus)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(trans = "log10", breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1","1", "10","100"),
                     limits=c(0.0001, 700)) +
  geom_text(data = SampleSize.IUCN,
            aes(x = reorder(IUCNstatus, desc(IUCNstatus)),
                label = paste0("n=", Freq), y = 120), hjust = 0) +
  xlab("") + ylab("population size in billions") +
  scale_fill_manual(name = "IUCN status", values = cols, na.value = "grey") +
  theme(axis.text = element_text(size = 12), legend.position = "") +
  coord_flip()

Box.IUCN.grouped <- ggplot(PopSizes[PopSizes$IUCNstatus.bin %in% c("Elevated","Low"),],
                    aes(x = IUCNstatus.bin, y = PopSizeMean), fill = white) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(trans = "log10", breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1","1", "10","100"),
                     limits = c(0.0001, 700)) +
  geom_text(data = SampleSize.IUCN.bin[
    SampleSize.IUCN.bin$IUCNstatus.bin %in% c("Elevated","Low"),],
    aes(x = reorder(IUCNstatus.bin, desc(IUCNstatus.bin)),
        label = paste0("n=", Freq), y = 120), hjust = 0) +
  xlab("") + ylab("population size in billions") +
  scale_fill_manual(name = "IUCN status", values = cols, na.value = "grey") +
  theme(axis.text = element_text(size = 12), legend.position = "") +
  coord_flip()

pdf("figures/FigS5_IUCNoverview.pdf", width = 6, height = 5)

(IUCN.main.plot <- cowplot::plot_grid(Box.IUCN, Box.IUCN.grouped, IUCN.effectsize.plot,
                                      ncol = 1, align = "v", labels = c("A","B","C"),
                                      rel_heights = c(.8, .45, .65)))

