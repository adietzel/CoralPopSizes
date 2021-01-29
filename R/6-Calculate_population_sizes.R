#### Calculate population sizes
# Code author: Andreas Dietzel
# Last modified: February 2, 2020

rm(list = ls())

# Load libraries and files ####
library(tidyverse)
library(sp)

load("RData/Input_Files.RData")
load("RData/OUT1_Habitat_proportions.RData")
load("RData/OUT2_Number_Corals.RData")
load("RData/OUT3_Relative_Species_Abundances.RData")
load("RData/OUT5_NrCoralsRange.RData")

# Load tree population size estimates from ter Steege et al. (2015)
TreesAmazon <- read.csv("data/Appendix_terSteege2015.csv",
                        strip.white = T, stringsAsFactors = T) %>%
  arrange(desc(Estimated.population)) %>%
  mutate(RelRank = row_number() / n()) %>%
  filter(!is.na(Estimated.population))

### IUCN color scheme
cols <- c("Least Concern" = "#60c659",
          "Near Threatened" = "#cce226",
          "Vulnerable" = "#f9e814",
          "Endangered" = "#FC7F3F",
          "Data Deficient" = "black",
          "Central and Eastern Indo-Pacific"="#22A884",
          "Species examined"="#7AD152")

# Posterior distributions of live-coral dominated crest, flat and slope
HPLC.crest <- HPLC.posterior_rescale$HPLC[
  HPLC.posterior_rescale$HABITAT == "Crest"]
HPLC.flat <- HPLC.posterior_rescale$HPLC[
  HPLC.posterior_rescale$HABITAT == "Flat/Back Reef"]
HPLC.slope <- HPLC.posterior_rescale$HPLC[
  HPLC.posterior_rescale$HABITAT == "Slope/Fore Reef"]

# Calculate observed abundances of each species in each region and habitat
# Join simulated true abundances by observed abundance per region and habitat
RelAbund <- Interc %>%
  group_by(Species, Region, Depth, RegHab) %>%
  summarise(Abund_Obs = n() / 12, .groups = "keep") %>%
  full_join(., Bootstrap.RelAbund[, -1],
            by = c("RegHab" = "RegHab", "Abund_Obs" = "Abund_Obs")) %>%
  filter(!is.na(Species))

### Summarize unified model parameters and community size
UnifiedModel.overview <- RelAbund %>%
  group_by(Region, Depth, RegHab, Species) %>%
  summarise(MeanRA = mean(RelAbund_True),
            .groups = "keep") %>%
  group_by(Region, Depth, RegHab) %>%
  summarise(Sum_RelAbund = sum(MeanRA),
            NrSp.Recorded = length(unique(Species)),
            .groups = "keep") %>%
  left_join(UnifiedModel.params, by = "RegHab") %>%
  mutate(Prop.Species = NrSp.Recorded/SpeciesPool) %>%
  dplyr::select(Region, Depth, m, sig, a, b, NrSp.Recorded, SpeciesPool,
                Prop.Species, MeanCommSize, Sum_RelAbund, RegHab) %>%
  arrange(RegHab) #%>%
  # write.csv("figures/UnifiedModelParameters.csv")

# Calculate population sizes for each simulation run, species and region
NC_Species.bootstrap <- data.frame()

nsim <- 10000

for (i in 1:length(NrCorals.Range$Species)) {

  # subset species and region
  Species = NrCorals.Range$Species[i]
  Region = NrCorals.Range$Region[i]

  # Subset range of relative abundances on crest, flat and slope
  RelAbund.crest = RelAbund$RelAbund_True[
    RelAbund$Species == Species & RelAbund$Region == Region & RelAbund$Depth == "crest"]
  RelAbund.flat = RelAbund$RelAbund_True[
    RelAbund$Species == Species & RelAbund$Region == Region & RelAbund$Depth == "flat"]
  RelAbund.slope = RelAbund$RelAbund_True[
    RelAbund$Species == Species & RelAbund$Region == Region & RelAbund$Depth == "slope"]

  # Calculate number of corals
  x = data.frame(
    # Randow draw from distribution of corals within region adjusted by species richness
    NC.SR = rnorm(
    n = nsim,
    mean = NrCorals.Range$Mean.NC.SR[i],
    sd = (NrCorals.Range$Upper.NC.SR[i] - NrCorals.Range$Lower.NC.SR[i]) / (2 * 1.96)),

    # Draw relative true abundances in habitats for each sim run
    RAT.crest = ifelse(length(RelAbund.crest) > 0,
                       sample(RelAbund.crest, size = nsim, replace = T) *
                         sample(HPLC.crest, nsim, replace = T), 0),
    RAT.flat = ifelse(length(RelAbund.flat) > 0,
                      sample(RelAbund.flat, size = nsim, replace = T) *
                        sample(HPLC.flat, nsim, replace = T), 0),
    RAT.slope = ifelse(length(RelAbund.slope) > 0,
                       sample(RelAbund.slope, size = nsim, replace = T) *
                         sample(HPLC.slope, nsim, replace = T), 0)) %>%
    mutate(NC.Species = NC.SR * RAT.crest + NC.SR * RAT.flat + NC.SR * RAT.slope)

  # Summarise and append to output data frame
  NC_Species.bootstrap <- bind_rows(
    NC_Species.bootstrap,
    data.frame(Species = Species, Region = Region,
               NC = x$NC.Species, Sim = 1:nsim))
}

# Calculate pop size for each species across regions and habitats
NC_Species.bootstrap2 <- NC_Species.bootstrap %>%
  group_by(Species,Sim) %>%
  summarise(NC.range = sum(NC))

# Save intermediate output
save(NC_Species.bootstrap2, NC_Species.bootstrap,
     file = "RData/OUT6_NC_BS.RData")

# # Load intermediate output
# load("RData/OUT6_NC_BS.RData")

# Calculate population size by region
PopSizes.Region <- NC_Species.bootstrap %>%
  group_by(Species,Region) %>%
  summarise(PopSizeMean = mean(NC) / 10 ^ 9,
            PopSizeLower = quantile(NC, .025) / 10 ^ 9,
            PopSizeUpper = quantile(NC, .975) / 10 ^ 9,
            .groups = "keep")

# Look-up table for regions
(Region.LookUp <- data.frame(Region = unique(PopSizes.Region$Region),
                             RegionCode = c(4,5,1,2,3)) %>%
    arrange(RegionCode))

# Find species in dataset that have population size = 0 in regions within their range
PopSizes.Region0 <- PopSizes.Region %>%
  filter(PopSizeMean == 0) %>%
  left_join(., Region.LookUp, by = "Region") %>%
  left_join(., PresAb, by = c("Species" = "Species",
                              "RegionCode" = "RegionCode")) %>%
  filter(PresAb.Either == 1) %>%
  mutate(RelRank = 0) %>%
  dplyr::select(Species, Region, PopSizeMean, RelRank)

# Calculate relative abundance ranks for each region and attach species with pop size = 0
PopSizes.Region <- PopSizes.Region %>%
  arrange(desc(PopSizeMean)) %>%
  filter(PopSizeMean > 0) %>%
  group_by(Region) %>%
  mutate(RankRegion = row_number(),
         RelRank = 1 - row_number() / n()) %>%
  ungroup() %>%
  bind_rows(., PopSizes.Region0) %>%
  group_by(Region) %>%
  mutate(RankRegion = ifelse(RelRank == 0, max(RankRegion, na.rm = T),
                             RankRegion))

# Number corals for comparison
(NrCorals.CEIP <- round(NumberCorals$NC.Scl.Mean[NumberCorals$Location == "CEIP"],
                        digits = 1))

# Calculate global population sizes as sums of regional pop sizes
PopSizes <- PopSizes.Region %>%
  group_by(Species) %>%
  summarise(PopSizeMean = sum(PopSizeMean, na.rm = T),
            PopSizeLower = sum(PopSizeLower, na.rm = T),
            PopSizeUpper = sum(PopSizeUpper, na.rm = T),
            MeanRelRankRegion = mean(RelRank),
            .groups = "keep") %>%
  ungroup() %>%
  left_join(., SupportData, by = "Species") %>%
  arrange(desc(PopSizeMean)) %>%
  mutate(CumPSMean = cumsum(PopSizeMean),
         Rel.CumPSMean = CumPSMean/max(CumPSMean),
         Rel.CumPSMean.CEIP = CumPSMean/NrCorals.CEIP,
         Rank = row_number()) %>%
  arrange(PopSizeMean) %>%
  mutate(CumPSMean.Inv = cumsum(PopSizeMean),
         Hyperdom = ifelse(Rel.CumPSMean < 0.5, "Hyperdominant", "Non-hyperdominant"))

# Set IUCN categories
levels(PopSizes$IUCNstatus) <- c("Least Concern", "Near Threatened", "Vulnerable",
                                 "Endangered", "Data Deficient")

# Calculate global abundance rank of each species
RelRank.Global <- PopSizes %>%
  arrange(desc(PopSizeMean)) %>%
  mutate(GlobalRelRank = 1 - row_number() / n(),
         PopSizeGlobal = PopSizeMean) %>%
  dplyr::select(Species, GlobalRelRank, PopSizeGlobal, Hyperdom)

# Join global rank to regional rank file
PopSizes.all <- PopSizes.Region %>%
  left_join(., RelRank.Global,by = "Species") %>%
  left_join(., SupportData, by = "Species") %>% ungroup()

# Set IUCN categories
levels(PopSizes.all$IUCNstatus) <-
  c("Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Data Deficient")

PopSizes.all$Region <- factor(PopSizes.all$Region,
                              levels(PopSizes.all$Region)[c(3,4,5,1,2)])

PopSizes.byRegion.Global <- PopSizes.all %>%
  mutate(Region = "Global") %>%
  bind_rows(., PopSizes.all) %>%
  mutate(Region = as.factor(Region))

PopSizes.byRegion.Global$Region <- factor(
  PopSizes.byRegion.Global$Region,
  levels(PopSizes.byRegion.Global$Region)[c(3,4,5,6,1,2)])

### Plotting

# Set font size
font_size = 12

# Set colours for corals and trees
cols2 <- c("Indo-Pacific coral species" = "#0070C0",
           "humans" = "#008CDA",
           "Amazonian tree species" = "#60c659")

# Calculate rank abundance distribution for corals and tree and combine
Coral.RAD <- PopSizes %>%
  arrange(desc(PopSizeMean)) %>%
  mutate(Rank = row_number(), RelRank = Rank / max(Rank),
         PopSize = PopSizeMean * 10 ^ 9,
         type = "Indo-Pacific coral species",
         hyperdom = ifelse(Rel.CumPSMean < 0.5, "hyperdominant", "non-hyperdominant"))

Tree.RAD <- TreesAmazon %>%
  arrange(desc(Estimated.population)) %>%
  mutate(Rank = row_number(), RelRank = Rank / max(Rank),
         PopSize = Estimated.population,
         type = "Amazonian tree species",
         hyperdom = ifelse(Rank < 228, "hyperdominant", "non-hyperdominant"))

Compare.RAD <- Tree.RAD[, c("PopSize","RelRank","type","hyperdom")] %>%
  bind_rows(., Coral.RAD[, c("PopSize","RelRank","type","hyperdom")])

# Calculate threshold of hyperdominant for plot
HD.min <- Compare.RAD %>%
  filter(hyperdom == "hyperdominant") %>%
  group_by(type) %>%
  summarise(MaxRank = max(RelRank), MinPS = min(PopSize),
            .groups = "keep")

# Plot rank abundance distributions
(RAD.plot <- ggplot(Compare.RAD, aes(x = RelRank, y = PopSize)) +
    scale_y_continuous(
      trans = "log10",
      breaks = scales::trans_breaks("log10", n = 4, function(x) 10 ^ (c(4:11))),
      labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
      limits = c(10 ^ 4, 10 ^ 10.5)) +
    geom_point(aes(colour = type)) +
    theme(legend.position = c(0.35,.2), axis.text = element_text(size = 12)) +
    scale_colour_manual(values = cols2) +
    scale_x_continuous(breaks = c(0, 1), labels = c("most common", "most rare"),
                       limits=c(0, 1.1)) +
    geom_segment(x = HD.min$MaxRank[1],
                 y = log10(HD.min$MinPS[1] * 0.5),
                 xend = HD.min$MaxRank[1],
                 yend = log10(HD.min$MinPS[1] * 2), size = 1) +
    geom_segment(x = HD.min$MaxRank[2],
                 y = log10(HD.min$MinPS[2] * 0.5),
                 xend = HD.min$MaxRank[2],
                 yend = log10(HD.min$MinPS[2] * 2), size = 1) +
    labs(x = "Species rank", y = "Population size", colour = ""))

# Calculate quantiles of pop size distributions for corals
(quantiles.Coral <- data.frame(
  quantiles = quantile(PopSizes$PopSizeMean * 10^9, probs = c(0.1,0.5,0.9)),
  label = as.character(c("10%", "50%", "90%"))))

# Histogram of coral population sizes
(Histo.Coral <- ggplot(PopSizes, aes(x = PopSizeMean*10^9)) +
    geom_histogram(bins = 15, fill = "#0070C0") +
    geom_vline(xintercept = quantiles.Coral$quantiles, colour = "black",
               linetype = 2) +
    geom_text(data = quantiles.Coral, aes(x = quantiles, label = label, y = 70),
              colour = "black", angle = 90, vjust = -.3, hjust = -.5) +
    scale_x_continuous(
      trans="log10",
      breaks = scales::trans_breaks("log10", function(x) 10^(c(4:11))),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(10^4,10^11)) +
    ylim(0,95) +
    theme(axis.text = element_text(size = 12)) +
    labs(x="Population size", y="Number of species"))

# Calculate number of coral species per IUCN category for labels
SampleSize.IUCN <- PopSizes %>%
  group_by(IUCNstatus) %>%
  summarise(Freq = n()) %>%
  filter(!is.na(IUCNstatus))

# Boxplot of population sizes per IUCN category
(Box.IUCN <- PopSizes %>%
    filter(!is.na(IUCNstatus)) %>%
    ggplot(aes(x = reorder(IUCNstatus, desc(IUCNstatus)),
               y = PopSizeMean * 10^9, fill = IUCNstatus)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(
      trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^(c(4:11))),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(10^4, 10^11.5)) +
    xlab("") + ylab("Population size") +
    scale_fill_manual(name = "IUCN status", values = cols, na.value = "grey") +
    geom_text(data = SampleSize.IUCN,
              aes(x = reorder(IUCNstatus, desc(IUCNstatus)),
                  label = paste0("n=", Freq), y = 10^10.8), hjust = 0) +
    theme(axis.text = element_text(size = 12), legend.position = "") +
    coord_flip())

### FIGURE 2
Figure2 <- cowplot::plot_grid(RAD.plot, Histo.Coral, Box.IUCN, ncol = 1,
                              labels = c("a","b","c"), align = "v", label_size = 14)

ggsave("figures/Fig2_PopSizesOverview.pdf", width = 12, height = 20, units = "cm")


# Calculate reef area within bbx vs. within total range of species
ReefArea.Species <- data.frame()
reefs.rast.CEIP2 <- as(reefs.rast.CEIP, 'SpatialGridDataFrame')
reefs.rast.IP2 <- as(reefs.rast.IP, 'SpatialGridDataFrame')

for (i in 1:length(unique(PopSizes$Species))) {

  Species = unique(PopSizes$Species[i])

  reefs.rast.CEIP2$PresAb <- sapply(over(
    reefs.rast.CEIP2,
    geometry(Ranges[Ranges@data$BINOMIAL == Species,]),
    returnList = T), length)

  reefs.rast.IP2$PresAb <- sapply(
    over(reefs.rast.IP2, geometry(Ranges[Ranges@data$BINOMIAL == Species,]),
         returnList = T), length)

  ReefArea.Sp = data.frame(
    Species = Species,
    ReefArea.CEIP = sum(reefs.rast.CEIP2$layer[reefs.rast.CEIP2$PresAb == 1], na.rm=T),
    ReefArea.IP = sum(reefs.rast.IP2$layer[reefs.rast.IP2$PresAb == 1], na.rm=T))

  ReefArea.Species <- bind_rows(ReefArea.Sp, ReefArea.Species)
}

ReefArea.Species$PropReefArea <-
  ReefArea.Species$ReefArea.CEIP / ReefArea.Species$ReefArea.IP

PopSizes <- PopSizes %>%
  left_join(., ReefArea.Species[, c("Species","PropReefArea")], by = "Species")

ED_Figure_1 <- ggplot(PopSizes, aes(PropReefArea)) +
  geom_histogram(bins = 10) +
  labs(x = "proportion of geographic range", y = "number of species")

ggsave("figures/ED_Fig_1_PropGeoRange.pdf", width = 6, height = 6, units = "cm")

# Calculate reef area within study domain relative to reef area in Indo-Pacific and world
sum(reefs.rast.CEIP@data@values, na.rm = T) / sum(reefs.rast.IP@data@values, na.rm = T)
sum(reefs.rast.CEIP@data@values, na.rm = T) / sum(reefs.rast.global@data@values, na.rm = T)

# Save population size table
PopSizes %>%
  dplyr::select(Species, PopSizeLower, PopSizeMean, PopSizeUpper, IUCNstatus) %>%
  arrange(desc(PopSizeMean)) %>%
  write.csv(file = "figures/DataS2_PopulationSizes.csv")

# Sensitivity analysis: calculate population sizes in countries were sampled
HPLC.mean <- HPLC.posterior_rescale %>%
  filter(type == "rescaled") %>%
  group_by(HABITAT) %>%
  summarise(Mean.HPLC = mean(HPLC),
            .groups = "keep") %>%
  ungroup() %>%
  mutate(Depth = ifelse(HABITAT == "Crest", "crest",
                        ifelse(HABITAT == "Slope/Fore Reef", "slope",
                               ifelse(HABITAT == "Flat/Back Reef", "flat",
                                      HABITAT)))) %>%
  dplyr::select(-HABITAT)

MeanRelAbund.RegHab <- RelAbund %>%
  group_by(Species, Region, Depth) %>%
  summarise(MeanRelAbund = mean(RelAbund_True),
            .groups = "keep")

PopSizes.Sens <- MeanRelAbund.RegHab %>%
  left_join(., ReefArea.byCountry, by = c("Region"="COUNTRY")) %>%
  left_join(., HPLC.mean, by = "Depth") %>%
  mutate(PopSize.RegHab = MeanRelAbund * ReefArea * Mean.HPLC * 28.9 * (10^6)) %>%
  group_by(Species) %>%
  summarise(PopSizeSens = sum(PopSize.RegHab),
            .groups = "keep")

# reef area in five countries were species abundance data were collected divided by
# reef area in study domain
sum(ReefArea.byCountry$ReefArea[
  ReefArea.byCountry$COUNTRY %in% c("Indonesia","Papua New Guinea",
                                    "American Samoa","Solomon Islands",
                                    "French Polynesia")]) / 107700

RAD.Sens <- Coral.RAD %>%
  dplyr::select(Species, RelRank, PopSize) %>%
  left_join(PopSizes.Sens, by = "Species")

Compare.RAD <- RAD.Sens %>%
  ggplot() +
  scale_y_continuous(
    trans = "log10",
    breaks = scales::trans_breaks("log10", n = 4, function(x) 10 ^ (c(4:11))),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
    limits=c(10 ^ 4, 10 ^ 10.5)) +
  geom_point(aes(x = RelRank, y = PopSize), colour = "blue", alpha = .5) +
  geom_point(aes(x = RelRank, y = PopSizeSens), colour = "orange", alpha = .5) +
  labs(x = "Species rank", y = "Population size", colour = "") +
  scale_x_continuous(breaks = c(0, 1), labels = c("most common", "most rare"),
                     limits = c(0, 1.1)) +
  theme_classic()

Compare.box <- RAD.Sens %>%
  gather(key, value, 3:4) %>%
  ggplot(aes(x = key, y = value, fill = key)) +
  scale_y_continuous(
    trans="log10",
    breaks = scales::trans_breaks("log10", n = 4, function(x) 10 ^ (c(4:11))),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
    limits = c(10 ^ 4, 10 ^ 10.5)) +
  geom_boxplot() +
  scale_fill_manual(values = c("blue","orange")) +
  labs(x = "", y = "Population size") +
  theme_classic() +
  theme(legend.position = "none")

cowplot::plot_grid(Compare.RAD, Compare.box, ncol = 2, labels = "auto",
                   rel_widths = c(1,.7))

ggsave("figures/Fig.PopSize_Sensitivity.pdf", width = 6, height = 2.5)

### Saving output files
save(PopSizes, PopSizes.Region, NC_Species.bootstrap2,
     NC_Species.bootstrap, file = "RData/OUT6_PopSizes.RData")
