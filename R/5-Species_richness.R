# Accounting for gradients in species richness across study domain
# Code author: Andreas Dietzel
# Last modified: February 2, 2020

rm(list = ls())

# Load libraries and files
library(tidyverse)
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(fields)

load("RData/Input_Files.RData")
load("RData/OUT2_Number_Corals.RData")

# Diversity Map of the Indo-Pacific
x.range <- as.integer(c(-20000000, 9500000))
y.range <- as.integer(c(-4500000, 4500000))
DiversityMap <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = res),
                            y = seq(from = y.range[1], to = y.range[2], by = res))
coordinates(DiversityMap) <- ~ x + y # assign coordinates
proj4string(DiversityMap) <- proj_CEA # project
DiversityMap$Species <- sapply(over(DiversityMap, geometry(Ranges),
                                    returnList = T), length)
gridded(DiversityMap) <- TRUE

# Finding regional centroids
Coords.Sites <- Locations@coords %>% as.data.frame() %>%
  mutate(RegionCode = rep(1:5, each = 3)) %>%
  group_by(RegionCode) %>%
  summarise_all(mean) %>%
  mutate(Region = unique(Locations@data$Region)) %>%
  as.data.frame()

RegCentroids <- SpatialPointsDataFrame(
  coords = Coords.Sites[, c(2, 3)], data = Coords.Sites,
  proj4string = CRS(proj4string(Locations)))

# Calculating distance matrix
coords.grd <- coordinates(reefs.rast.IP)
coords.Sites <- coordinates(RegCentroids)

distance.matrix <- as.data.frame(
  rdist.earth(coords.grd, coords.Sites, miles = F, R = NULL))

colnames(distance.matrix) <- unique(Locations@data$Region)
distance.matrix$CellNo <- seq(1 : length(distance.matrix$Indonesia))

# Presence and absence of species at islands according to range
PresAb.Range <- data.frame(Species = as.character(),
                           Region = as.double())

for (i in 1:length(Ranges@data$BINOMIAL)) {
  Species <- Ranges@data$BINOMIAL[i]

  PresAb.Range.Sp <- as.data.frame(
    gWithin(RegCentroids, Ranges[Ranges@data$BINOMIAL == Species, ], byid = TRUE))

  PresAb.Range.Sp$Species <- Species

  PresAb.Range <- rbind(PresAb.Range, PresAb.Range.Sp)
}

PresAb.Region <- PresAb.Range %>%
  tidyr::gather(RegionCode, presab, -Species) %>%
  arrange(desc(Species)) %>%
  mutate(RegionCode = rep.int(1:5, length(PresAb.Range$Species)),
         PresAb.Range = as.integer(presab)) %>%
  dplyr::select(-presab)

# Presence and absence of species at islands according to abundance dataset
PresAb <- Interc %>%
  mutate(RegionCode = as.numeric(substr(IslandCode, 0, 1))) %>%
  mutate(RegionCode = ifelse(RegionCode > 3, RegionCode - 1, RegionCode)) %>%
  group_by(Species, RegionCode) %>%
  dplyr::summarise(nr.int = n(), .groups = "keep") %>%
  ungroup() %>%
  tidyr::complete(Species, RegionCode, fill = list(nr.int = 0)) %>%
  left_join(., PresAb.Region,
            by=c("Species"="Species", "RegionCode"="RegionCode")) %>%
  mutate(PresAb.Either = ifelse(PresAb.Range > 0 | nr.int > 0, 1, 0)) %>%
  mutate(PresAb.Data = ifelse(nr.int > 0, 1, 0)) %>%
  dplyr::select(-nr.int)

# Species richness in cells with reef area > 0 (to reduce computation)
DiversityMap.df <- data.frame(
  NrSpecies = DiversityMap@data$Species,
  CellNo = 1:length(DiversityMap@data$Species))

# Join to df with confidence intervals of number of corals per grid cell
RasterCells.NC <- RasterCells.NC %>%
  filter(Mean > 0) %>%
  left_join(., DiversityMap.df, by = "CellNo")

distmat.long <- distance.matrix %>%
  tidyr::gather(., Region, distance, c(1:5)) %>%
  filter(CellNo %in% RasterCells.NC$CellNo) %>%
  left_join(., RegCentroids@data[,c(1,4)], by = c("Region" = "Region"))

# Number of species at regional centroids
SpeciesRichness.Regions <- data.frame(
  Region = c("Indonesia", "Papua New Guinea", "Solomon Islands",
             "American Samoa", "French Polynesia"),
  RegionCode = 1:5,
  NrSpecies.Region = sapply(over(RegCentroids, geometry(Ranges),
                                 returnList = T), length))

# Calculate for each species the total number of corals within its range, adjusted by changes in species richness
NrCorals.Range <- lapply(

  seq_along(unique(Interc$Species)),

  function(i) {

    Species <- unique(Interc$Species)[i]

    Cells.Range <- over(
      grd.IP,
      geometry(Ranges[as.character(Ranges$BINOMIAL) == Species,]))

    Cells.Range2 <- data.frame(InRange = Cells.Range,
                               CellNo = 1:length(Cells.Range))

    Cells.Range3 <- filter(Cells.Range2, InRange > 0)

    NrCorals.Ran <- distmat.long %>%
    left_join(., PresAb[PresAb$Species == Species,], by = "RegionCode") %>%
    filter(PresAb.Either > 0) %>%
    group_by(CellNo) %>%
    arrange(distance, .by_group = T) %>%
    top_n(1, distance) %>%
    dplyr::select(CellNo, RegionCode, Species) %>%
    filter(CellNo %in% Cells.Range3$CellNo) %>%
    left_join(., SpeciesRichness.Regions, by = "RegionCode") %>%
    left_join(., RasterCells.NC, by = "CellNo") %>%
    group_by(Species, Region) %>%
    summarise(Lower.NC.SR = round(sum(Lower * NrSpecies.Region / NrSpecies)),
              Mean.NC.SR = round(sum(Mean * NrSpecies.Region / NrSpecies)),
              Upper.NC.SR = round(sum(Upper * NrSpecies.Region / NrSpecies)),
              Lower.NC = round(sum(Lower)),
              Mean.NC = round(sum(Mean)),
              Upper.NC = round(sum(Upper)),
              NrCells = n(),
              .groups = "keep")
  })

NrCorals.Range <- bind_rows(NrCorals.Range)

# Saving output files
save(NrCorals.Range, PresAb,
     SpeciesRichness.Regions,
     file = "RData/OUT5_NrCoralsRange.RData")
