#### Composition of Indo-Pacific coral fauna by family and growth form
#### Andreas Dietzel, August 2018

rm(list = ls())

library(tidyverse)

load("RData/Input_Files.RData")
load("RData/OUT6_PopSizes.RData")

#### ABUNDANCES BY FAMILY AND GROWTH FORM
SupportData$Morpho <-
  ifelse(grepl('tables', SupportData$GrowthForm_CTD, ignore.case = T), "tabular",
         ifelse(grepl('massive', SupportData$GrowthForm_CTD, ignore.case = T), "massive",
                ifelse(grepl('encrusting', SupportData$GrowthForm_CTD, ignore.case = T), "encrusting",
                       ifelse(grepl('branching', SupportData$GrowthForm_CTD, ignore.case = T),
                              "branching",
                              ifelse(grepl('corymbose|digitate|hispidose', SupportData$GrowthForm_CTD,
                                           ignore.case = T),
                                     "corymbose|digitate|hispidose",
                                     as.character(SupportData$GrowthForm_CTD))))))

# Fill in gaps of growth form data not available on Coral Traits Database

SupportData$Morpho[SupportData$Species == "Acropora muricata"] <- "branching"
SupportData$Morpho[SupportData$Species == "Goniastrea aspera"] <- "massive"
SupportData$Morpho[SupportData$Species == "Goniastrea australensis"] <- "massive"
SupportData$Morpho[SupportData$Species == "Favia stelligera"] <- "columnar"
SupportData$Morpho[SupportData$Species == "Favia pallida"] <- "massive"
SupportData$Morpho[SupportData$Species == "Favia danae"] <- "massive"
SupportData$Morpho[SupportData$Species == "Favia lizardensis"] <- "massive"
SupportData$Morpho[SupportData$Species == "Favia maritima"] <- "massive"
SupportData$Morpho[SupportData$Species == "Favia maxima"] <- "massive"
SupportData$Morpho[SupportData$Species == "Favia rotundata"] <- "massive"
SupportData$Morpho[SupportData$Species == "Favia matthaii"] <- "massive"
SupportData$Morpho[SupportData$Species == "Favia favus"] <- "massive"
SupportData$Morpho[SupportData$Species == "Caulastraea curvata"] <- "plocoid"
SupportData$Morpho[SupportData$Species == "Caulastraea tumida"] <- "plocoid"
SupportData$Morpho[SupportData$Species == "Favites russelli"] <- "massive"
SupportData$Morpho[SupportData$Species == "Montastrea annuligera"] <- "encrusting"
SupportData$Morpho[SupportData$Species == "Montastrea magnistellata"] <- "massive"
SupportData$Morpho[SupportData$Species == "Montastrea salebrosa"] <- "massive"
SupportData$Morpho[SupportData$Species == "Fungia danai"] <- "solitary"
SupportData$Morpho[SupportData$Species == "Herpolitha weberi"] <- "solitary"
SupportData$Morpho[SupportData$Species == "Fungia horrida"] <- "solitary"
SupportData$Morpho[SupportData$Species == "Paraclavarina triangularis"] <- "branching"

# Calculate sum of population sizes by growth form and family

Total.Family <- PopSizes %>%
  group_by(family_morphology) %>%
  summarise(Total = sum(PopSizeMean),
            Freq = n(),
            .groups = "keep") %>%
  mutate(Type = "morphological", family = family_morphology)

Total.FamilyMol <- PopSizes %>%
  group_by(family_molecules) %>%
  summarise(Total = sum(PopSizeMean),
            Freq = n(),
            .groups = "keep") %>%
  mutate(Type = "molecular",
         family = family_molecules)

Total.Morpho <- PopSizes %>%
  group_by(GrowthForm_CTD) %>%
  summarise(Total = sum(PopSizeMean),
            Freq = n(),
            .groups = "keep") %>%
  filter(!is.na(GrowthForm_CTD))

Total.MorphoBin <- PopSizes %>%
  left_join(.,SupportData[,c("Species","Morpho")], by = "Species") %>%
  group_by(Morpho) %>%
  summarise(Total = sum(PopSizeMean),
            Freq = n(),
            .groups = "keep") %>%
  filter(!is.na(Morpho)) %>%
  mutate(Morpho = gsub("_"," ",Morpho))

Total.FamilyCombined <- bind_rows(Total.Family[,c("Type","Freq","family","Total")],
                                  Total.FamilyMol[,c("Type","Freq","family","Total")]) %>%
  complete(Type, family, fill = list(Freq = 0, Total = 0))

Total.FamilyBin <- Total.Family %>%
  mutate(Fam.bin =
           ifelse(family %in%
                    c("Acroporidae","Poritidae","Pocilloporidae",
                      "Faviidae","Merulinidae","Agariciidae"),
                  as.character(family), "Others")) %>%
  group_by(Fam.bin,Type) %>%
  summarise(Total = sum(Total),
            Freq = sum(Freq),
            .groups = "keep")

(Total.FamilyBin.plot <- ggplot(Total.FamilyBin, aes(x=reorder(Fam.bin,Total),y=Total)) +
    geom_bar(stat="identity") + coord_flip() +
    labs(y = "Billion colonies", x = "") +
    geom_text(aes(Fam.bin, label = paste0(Freq), y = 125),
              position = position_dodge(width = 1), hjust = 0) +
    ylim(0, 140) + theme(legend.position = "top"))

Total.Morpho$GrowthForm_CTD <- gsub("_"," ", Total.Morpho$GrowthForm_CTD)

(Total.Morpho.plot <- ggplot(Total.Morpho,
                             aes(x = reorder(GrowthForm_CTD, Total), y = Total)) +
    geom_bar(stat = "identity") + coord_flip() +
    labs(y = "Billion colonies", x = "") +
    geom_text(aes(x = GrowthForm_CTD, label = paste0(Freq), y = 70),
              colour = "black", hjust = 0) +
    ylim(0, 80))

cowplot::plot_grid(Total.FamilyBin.plot, Total.Morpho.plot,
                   ncol = 1, labels = c("a","b"),
                   rel_heights = c(.6,.8), align = "v")

ggsave("figures/FigS6_GrowthForm_Family.pdf", width = 4.5, height = 5)


