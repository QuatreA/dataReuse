
#STAVIRO and BRUVs dataset

# Written by jacquomo.monk@utas.edu.au

## Install required packages
# install.packages(c("googlesheets4"))

# Load libraries
library(janitor)
library(tidyverse)
library(iNEXT)
library(cowplot)

# unique(sppR$scientific)
setwd("G:\\.shortcut-targets-by-id\\17M0oqtCKud0SoZZnJLNFGXemSxFsTlw0\\QuatreA\\DATA EXCHANGE\\AMBIO STAVIRO DATA")
# Read the stavio data from the CSV file
stavio.fish.dat <- read.csv("UnitobsEspeceMetriques.csv")%>%
  filter(!str_detect(scient.name, "sp.")) %>%
  glimpse()

setwd("G:\\.shortcut-targets-by-id\\17M0oqtCKud0SoZZnJLNFGXemSxFsTlw0\\QuatreA\\DATA EXCHANGE\\IMAS_EMR_BRUV")
# Read the stavio data from the CSV file
bruv.fish.dat <- read.csv("EMR_2020_maxn.summary_2022-02-07.csv")%>%
  filter(!str_detect(scientific, "sp.")) %>%
  filter(!str_detect(scientific, "sp")) %>%
  filter(!str_detect(scientific, "Unknown")) %>%
  filter(!str_detect(scientific, "spp")) %>%
  glimpse()
      
# Calculate species richness
stavio.sppR <- stavio.fish.dat %>%
  mutate(ID = paste(observation.unit, sep = " ")) %>%
  group_by(ID) %>%
  summarise(counts = sum(pres.abs > 0, na.rm = TRUE))%>%
  glimpse()

bruv.sppR <- bruv.fish.dat %>%
  mutate(ID = paste(sample, sep = " ")) %>%
  mutate(pres.abs= ifelse(maxn > 0, 1, 0),)%>%
  group_by(ID) %>%
  summarise(counts = sum(pres.abs > 0, na.rm = TRUE))%>%
  glimpse()
# names(stavio.fish.dat)

# Make data wide
#Stavio
stavio.fish.w.dat<-stavio.fish.dat%>%
  # mutate(ID = paste(observation.unit, sep = " ")) %>%
  dplyr::select(!c("species.code","number","number.sd","mean.length","weight","mean.weight","density","density.max","density.sd",
                   "biomass","biomass.max","biomass.sd","pres.abs","year","site","protection.status","biotop","latitude","longitude",
                   "habitat1","habitat2","habitat3" ,"family", "genus","species"))%>%
  group_by(observation.unit) %>%
  pivot_wider(names_from = scient.name, values_from = number.max)%>%
  column_to_rownames(var = "observation.unit") %>%
  janitor::clean_names()%>%
  glimpse()

#Bruvs
names(bruv.fish.dat)
bruv.fish.w.dat<-bruv.fish.dat%>%
  dplyr::select(!c("family","genus","species","marine.region", "latitude","longitude","date","time","site",
                   "location","status","depth","successful.count","successful.length","NetName","ResName",
                   "ResArea","zone","ZoneIUCN","NatLegend","Area_KM2","GAZ_DATE","LATEST_GAZ","campaignid",
                   "project","id"))%>%
  group_by(sample) %>%
  pivot_wider(names_from = scientific, values_from = maxn)%>%
  column_to_rownames(var = "sample") %>%
  janitor::clean_names()%>%
  glimpse()


## convert into P/A data ##
stavio.fish.w.dat.pa <- decostand(stavio.fish.w.dat, method = "pa")
bruv.fish.w.dat.pa <- decostand(bruv.fish.w.dat, method = "pa")


## transpose the matrix to obtain a species by site matrix ##
stavio.fish.w.dat.pa.t <- t(stavio.fish.w.dat.pa)
bruv.fish.w.dat.pa.t <- t(bruv.fish.w.dat.pa)
  


## calculate the frequency (no. of sampling units in which taxa "X" is found => its frequency of occurrence##
## the first row of the matrix will correspond to the number of samples in the considered dive ##
stavio.fish_incfreq <- as.matrix(as.incfreq(stavio.fish.w.dat.pa.t))
bruv.fish_incfreq <- as.matrix(as.incfreq(bruv.fish.w.dat.pa.t))


## create the list of dives (or transects) for the iNEXT calculation ##
list<- list(stavio.fish_incfreq, bruv.fish_incfreq)
names(list) <- c("Stavio", "BRUV")

## create an object containing the list of colours per dive to be implemented in the ggplot  ##
color_dive <- c('#1b9e77','#d95f02')


## calculation of the rarefaction curve ##
#Stavro = 109 sampling units (drops), BRUV = 170 sampling units (drops)
inext_pa <-iNEXT(list, q = 0, datatype = "incidence_freq", endpoint = 170*2, nboot = 999) #endpoint = 340 sampling units corresponding to 
#twice the size of the referecnce sample size (170 images); hill number q = 0 (richness), q = 1 (Shannon), q = 2 (Simpson)

## representation of the curve with ggplot ##
gg_raref <- ggiNEXT(inext_pa, type = 1, se = T) + 
  scale_color_manual(values = color_dive,
                     labels = c("Stavio", "BRUV")) +
  scale_fill_manual(values = color_dive,
                    labels = c("Stavio", "BRUV"))  +
  scale_shape_manual(values = c(19, 19))+
  #theme(legend.position = "right") +
  labs(y = "Taxonomic richness", x = "Number of sampling units") +
  theme_bw() +
  guides(shape = "none", 
         color = guide_legend("Fish"),
         fill = guide_legend("Fish")) +
  theme(panel.grid = element_blank(),  #get rid of grey backgound
        axis.title.x = element_text(size=18), # change the size of the x-axis title
        axis.text.x = element_text(size = 16), # change the size of the x-labels
        axis.title.y = element_text(size = 18), # change the size of the ordinate title
        axis.text.y = element_text(size = 16),
        panel.border = element_blank(),  # Remove the panel border
        axis.line.x = element_line(),   # Add back the x-axis line
        axis.line.y = element_line(),    # Add back the y-axis line
        legend.position = "none",
        # legend.text = element_text(size = 12), # change the size of the caption labels
        # legend.title = element_text(size = 14) # change the size of the caption titles
        )
gg_raref


## Save the plot ##
setwd("G:\\.shortcut-targets-by-id\\17M0oqtCKud0SoZZnJLNFGXemSxFsTlw0\\QuatreA\\EBV article")

ggsave(filename = "raref_sample-based-FishpaMATRIX.png", plot = gg_raref, width = 9, height = 10, dpi = 300)





#Calculate total number of species
stavio.sppR.count <- stavio.fish.dat %>%
  mutate(ID = paste(scient.name, sep = " ")) %>%
   group_by(ID) %>%
  summarise(counts = sum(pres.abs > 0, na.rm = TRUE))%>%
  glimpse()


stavio.total.filter <- stavio.fish.dat %>%
  glimpse()

sum(stavio.total.filter$number.max)# total abundance



#Calculate total number of species
bruv.sppR.count <- bruv.fish.dat %>%
  mutate(ID = paste(scientific, sep = " ")) %>%
  mutate(pres.abs= ifelse(maxn > 0, 1, 0),)%>%
  group_by(ID) %>%
  summarise(counts = sum(pres.abs > 0, na.rm = TRUE))%>%
  glimpse()

bruv.total.filter <- bruv.fish.dat %>%
  glimpse()

sum(bruv.total.filter$maxn)# total abundance

#Calculate margalef diversity

# Calculate the total number of individuals (N)
total_individuals <- sum(stavio.fish.dat$number.max)

# Calculate the number of unique species (S)
unique_species <- length(unique(stavio.fish.dat$scient.name))
unique_species
# Calculate the Margalef diversity index
margalef_index <- (unique_species - 1) / log(total_individuals)

# margalef_index <-log(unique_species) / log(total_individuals)


# Print the Margalef diversity index
print(margalef_index)


# Calculate the total number of individuals (N)
total_individuals <- sum(bruv.fish.dat$maxn)

# Calculate the number of unique species (S)
unique_species <- length(unique(bruv.fish.dat$scientific))
unique_species
# Calculate the Margalef diversity index
margalef_index <- (unique_species - 1) / log(total_individuals)

# margalef_index <-log(unique_species) / log(total_individuals)


# Print the Margalef diversity index
print(margalef_index)


# Calculate Shannon diversity index using vegan
shannon_diversity <- diversity(stavio.fish.dat$number.max)
print(shannon_diversity)
simpson_diversity <- diversity(stavio.fish.dat$number.max, index = "simpson")
print(simpson_diversity)


shannon_diversity <- diversity(bruv.fish.dat$maxn)
print(shannon_diversity)
simpson_diversity <- diversity(bruv.fish.dat$maxn, index = "simpson")
print(simpson_diversity)




# Calculate the relative abundance
relative_abundance <- bruv.fish.dat$maxn / total_individuals

# Calculate the Simpson diversity index
simpson_index <- sum(relative_abundance^2)


#Number of samples
unique_count <- length(unique(as.vector(as.matrix(bruv.fish.dat$sample))))
unique_count <- length(unique(as.vector(as.matrix(stavio.fish.dat$observation.unit))))
unique_count


#Calculate abundance of fish > 20cm BRUV
bruv.a20<- read.csv("EMR_2020_length.summary_2022-02-07.csv")%>%
  tidyr::uncount(number)%>%
  filter(!str_detect(scientific, "sp.")) %>%
  filter(!str_detect(scientific, "sp")) %>%
  filter(!str_detect(scientific, "Unknown")) %>%
  filter(!str_detect(scientific, "spp")) %>%
  dplyr::group_by(sample,scientific)%>%
  dplyr::filter(length>=200)%>%
  dplyr::summarise(value = n())%>%
  dplyr::group_by(sample)%>%
  dplyr::ungroup()%>%
  dplyr::summarise(mean.a20 = mean(value))%>%
  glimpse()

#Calculate biomass of fish > 20cm BRUV
bruv.b20<- read.csv("EMR_2020_mass.summary_2022-02-07.csv")%>%
  tidyr::uncount(number)%>%
  filter(!str_detect(scientific, "sp.")) %>%
  filter(!str_detect(scientific, "sp")) %>%
  filter(!str_detect(scientific, "Unknown")) %>%
  filter(!str_detect(scientific, "spp")) %>%
  dplyr::group_by(sample,scientific)%>%
  dplyr::filter(length>=200)%>%
  dplyr::summarise(value = sum(mass.kg*1000))%>%
  dplyr::group_by(sample)%>%
  dplyr::ungroup()%>%
  dplyr::summarise(mean.b20 = mean(value))%>%
  glimpse()


#Calculate abundance of fish > 20cm stavio
# unique(sppR$scientific)
setwd("G:\\.shortcut-targets-by-id\\17M0oqtCKud0SoZZnJLNFGXemSxFsTlw0\\QuatreA\\DATA EXCHANGE\\AMBIO STAVIRO DATA")
# Read the stavio data from the CSV file
stavio.fish.dat.a20 <- read.csv("UnitobsEspeceClassetailleMetriques.csv")%>%
  filter(!str_detect(scient.name, "sp.")) %>%
  tidyr::uncount(number.max)%>%
  dplyr::group_by(observation.unit,scient.name)%>%
  dplyr::filter(size.class %in%c("G"))%>%
  dplyr::summarise(value = n())%>%
  dplyr::group_by(observation.unit)%>%
  dplyr::ungroup()%>%
  dplyr::summarise(mean.a20 = mean(value))%>%
  glimpse()



#Calculate biomass of fish > 20cm stavio
# unique(sppR$scientific)
setwd("G:\\.shortcut-targets-by-id\\17M0oqtCKud0SoZZnJLNFGXemSxFsTlw0\\QuatreA\\DATA EXCHANGE\\AMBIO STAVIRO DATA")
# Read the stavio data from the CSV file
stavio.fish.dat.a20 <- read.csv("UnitobsEspeceClassetailleMetriques.csv")%>%
  filter(!str_detect(scient.name, "sp.")) %>%
  tidyr::uncount(number.max)%>%
  dplyr::group_by(observation.unit,scient.name)%>%
  dplyr::filter(size.class %in%c("G"))%>%
  dplyr::summarise(value = sum(mass.kg*1000))%>%
  dplyr::group_by(sample)%>%
  dplyr::ungroup()%>%
  dplyr::summarise(mean.b20 = mean(value))%>%
  glimpse()



