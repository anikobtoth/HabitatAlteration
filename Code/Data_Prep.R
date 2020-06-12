## Code to prepare raw data in ./Raw_Data folder 
## Outputs will be placed in ./Data folder and used in Analysis_Script.R
library(reshape2)
library(tidyverse)
#source('./Code/HelperFunctions.R')

l <- list.files("../Raw_Data", ".txt", full.names = T)

dat <- lapply(l, read.delim) %>% setNames(list.files("../Raw_Data", ".txt"))


#### Species metadata ####
spp <- dat[grep("species", names(dat))] %>% bind_rows() %>% 
  select(species, life.form, individuals, samples, diet.1, diet.2, mass_g = mass..g.)

rownames(spp) <-  gsub(" ", "_", spp$species)

# redo guild categories ####
spp$diet.2[spp$diet.2 == ""] <- spp$diet.1[spp$diet.2 == ""]
spp$guild[spp$diet.1 == "invertivore" & spp$diet.2 == "invertivore"] <- "I"
spp$guild[spp$diet.1 == "invertivore" & spp$diet.2 == "frugivore"] <- "FI"
spp$guild[spp$diet.1 == "frugivore" & spp$diet.2 == "invertivore"] <- "FI"
spp$guild[spp$diet.1 == "nectarivore" & spp$diet.2 == "frugivore"] <- "FN"
spp$guild[spp$diet.1 == "frugivore" & spp$diet.2 == "nectarivore"] <- "FN"
spp$guild[spp$diet.1 == "nectarivore" & spp$diet.2 == "nectarivore"] <- "N"
spp$guild[spp$diet.1 == "frugivore" & spp$diet.2 == "frugivore"] <- "F"
spp$guild[spp$diet.1 == "invertivore" & spp$diet.2 == "carnivore"] <- "CI"
spp$guild[spp$diet.1 == "carnivore" & spp$diet.2 == "invertivore"] <- "CI"
spp$guild[spp$diet.1 == "carnivore" & spp$diet.2 == "carnivore"] <- "C"
spp$guild[spp$diet.1 == "invertivore" & spp$diet.2 == "piscivore"] <- "CI"
spp$guild[spp$diet.1 == "piscivore" & spp$diet.2 == "piscivore"] <- "C"
spp$guild[spp$diet.1 == "nectarivore" & spp$diet.2 == "invertivore"] <- "IN"
spp$guild[spp$diet.1 == "invertivore" & spp$diet.2 == "nectarivore"] <- "IN"
spp$guild[spp$diet.1 == "sanguinivore" & spp$diet.2 == "sanguinivore"] <- "S"
spp$guild[spp$diet.1 == "frugivore" & spp$diet.2 == "granivore"] <- "FG"
spp$guild[spp$diet.1 == "granivore" & spp$diet.2 == "frugivore"] <- "FG"
spp$guild[spp$diet.1 == "granivore" & spp$diet.2 == "granivore"] <- "G"
spp$guild[spp$diet.1 == "invertivore" & spp$diet.2 == "granivore"] <- "IG"
spp$guild[spp$diet.1 == "granivore" & spp$diet.2 == "invertivore"] <- "IG"
spp$guild[spp$diet.1 == "nectarivore" & spp$diet.2 == "granivore"] <- "NG"



##### Site metadata ####
sitedat <- dat[grep("samples", names(dat))] %>% setNames(c("bat", "bird")) %>% bind_rows(.id = "taxon") %>% 
  select(taxon, sample.no, sample.name, country, ecozone, latitude, longitude, habitat, altered.habitat, MAT, MAP, richness, fragment.size)
sitedat$status[!sitedat$altered.habitat == ""] <- "Altered" 
sitedat$status[sitedat$altered.habitat == "" | sitedat$fragment.size == "1000 - 10000 ha"] <- "Unaltered" 

sitedat$siteid <- paste0("X", sitedat$sample.no)

#### Occurrences ####
PAn <- dat[grep("register", names(dat))] %>% 
  purrr::map(~filter(.,!species %in% c("sp.", "spp.", "sp. l", "indet")) %>% 
              mutate(name = paste(genus, species, sep = "_")) %>% 
              dcast(name~sample.no, value.var = "count", fill = 0) %>% 
              namerows()) %>%
  setNames(c("bat", "bird"))


#### Calculate habitat preference ####
unalt_sites <- sitedat$sample.no[sitedat$status == "Unaltered"]

nsamp <- table(sitedat$taxon, sitedat$status)

u = purrr::map(PAn, ~apply(., 1, function(x) names(which(x > 0))) %>% 
             sapply(function(x) length(which(x %in% (unalt_sites))))) 
a = purrr::map(PAn, ~apply(., 1, function(x) names(which(x > 0))) %>% 
             sapply(function(x) length(which(!x %in% (unalt_sites)))))
cosmo <- map2(u, a, data.frame) %>% 
  map(~mutate(.,name = rownames(.))) %>% bind_rows(.id = 'taxon') %>% 
  setNames(c("taxon", "unaltered", "altered", "name")) %>% 
  mutate(
    occurrences = unaltered+altered, 
    Talt = nsamp[taxon,"Altered"],
    Tunalt = nsamp[taxon, "Unaltered"],
    cosmopolitan = FETmP_(Talt, Tunalt, altered, unaltered)
  )

cosmo$group[qnorm(cosmo$cosmopolitan) < -1] <- "restricted"
cosmo$group[qnorm(cosmo$cosmopolitan) > 1] <- "synanthropic"
cosmo$group[is.na(cosmo$group)] <- "cosmopolitan"

cosmo$group_abbr <- substr(cosmo$group, start = 1, stop = 5)
rownames(cosmo) <- cosmo$name

