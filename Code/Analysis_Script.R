##    Habitat alteration reduces food competition and local    ###
##      functional diversity in Neotropical bats and birds.    ### 
##                     Anikó B. Tóth                           ##

# Submitted 1 November 2020 
# Analysis script 

library(tidyverse)
library(vegan)
library(reshape2)
library(stringi)
library(parallel)

## Load helper functions
source('./Code/HelperFunctions.R')

## Prep raw data
source('./Code/Data_Prep.R')

###### Richness - run using raw abundance data ####
# can only be run on raw data because it requires a singleton count of abundances.
# cJ1
cJ1 <- map(PAn, map_dbl, cJ1rich) %>% map(~split(., names(.) %in% unalt_sites)) %>% 
  map(map, cbind) %>% map(map, data.frame) %>% map(bind_rows, .id = "status") %>% 
  bind_rows(.id = "taxon") %>% setNames(c("taxon", "status", "richness"))
cJ1$status <- plyr::revalue(cJ1$status, c("FALSE" = "Altered", "TRUE" = "Unaltered"))
# summary stats
cJ1 %>% group_by(taxon, status) %>% summarise(mean.rich = mean(richness), median.rich = median(richness))

# test for significant difference between altered and unaltered
cJ1 %>% group_by(taxon) %>% summarise(w=wilcox.test(richness~status, paired=FALSE)$p.value)

# CHAO 1
chao <- map(PAn, map_dbl, chao1) %>% map(~split(., names(.) %in% unalt_sites)) %>% 
  map(map, cbind) %>% map(map, data.frame) %>% map(bind_rows, .id = "status") %>% 
  bind_rows(.id = "taxon") %>% setNames(c("taxon", "status", "richness"))
chao$status <- plyr::revalue(chao$status, c("FALSE" = "Altered", "TRUE" = "Unaltered"))
#summary status
chao %>% group_by(taxon, status) %>% summarise(mean.rich = mean(richness, na.rm = T), median.rich = median(richness, na.rm = T))
#significance test
chao %>% group_by(taxon) %>% summarise(w=wilcox.test(richness~status, paired=FALSE)$p.value)


#### Match Biogeography ####
PAnu <- PAn %>% map(~return(.[,colnames(.) %in% unalt_sites])) # separate unaltered
PAna <- PAn %>% map(~return(.[,!colnames(.) %in% unalt_sites])) # and altered sites

coords <- sitedat %>% select(taxon, p.sample, latitude, longitude) %>% distinct() %>% split(.$taxon)

# Remove unaltered or altered sites that are not near a site of the other type. This is done to ensure the two sets have similar biogoegraphical distributions.
coords1 <- purrr::map2(coords, PAnu, function(x, y) x[x$p.sample %in% colnames(y),])
coords2 <- purrr::map2(coords, PAna, function(x, y) x[x$p.sample %in% colnames(y),])
keep <- map2(coords1, coords2, matchbiogeo) %>% map(unlist)

PAn <- map2(PAn, keep, function(x, y) return(x[,as.character(y)]))
PAn <- map(PAn, clean.empty, minrow = 1) # remove any species that now have no occurrences

# recalculate alt/unalt split from new PAn
PAnu <- PAn %>% map(~return(.[,colnames(.) %in% unalt_sites]))
PAna <- PAn %>% map(~return(.[,!colnames(.) %in% unalt_sites]))

# coords of sites we kept, for plotting later.
coords1.keep <- purrr::map2(coords, PAnu, function(x, y) x[x$p.sample %in% colnames(y),])
coords2.keep <- purrr::map2(coords, PAna, function(x, y) x[x$p.sample %in% colnames(y),])


#### Beta diversity and composition analyses #####
PAnb<- tobinary(PAn)

PA <- PAn  ## can run analyses with binary or abudance data

# with bray-curtis index
dist <- map(PA, ~t(.)) %>% map(vegdist, method = "bray")
# with Ochiai index
dist <- map(PA, ochiaiMatrix) %>% map(as.dist, upper = F) %>% map(~return(1-.))

# beta diversity using betadisper
beta <- map2(dist, PAn, function(x, y) betadisper(x, group = colnames(y) %in% unalt_sites) %>% anova)

# compositional change using adonis
data <- map(PA, ~t(.x) %>% data.frame() %>% rownames_to_column("p.sample") %>% mutate(ID = colnames(.x) %in% unalt_sites))
comp1 <- map2(dist, data, function(x, y) adonis(x~ID, data = y, permutations = 10000))

# compositional change using canonical correspondence analysis
dat <- map(data, select, contains("_"))
comp2 <- map2(dat, data, function(x, y) cca(x~ID, data = y) %>% anova.cca)

##### occupancy calculations (used later) ######
occ.count <- map(PAn, t) %>% map(data.frame) %>% map(~split(., rownames(.) %in% unalt_sites)) %>% 
  map(map, ~colSums(.)) %>% map(map, data.frame) %>% map(~merge(.[[1]], .[[2]], by = 0)) %>% map(setNames, c("name", "altered", "unaltered")) %>% bind_rows(.id = "taxon")
occ.count$name <- stri_replace_all_regex(occ.count$name, "[.]", " ")

uniquesp <- occ.count$name[which(occ.count$altered == 0 | occ.count$unaltered == 0)]
sharedsp <- occ.count$name[which(!occ.count$name %in% uniquesp)]

occ <- map(PAn, t)  %>% map(data.frame) %>% map(~split(., rownames(.) %in% unalt_sites)) %>% 
  map(map, ~colSums(.)/nrow(.)) %>% map(map, data.frame) %>% map(~merge(.[[1]], .[[2]], by = 0)) %>% map(setNames, c("name", "altered", "unaltered")) %>% bind_rows(.id = "taxon")

occ %>% group_by(taxon) %>% summarise(
  declined = length(which(altered < unaltered)), 
  increased = length(which(unaltered < altered)), 
  extirpated = length(which(altered == 0 & unaltered > 0)), 
  appeared = length(which(altered > 0 & unaltered == 0)), 
  shared = length(which(altered > 0 & unaltered > 0)), 
  total = length(altered))

#### Interaction Co-occurrence analysis #####

# remove singletons
#PAn.ns <- map(PAn, clean.empty, minrow = 2)

# Format data (full)
tables <- PAnb %>% map(~t(.)) %>% map(as.data.frame) %>% map(~split(., f = rownames(.) %in% unalt_sites)) %>% 
  map(map, ~t(.)) %>% map(map, clean.empty) %>% purrr::map(setNames, c("altered", "unaltered"))
# shared (common-only)
 shared_bats <- rownames(tables[[1]][[1]])[which(rownames(tables[[1]][[1]]) %in% rownames(tables[[1]][[2]]))]
 shared_birds <- rownames(tables[[2]][[1]])[which(rownames(tables[[2]][[1]]) %in% rownames(tables[[2]][[2]]))]
 tables <- map2(tables, list(shared_bats, shared_birds), function(x, y) map(x, function(z) return(z[y,])))

# Contingency table
contables <- map(tables, map, cont_table) %>% map(bind_rows, .id = "status") %>% 
  bind_rows(.id = "taxon") %>% diet_cat(spp, related = TRUE) %>% na.omit()

source('./Code/omega.R')
# Calculate omega 
## interaction, no bagging, no related, clipped to median
out_bat <- omega(contables, tax = "bat", type = "interaction") 
out_bird <- omega(contables, tax = "bird", type = "interaction")

##  interaction, no bagging, no related, all species
out_bat <- omega(contables, tax = "bat", type = "interaction", medn = FALSE)
out_bird <- omega(contables, tax = "bird", type = "interaction", medn = FALSE)


# Code for plotting omega function output
ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch)) + 
  geom_density() + facet_wrap(parameter~., scales = "free_x")

  
#### Analysis of co-occurrence at altered habitats ####
  # bats
  gld <- c("N", "I", "F", "CI")  
  tab <- tables[[1]][[1]]
  tab <- merge(tables[[1]][[1]], tables[[1]][[2]], by = 0) %>% namerows
  pr.bat <- percent.occupancy.by.guild(tab, gld, "bat", sitedat) 

  # birds
  gld <- c("N", "I", "F", "FG", "FI", "G", "IN") 
  tab <- tables[[2]][[1]]
  tab <- merge(tables[[2]][[1]], tables[[2]][[2]], by = 0) %>% namerows
  pr.bird <- percent.occupancy.by.guild(tab, gld, "bird", sitedat) 

#### Exploratory analysis of co-occurring competitors ####
comp_nc_ratio <- function(x, spp){
  spplists <- apply(x, 2, function(y) rownames(x)[which(y==1)])
  coocc_counts <- lapply(spplists, combn, 2) %>% lapply(t) %>% lapply(data.frame) %>% bind_rows(.id = "site") %>% group_by(X1, X2) %>% summarise(n = length(X1))
  coocc_counts$X1.diet <- spp[coocc_counts$X1,"guild"]
  coocc_counts$X2.diet <- spp[coocc_counts$X2,"guild"]
  coocc_counts$comp <- coocc_counts$X1.diet == coocc_counts$X2.diet
  #return(table(coocc_counts$comp))
  return(coocc_counts %>% split(.$comp) %>% sapply(function(x) sum(x$n)))
  #dietlists <- lapply(spplists, function(y) spp[y,"guild"])
  #dietcounts <- lapply(dietlists, table)
  #competing_coocc <- lapply(dietcounts, function(y) y*(y-1)/2) %>% sapply(sum)
  #occ <-colSums(x)
  #total_coocc <- occ*(occ-1)/2
  #nc_coocc <- total_coocc - competing_coocc
  #return(competing_coocc / nc_coocc)
}
  ratios <- purrr::map(tables, purrr::map, comp_nc_ratio, spp)
  purrr::map(ratios, purrr::map_dbl, function(x) x[2]/x[1])
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#~~~~~~~~~~~~~~~~~~END OF SCRIPT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#