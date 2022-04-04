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
library(lsa)
## Load helper functions
source('./Code/HelperFunctions.R')

## Prep raw data
source('./Code/Data_Prep.R')


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


###### Richness - run using raw abundance data ####
# can only be run on raw data because it requires a singleton count of abundances.
rich <- list(cJ1 = map(PAn, map_dbl, cJ1rich) %>% map(~split(., names(.) %in% unalt_sites)) %>% # corrected 1st order jackknife
               map(map, cbind) %>% map(map, data.frame), 
             chao = map(PAn, map_dbl, chao1) %>% map(~split(., names(.) %in% unalt_sites)) %>%  # chao1
               map(map, cbind) %>% map(map, data.frame) , 
             fa = map(PAn, map_dbl, fisher.alpha) %>% map(~split(., names(.) %in% unalt_sites)) %>%  # fisher's alpha
               map(map, cbind) %>% map(map, data.frame)) %>% 
  map(map, bind_rows, .id = "status") %>% map(bind_rows, .id = "taxon") %>% bind_rows(.id = "metric") %>% 
  setNames(c("metric", "taxon", "status", "richness")) %>% 
  mutate(status = recode(status, `FALSE` = "Altered", `TRUE` = "Intact"))

# Summary
rich %>% group_by(metric, taxon, status) %>% summarise(mean.rich = mean(richness, na.rm = T), median.rich = median(richness, na.rm = T))
# Significance
rich %>% group_by(metric, taxon) %>% summarise(p=wilcox.test(richness~status, paired=FALSE)$p.value, 
                                               W=wilcox.test(richness~status, paired=FALSE)$statistic)

#### Beta diversity and composition analyses #####
PAnb<- tobinary(PAn)

PA <- PAnb  ## can run analyses with binary or abudance data

# with bray-curtis index
dist <- map(PA, ~t(.)) %>% map(vegdist, method = "jaccard")
# with Ochiai index or cosine similarity
dist <- map(PA, as.matrix) %>% map(cosine) %>% map(as.dist, upper = F) %>% map(~return(1-.))

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
#source('./Code/stan.R')

# Format data (full)
tables <- PAnb %>% map(~t(.)) %>% map(as.data.frame) %>% map(~split(., f = rownames(.) %in% unalt_sites)) %>% 
  map(map, ~t(.)) %>% map(map, clean.empty) %>% purrr::map(setNames, c("altered", "unaltered"))

# Contingency table
contables <- map(tables, map, cont_table) %>% map(bind_rows, .id = "status") %>% 
  bind_rows(.id = "taxon") %>% diet_cat(spp, related = TRUE) %>% na.omit()

# Calculate theta 
## Bats full
stan_data <- stan_data_fun(filter(contables, taxon == "bat"), medn = F)
bat_full <- stan_theta(stan_data)

## Birds full
stan_data <- stan_data_fun(filter(contables, taxon == "bird"), medn = F)
bat_noTurnover <- stan_theta(stan_data)

## No-turnover models ####
shared <- tables %>% map(~.x %>% map(rownames) %>% reduce(match_val))
tables <- map2(tables, shared, function(x, y) map(x, function(z) return(z[y,])))

# Contingency table
contables <- map(tables, map, cont_table) %>% map(bind_rows, .id = "status") %>% 
  bind_rows(.id = "taxon") %>% diet_cat(spp, related = TRUE) %>% na.omit()

# Calculate theta 
## Bats no turnover
stan_data <- stan_data_fun(filter(contables, taxon == "bat"), medn = F)
bird_full <- stan_theta(stan_data)

## Birds no turnover
stan_data <- stan_data_fun(filter(contables, taxon == "bird"), medn = F)
bird_noTurnover <- stan_theta(stan_data)


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