### Effects of phylogenetic distance, niche overlap and habitat ###
### alteration on spatial co-occurrence patterns in Neotropical ###
### bats and birds.  

### By Aniko B. Toth, John Alroy, S Kathleen Lyons, and Andrew P. Allen

# Submitted 15 July 2024 

# Instructions #
# Create a working directory. 
# In your working directory, create the following folders:
# Raw_Data (place raw data files from Dryad in here)
# Code (place script files in here)
# Results (leave empty)
# stan (leave empty)

# Analysis script 

library(tidyverse) 
library(vegan) 
library(stringi) 
library(parallel)
library(lsa) 
library(rstan)

## Load helper functions
library(sp)
library(rlang)
library(ape)
library(rotl)
source('./Code/HelperFunctions.R')

## Prep raw data
source('./Code/Data_Prep.R')

#### Interaction Co-occurrence analysis #####
library(brms)
library(loo)
source('./Code/brms_fxns.R')

# remove abundance data
PAnb<- tobinary(PAn)

# Format data (full)
tables <- PAnb %>% map(~t(.)) %>% map(as.data.frame) %>% map(~split(., f = rownames(.) %in% unalt_sites)) %>% 
  map(map, ~t(.)) %>% map(map, clean.empty) %>% purrr::map(setNames, c("altered", "unaltered"))

# Contingency table
contables <- map(tables, map, cont_table) %>% map(bind_rows, .id = "status") %>% 
  bind_rows(.id = "taxon") %>% diet_cat(pull(spp, guild) %>% setNames(spp$unique_name), related = TRUE) %>%
  left_join(bind_rows(tax_dist)) %>% na.omit()

contables_full <- contables #save for later suppmat analyses

# Data dimensions
# species counts
contables %>% split(.$taxon) %>% 
  map(~.x %>% select(Sp1, Sp2) %>% unlist() %>% unique() %>% length())
# site and pair counts
contables %>% group_by(taxon, status) %>% 
  summarise(nsite = first(samples), npairs = n())

### Model fitting code #####
### Minimum Computing requirements: ~32GB of RAM and 4 cores.

# fit bat data
tax <- "bat"
bat_data <- stan_data_fun(filter(contables, taxon == tax))[[1]]
bat_FULL_winner <- find_best_model(tax=tax, bat_data) # runs all models, saves them and saves the loo object in current wd
summary(bat_FULL_winner)
saveRDS(bat_FULL_winner, "./Results/bat_FULL_winner.rds") # save a copy of best model to results

#bat_FULL_winner <- singlerun(bat_data, frm = 'sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ PhyloD + DietOvlp + (1 | gr(DietPair:PhyloD:Habitat:OccPair, by = DietOvlp))')

tax <- "bird"
bird_data <- stan_data_fun(filter(contables, taxon == tax))[[1]]
bird_FULL_winner <- find_best_model(tax= tax, bird_data)
summary(bird_FULL_winner)
saveRDS(bird_FULL_winner, "./Results/bird_FULL_winner.rds")

#bird_FULL_winner <- singlerun(bird_data, frm = 'sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ PhyloD + (1 | gr(DietPair:PhyloD:Habitat:OccPair, by = Habitat))')

## No-turnover models ####
shared <- tables %>% map(~.x %>% map(rownames) %>% reduce(intersect))
tables <- map2(tables, shared, function(x, y) map(x, function(z) return(z[y,])))

# Contingency table
contables <- map(tables, map, cont_table) %>% map(bind_rows, .id = "status") %>% 
  bind_rows(.id = "taxon") %>% diet_cat(pull(spp, guild) %>% setNames(spp$unique_name), related = TRUE) %>%
  left_join(bind_rows(tax_dist)) %>% na.omit()

# same code as above with new data shared tables.
tax <- "bat"
bat_data <- stan_data_fun(filter(contables, taxon == tax))[[1]]
bat_RP_winner <- find_best_model(tax, bat_data)
summary(bat_RP_winner)
saveRDS(bat_RP_winner, "./Results/bat_RP_winner.rds")

#bat_RP_winner <- singlerun(bat_data, frm = 'sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ PhyloD + DietOvlp + Habitat + (1 | gr(DietPair:PhyloD:Habitat:OccPair, by = DietOvlp)) ')


tax <- "bird"
bird_data <- stan_data_fun(filter(contables, taxon == tax))[[1]]
bird_RP_winner <- find_best_model(tax, bird_data). # runs all models, saves them and saves the loo object in current wd.
summary(bird_RP_winner)
saveRDS(bird_RP_winner, "./Results/bird_RP_winner.rds")

#bird_RP_winner <- singlerun(bird_data, frm = 'sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ PhyloD + Habitat + PhyloD:Habitat + (1 | gr(DietPair:PhyloD:Habitat:OccPair, by = DietOvlp)) ')


## Randomisations ####
library(furrr)
library(future)

# FULL MODELS
tables <- PAnb %>% map(~t(.)) %>% map(as.data.frame) %>% map(~split(., f = rownames(.) %in% unalt_sites)) %>% 
  map(map, ~t(.)) %>% map(map, clean.empty) %>% purrr::map(setNames, c("altered", "unaltered"))

plan(multicore) ## Mac/Linux
#plan(multisession) ## Windows

tables_rand <- future_map(1:100, function(y) map(tables, ~map(.x, rand_mat, i = y)) %>% 
                                               flatten() %>% setNames(c("bat$altered", "bat$unaltered", "bird$altered", "bird$unaltered")))

# list of 100 randomised contingency tables
contables_rand <- future_map(tables_rand, map, cont_table) %>% 
                      map(bind_rows, .id = "taxon_status") %>% 
                      map(~.x %>% separate(taxon_status, into = c("taxon", "status"), sep = "$") %>% 
                      diet_cat(pull(spp, guild) %>% setNames(spp$unique_name), related = TRUE) %>%
                      left_join(bind_rows(tax_dist)) %>% na.omit())


formulas <- c(bat = 'sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ PhyloD + DietOvlp + (1 | gr(DietPair:PhyloD:Habitat:OccPair, by = DietOvlp))',
              bird ='sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ PhyloD + (1 | gr(DietPair:PhyloD:Habitat:OccPair, by = Habitat))')

summary <- list()
for(i in seq_along(contables_rand)){
  summary[[i]] <- contables_rand[[i]] %>% split(.$taxon) %>% 
    map(~stan_data_fun(.x)[[1]]) %>%
    map2(formulas, ~singlerun(.x, .y) %>% summary() %>% `$`(fixed))
}

rand_out <- summary %>% map(~.x %>% map(rownames_to_column) %>% bind_rows(.id = "taxon")) %>% bind_rows(.id = "rep")
rand_out %>% group_by(taxon, rowname) %>% summarise(mean = mean(Estimate), median = median(Estimate), lowQ = quantile(Estimate, 0.25), highQ = quantile(Estimate, 0.975))

## RESTRICTED POOL model randomisations
tables <- map2(tables, shared, function(x, y) map(x, function(z) return(z[y,])))

plan(multicore) ## Mac/Linux
#plan(multisession) ## Windows

tables_rand <- future_map(1:100, function(y) map(tables, ~map(.x, rand_mat, i = y)) %>% 
                            flatten() %>% setNames(c("bat$altered", "bat$unaltered", "bird$altered", "bird$unaltered")))

# list of 100 randomised contingency tables
contables_rand <- future_map(tables_rand, map, cont_table) %>% 
  map(bind_rows, .id = "taxon_status") %>% 
  map(~.x %>% separate(taxon_status, into = c("taxon", "status"), sep = "$") %>% 
        diet_cat(pull(spp, guild) %>% setNames(spp$unique_name), related = TRUE) %>%
        left_join(bind_rows(tax_dist)) %>% na.omit())

formulas <- c(bat = 'sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ PhyloD + DietOvlp + (1 | gr(DietPair:PhyloD:Habitat:OccPair, by = DietOvlp))',
              bird ='sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ PhyloD + Habitat + (1 | gr(DietPair:PhyloD:Habitat:OccPair, by = DietOvlp))')

summary <- list()
for(i in seq_along(contables_rand)){
  summary[[i]] <- contables_rand[[i]] %>% split(.$taxon) %>% 
    map(~stan_data_fun(.x)[[1]]) %>%
    map2(formulas, ~singlerun(.x, .y) %>% summary() %>% `$`(fixed))
}

rand_out <- summary %>% map(~.x %>% map(rownames_to_column) %>% bind_rows(.id = "taxon")) %>% bind_rows(.id = "rep")
rand_out %>% group_by(taxon, rowname) %>% summarise(mean = mean(Estimate), median = median(Estimate), lowQ = quantile(Estimate, 0.25), highQ = quantile(Estimate, 0.975))


###### Supplement: classical analyses ########

###### Richness - run using raw abundance data ####
# only include species that were in the main analysis: 
PAn <- PAn %>% map(~.x[which(rownames(.x) %in% c(contables_full$Sp1, contables_full$Sp2)),])

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
rich %>% group_by(metric, taxon, status) %>% 
  summarise(median.rich = median(richness, na.rm = T)) %>% 
  pivot_wider(names_from = "status", values_from = "median.rich")

# Significance
rich %>% group_by(metric, taxon) %>% summarise(p=wilcox.test(richness~status)$p.value, 
                                               W=wilcox.test(richness~status)$statistic)

#### Beta diversity and composition analyses #####
classical_analysis <- list()

data = c("binary", "abundance")
method = c("Jaccard", "cosine")

for(d in data){
  for(m in method){
    if(d == "binary") PA <- tobinary(PAn) else if(d == "abundance") PA <- PAn
    if(m == "Jaccard"){ 
      dist <- map(PA, ~t(.)) %>% map(vegdist, method = "jaccard") # with bray-curtis index
    }else if(m == "cosine"){
      dist <- map(PA, as.matrix) %>% map(cosine) %>% map(as.dist, upper = F) %>% map(~return(1-.)) # with Ochiai index or cosine similarity
    }
    # beta diversity using betadisper
    beta <- map2(dist, PAn, function(x, y) betadisper(x, group = colnames(y) %in% unalt_sites) %>% anova)
    
    # compositional change using adonis
    data <- map(PA, ~t(.x) %>% data.frame() %>% rownames_to_column("p.sample") %>% mutate(ID = colnames(.x) %in% unalt_sites))
    comp1 <- map2(dist, data, function(x, y) adonis2(x~ID, data = y, permutations = 10000))
    
    # compositional change using canonical correspondence analysis
    dat <- map(data, ~select(.x, -p.sample, -ID))
    comp2 <- map2(dat, data, function(x, y) cca(x~ID, data = y) %>% anova.cca)
    
    classical_analysis[[paste(d, m, sep = "_")]] <- list(`beta dispersion` = beta, adonis = comp1, cca = comp2) %>% 
      map(~.x %>% map(data.frame) %>% map(rownames_to_column) %>% bind_rows(.id = "taxon") %>% na.omit() %>% 
            select(1, last_col(offset = 1), last_col()) %>% setNames(c("taxon", "F statistic", "p-value"))) %>% 
      bind_rows(.id = "analysis") %>% 
      pivot_wider(names_from = taxon, values_from = c(`F statistic`, "p-value"))
    
  }
}

classical_analysis %>% bind_rows(.id = "data_metric") %>% 
  separate(data_metric, into = c("data", "metric")) %>% 
  arrange(analysis, data) %>% relocate(analysis, .before = metric) %>% relocate(data, .before = metric) %>%
  write_csv("./Results/classical_analyses.csv")

##### occupancy calculations ######
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


#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#~~~~~~~~~~~~~~~~~~END OF SCRIPT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#