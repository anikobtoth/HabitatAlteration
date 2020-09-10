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

#### Change to presence-absence ####
PAn<- tobinary(PAn)
PAnu <- PAn %>% map(~return(.[,colnames(.) %in% unalt_sites])) # separate unaltered
PAna <- PAn %>% map(~return(.[,!colnames(.) %in% unalt_sites])) # and altered sites

#### Match Biogeography ####
meta <- sitedat %>% split(.$taxon)
coords <- map(meta, select, c(sample.no, latitude, longitude)) 

# Remove unaltered or altered sites that are not near a site of the other type. This is done to ensure the two sets have similar biogoegraphical distributions.
coords1 <- purrr::map2(coords, PAnu, function(x, y) x[x$sample.no %in% colnames(y),])
coords2 <- purrr::map2(coords, PAna, function(x, y) x[x$sample.no %in% colnames(y),])
keep <- map2(coords1, coords2, matchbiogeo) %>% map(unlist)

PAn <- map2(PAn, keep, function(x, y) return(x[,as.character(y)]))
PAn <- map(PAn, clean.empty, minrow = 1) # remove any species that now have no occurrences

# recalculate alt/unalt split from new PAn
PAnu <- PAn %>% map(~return(.[,colnames(.) %in% unalt_sites]))
PAna <- PAn %>% map(~return(.[,!colnames(.) %in% unalt_sites]))

# coords of sites we kept, for plotting later.
coords1.keep <- purrr::map2(coords, PAnu, function(x, y) x[x$sample.no %in% colnames(y),])
coords2.keep <- purrr::map2(coords, PAna, function(x, y) x[x$sample.no %in% colnames(y),])

### Beta - run using presence-absence, biogeo matched data ####

beta <- beta.types(PAn, unalt_sites)

# summary stats
b <- beta %>% group_by(taxon, unalt.pair) %>% summarise(mean = mean(Z.Score), median = median(Z.Score))
# p - values
beta %>% filter(unalt.pair != "Unaltered-Altered") %>% 
  mutate(similarity = qnorm(Z.Score)) %>% split(.$taxon) %>% 
  purrr::map(~wilcox.test(data = ., similarity~unalt.pair, paired=FALSE))

#occupancy calculations (used later)
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

#### Composition - run using presence-absence, biogeo matched data ####

dist <- map(PAn, ochiaiMatrix) %>% map(as.dist, upper = F) %>% map(~return(1-.)) 
ord <- map(dist, cmdscale, k = 2) %>% map(data.frame)
ord <- map(ord, ~merge(., y = sitedat, all.x = T, all.y = F, by.x = 0, by.y = "sample.no"))

w <- numeric()
for(j in 1:length(ord)){
  i <- ord[[j]]
  d <- numeric()
  for(r in 1:1000){
    i$rand = round(runif(n = nrow(i), min = 0, max = 1), 0)
    cent <- i %>% group_by(rand) %>% dplyr::summarise(x = mean(X1), y = mean(X2))
    d[r] <- spDists(as.matrix(cent[,c(2:3)]))[2]
  }
  cent <- i %>% group_by(status) %>% dplyr::summarise(x = mean(X1), y = mean(X2))
  obs <- spDists(as.matrix(cent[,c(2:3)]))[2]
  w[j] <- wilcox.test(d, obs)$p.value
}

# P-value of PCoA centroid shift
w


#### Format Data for Co-occurrence analysis #####
  # Remove empty rows from altered and unaltered tables
PAna <- map(PAna, clean.empty) 
PAnu <- map(PAnu, clean.empty)

# no singletons
PAn.ns <- map(PAn, clean.empty, minrow = 2)

# Format data
tables <- PAn %>% map(~t(.)) %>% map(as.data.frame) %>% map(~split(., f = rownames(.) %in% unalt_sites)) %>% 
  map(map, ~t(.)) %>% map(map, clean.empty) %>% purrr::map(setNames, c("altered", "unaltered"))

# Contingency table
contables <- map(tables, map, cont_table) %>% map(bind_rows, .id = "status") %>% 
  bind_rows(.id = "taxon") %>% diet_cat(spp, related = TRUE) %>% na.omit()

# Calculate omega 
## interaction, no bagging, no related, clipped to median
out_bat <- omega(contables, tax = "bat", type = "interaction") 
out_bird <- omega(contables, tax = "bird", type = "interaction")

##  interaction, no bagging, no related, all species
out_bat <- omega(contables, tax = "bat", type = "interaction", medn = FALSE)
out_bird <- omega(contables, tax = "bird", type = "interaction", medn = FALSE)

## interaction, bagging, no related, clipped to median
out_bat_bag <- omega(contables, tax = "bat", type = "interaction", bagging = TRUE) 
out_bird_bag <- omega(contables, tax = "bird", type = "interaction", bagging = TRUE)

## interaction, bagging, no related, all species
out_bat_bag <- omega(contables, tax = "bat", type = "interaction", bagging = TRUE, medn = FALSE) 
out_bird_bag <- omega(contables, tax = "bird", type = "interaction", bagging = TRUE, medn = FALSE)

# Code for plotting omega function output
ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch)) + 
  geom_density() + facet_wrap(parameter~., scales = "free_x")

##### Functional FETmP analysis ####
tables <- PAn.ns %>% map(~t(.)) %>% map(as.data.frame) %>% map(~split(., f = rownames(.) %in% unalt_sites)) %>% 
  map(map, ~t(.)) %>% map(map, clean.empty) %>% purrr::map(setNames, c("altered", "unaltered"))

# observed
obs <- purrr::map(tables, ~purrr::map(., function(x) simpairs(x) %>% dist2edgelist(x))) %>%
  purrr::map(bind_rows, .id = 'status') %>% bind_rows(.id = "taxon") %>% mutate(Taxon_status = paste(taxon, status, sep = "_"))

# expected
ss1 <- map(tables, map_int, ncol)  
rps <- 100

ss.bat <- lapply(1:rps, function(x) sample(c(rep(1, ss1$bat["altered"]), rep(2, ss1$bat["unaltered"]))))
ss.bird <- lapply(1:rps, function(x) sample(c(rep(1, ss1$bird["altered"]), rep(2, ss1$bird["unaltered"]))))

# no replacement site subsampling
bat.a <- lapply(ss.bat, function(x) PAn.ns[[1]][,which(x == 1)]) %>% lapply(clean.empty)
bat.u <- lapply(ss.bat, function(x) PAn.ns[[1]][,which(x == 2)]) %>% lapply(clean.empty)
bird.a <- lapply(ss.bird, function(x) PAn.ns[[2]][,which(x == 1)]) %>% lapply(clean.empty)
bird.u <- lapply(ss.bird, function(x) PAn.ns[[2]][,which(x == 2)]) %>% lapply(clean.empty)


input <- list(bat.a, bat.u, bird.a, bird.u) %>% setNames(c("bat_altered", "bat_unaltered", "bird_altered", "bird_unaltered"))

  ## CAUTION: The following code is the rate-limiting step when there are many species in the input tables (e.g. more than ~150).
  ## CAUTION: The following code will also hit the RAM limit on most computers.
  ## See below for two alternative ways to explore the results.
    # out <- map(y, map, map, simpairs)  # FETmP calculations
  
  ## The rate-limiting step can be parallelized on machines with multiple cores ##     
     cl <- makeCluster(detectCores()) 
     clusterExport(cl, c("simpairs", "dist2edgelist"))
     clusterEvalQ(cl, library(tidyverse))
  ## 
    
   exp <- purrr::map(input, ~parLapply(cl, ., fun = function(x) simpairs(x) %>% dist2edgelist(x)) %>% 
      bind_rows(.id = "subsample")) %>% bind_rows(.id = "Taxon_status")
   exp$taxon <- word(exp$Taxon_status, 1, 1, sep = "_")
   exp$status <- word(exp$Taxon_status, 2, 2, sep = "_")
   
   out <- bind_rows(list(expected = exp, observed = obs), .id = "type")
  
#### Add categorical data #####
  out <- diet_cat(out, spp, related = FALSE)
   
  out$cat.Sp1[out$Sp1 %in% uniquesp] <- "Unique"
  out$cat.Sp1[out$Sp1 %in% sharedsp] <- "Shared"
  out$cat.Sp2[out$Sp2 %in% uniquesp] <- "Unique"
  out$cat.Sp2[out$Sp2 %in% sharedsp] <- "Shared"
  
  out$cat.pair <- paste(out$cat.Sp1, out$cat.Sp2, sep = "-")
  out$cat.pair[out$cat.pair == "Shared-Unique"] <- "Unique-Shared"
  
  out$cat.group[out$cat.pair == "Unique-Unique"] <- "Unique"
  out$cat.group[out$cat.pair == "Unique-Shared"] <- "Unique"
  out$cat.group[out$cat.pair == "Shared-Shared"] <- "Shared"
  
### Cosmopolitan/restricted groupings #  
  out$cosmo.Sp1 <- cosmo[out$Sp1,"group_abbr"]
  out$cosmo.Sp2 <- cosmo[out$Sp2,"group_abbr"]
  out$cosmo.pair <- paste(out$cosmo.Sp1, out$cosmo.Sp2, sep = "-")
  
  out$cosmo.pair[out$cosmo.pair == "cosmo-restr"] <- "restr-cosmo"
  out$cosmo.pair[out$cosmo.pair == "synan-cosmo"] <- "cosmo-synan"
  out$cosmo.pair[out$cosmo.pair == "synan-restr"] <- "restr-synan"
  
#### Results collated as summary tables from output ####  
  
# d1: overall [NOT PRESENTED IN MANUSCRIPT] #####
  d1.prop.all <- out %>% group_by(subsample, Taxon_status, type) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), sd = sd(Z.Score), count = length(Z.Score))
  d1.prop.cat <- out %>% group_by(subsample, Taxon_status, cat.group, type) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), sd = sd(Z.Score), count = length(Z.Score))
  d1.mag.all <- out %>% group_by(subsample, Taxon_status, pnz = posnegzero(Z.Score), type) %>% 
    summarise(`Mean magnitude` = mean(Z.Score), count = length(Z.Score)) %>% filter(!pnz == "ZERO")
  d1.mag.cat <- out %>% group_by(subsample, Taxon_status, cat.group, pnz = posnegzero(Z.Score), type) %>% 
    summarise(`Mean magnitude` = mean(Z.Score), count = length(Z.Score)) %>% filter(!pnz == "ZERO")
  
# d2: Diet.match (competing and non-competing pairs)  ####
  # MAGNITUDES [PRESENTED IN MAIN TEXT]
    # magnitude of agg & seg, overall by diet group 
    d2.mag.all <- out %>% group_by(subsample, Taxon_status, diet.match, pnz = posnegzero(Z.Score), type) %>% 
      summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!pnz == "ZERO")
    
    # magnitude of agg & seg, by diet group and shared/unique    
    d2.mag.catp <- out %>% group_by(subsample, Taxon_status, diet.match, cat.pair, pnz = posnegzero(Z.Score), type) %>% 
      summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!pnz == "ZERO")
    
    # magnitude of agg & seg, by diet group and synanthropic/cosmpolitan/restricted  
    d5.mag.catp <- out %>% group_by(subsample, Taxon_status, diet.match, cosmo.pair, pnz = posnegzero(Z.Score), type) %>% 
      summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!pnz == "ZERO")
    
   # PROPORTIONS [PRESENTED IN SUPPLEMENT]  
    # proportion of agg:seg, overall by diet group
    d2.prop.all <- out %>% group_by(subsample, Taxon_status, diet.match, type) %>% 
      summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
    
    # proportion of agg:seg, by diet group and shared/unique
    d2.prop.catp <- out %>% group_by(subsample, Taxon_status, diet.match, cat.pair, type) %>% 
      summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
    
    # proportion of agg:seg, by diet group and synanthropic/cosmopolitan/restricted
    d5.prop.catp <- out %>% group_by(subsample, Taxon_status, diet.match, cosmo.pair, type) %>% 
      summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
    
          
# d3: Within-guild analysis: Diet.pair == "Same" [NOT PRESENTED IN MANUSCRIPT] ####
  d3.mag.all <- out[out$diet.match == "Same",]  %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type, pnz = posnegzero(Z.Score)) %>% 
    summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!pnz == "ZERO")
  d3.mag.catp <- out[out$diet.match == "Same",]  %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type, pnz = posnegzero(Z.Score), cat.pair) %>% 
    summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!pnz == "ZERO")
    
  d3.prop.all <- out[out$diet.match == "Same",] %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
  d3.prop.catp <- out[out$diet.match == "Same",] %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type, cat.pair) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
  
  
  
# GLM approach ####
  gd <- out %>% mutate(score = pnorm(Z.Score)) %>% 
    select(type, subsample, Taxon_status, id, score, diet.match, cat.pair, cosmo.pair) %>% 
    obsDexp( split.var = "type", data.var = "score", Taxon_status, id, diet.match, cat.pair, cosmo.pair)
  
  gd <- gd %>% pivot_longer(cols = c("altered", "unaltered"), names_to = "status", values_to = "score")
  bat <- gd %>% filter(taxon == "bat") %>% mutate(scaled.score = scale(score, center= FALSE), 
                                                  log.sc.score = log(scaled.score)) %>% na.omit()
  bird <- gd %>% filter(taxon == "bird") %>% mutate(scaled.score = scale(score, center= FALSE),
                                                  log.sc.score = log(scaled.score)) %>% na.omit()
  
  interaction.plot(bat$status, bat$diet.match, response = bat$log.sc.score)
  interaction.plot(bird$status, bird$diet.match, response = bird$log.sc.score)
  
  mbat <- glm(log.sc.score~diet.match*status, data = bat)
  mbird <- glm(log.sc.score~diet.match*status, data = bird)
  
  #  Anikos GLM approach (not needed anymore?) ####
  out$score <- pnorm(out$Z.Score)
  x <- out %>% group_by(subsample, Taxon_status, diet.match, cosmo.pair, type) %>%
     summarise(avmag = mean(score), count = length(score))

  obsexp <- obsDexp(x, split.var = "type", data.var = "avmag", Taxon_status, diet.match, cosmo.pair)
  glmdat <- obsexp %>% melt() %>% dcast(taxon+cosmo.pair+subsample+variable~diet.match, value.var = "value")
  
  glmdat$comp.diff <- glmdat$Same - glmdat$Different
  bat <- glmdat %>% filter(taxon == "bat")
  bird <- glmdat %>% filter(taxon == "bird")
  
  interaction.plot(bat$variable, bat$cosmo.pair, response = bat$comp.diff)
  interaction.plot(bird$variable, bird$cosmo.pair, response = bird$comp.ratio)
  
  mbat <- glm(comp.diff~variable*cosmo.pair, data = bat)
  mbird <- glm(comp.diff~variable*cosmo.pair, data = bird)
  
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
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#~~~~~~~~~~~~~~~~~~END OF SCRIPT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#