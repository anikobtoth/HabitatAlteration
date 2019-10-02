##    Habitat alteration reduces food competition and local    ###
##      functional diversity in Neotropical bats and birds.    ### 
##                     Anikó B. Tóth                           ##

# Submitted 18 July 2019 
# Analysis script 

source('./Code/HelperFunctions.R')
source('./Code/Forbes_Index.R') # Alroy 2018

library("tidyverse")
library("vegan")
library("reshape2")
library("stringi")
library("parallel")

# Load Datasets ####
# settings
options(stringsAsFactors = FALSE)
# abundance data (Neotropical bats and birds, raw abundance data)
load("./Data/PA1_raw_abund_all.RData")
PAn <- map(PAn, clean.empty) # remove empty rows.
# species metadata
load("./Data/Species_metadata.RData")
# site metadata 
load("./Data/sitedat_metadata.Rdata")

# Site data #
#load("~/Desktop/MacQuarie_PhD/Thesis/Chapter_3/R_Files/Data/sitedat_meta_clim.RData")

neotropics <- sitedat[sitedat$ecozone %in% c("Neotropic"),]$siteid
unalt_sites <- sitedat[sitedat$altered_habitat == "",]$sample_no 
unalt_sites <- paste("X", unalt_sites, sep = "")
sitedat$status[sitedat$altered_habitat == ""] <- "Unaltered" 
sitedat$status[!sitedat$altered_habitat == ""] <- "Altered" 

  ##### Guild Categories #####
  
  # redo guild categories
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
  
  
###### Richness - run using raw abundance data ####
# can only be run on raw data because it requires a singleton count of abundances.
# SQUARES
squares <- map(PAn, map_dbl, rscale) %>% map(~split(., names(.) %in% unalt_sites)) %>% 
  map(map, cbind) %>% map(map, data.frame) %>% map(bind_rows, .id = "status") %>% 
  bind_rows(.id = "taxon") %>% setNames(c("taxon", "status", "richness"))
squares$status <- plyr::revalue(squares$status, c("FALSE" = "Altered", "TRUE" = "Unaltered"))
# summary stats
squares %>% group_by(taxon, status) %>% summarise(mean.rich = mean(richness), median.rich = median(richness))

# test for significant difference between altered and unaltered
squares %>% group_by(taxon) %>% summarise(w=wilcox.test(richness~status, paired=FALSE)$p.value)

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
coords <- map(meta, select, c(siteid, latitude, longitude)) 

# Remove unaltered or altered sites that are not near a site of the other type. This is done to ensure the two sets have similar biogoegraphical distributions.
coords1 <- purrr::map2(coords, PAnu, function(x, y) x[x$siteid %in% colnames(y),])
coords2 <- purrr::map2(coords, PAna, function(x, y) x[x$siteid %in% colnames(y),])
keep <- map2(coords1, coords2, matchbiogeo) %>% map(unlist)

PAn <- map2(PAn, keep, function(x, y) return(x[,y]))
PAn <- map(PAn, clean.empty, minrow = 1) # remove any species that now have no occurrences

# recalculate alt/unalt split from new PAn
PAnu <- PAn %>% map(~return(.[,colnames(.) %in% unalt_sites]))
PAna <- PAn %>% map(~return(.[,!colnames(.) %in% unalt_sites]))

# coords of sites we kept, for plotting later.
coords1.keep <- purrr::map2(coords, PAnu, function(x, y) x[x$siteid %in% colnames(y),])
coords2.keep <- purrr::map2(coords, PAna, function(x, y) x[x$siteid %in% colnames(y),])

### Beta - run using presence-absence, biogeo matched data ####

beta <- beta.types(PAn, unalt_sites)

# summary stats
b <- beta %>% group_by(taxon, unalt.pair) %>% summarise(mean = mean(Z.Score), median = median(Z.Score))
# p - values
beta %>% filter(unalt.pair != "Unaltered-Altered") %>% group_by(taxon) %>% summarise(w=wilcox.test(Z.Score[unalt.pair=="Altered-Altered"], Z.Score[unalt.pair=="Unaltered-Unaltered"], paired=FALSE)$p.value)

#occupancy calculations (used later)
occ.count <- map(PAn, t)  %>% map(data.frame) %>% map(~split(., rownames(.) %in% unalt_sites)) %>% 
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

dist <- map(PAn, forbesMatrix) %>% map(as.dist, upper = F) %>% map(~return(1-.)) 
ord <- map(dist, cmdscale, k = 2) %>% map(data.frame)
ord <- map(ord, ~merge(., y = sitedat[,c(8:16)], all.x = T, all.y = F, by.x = 0, by.y = "siteid"))

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


#### Select shared species #####
  # Remove empty rows from altered and unaltered tables
PAna <- map(PAna, clean.empty) 
PAnu <- map(PAnu, clean.empty)

  # Limit PAn to shared species. Uncomment to run analysis on shared pairs only. 
#PAn  <- map2(PAn, map2(lapply(PAna, rownames), lapply(PAnu, rownames), function(x, y) x[which(x %in% y)]), function(x,y) return(x[y,]))
#PAnu <- PAn %>% map(~return(.[,colnames(.) %in% unalt_sites]))
#PAna <- PAn %>% map(~return(.[,!colnames(.) %in% unalt_sites]))

##### Functional Altered/unaltered analysis ####
# Format data
tables <- PAn %>% map(~t(.)) %>% map(as.data.frame) %>% map(~split(., f = rownames(.) %in% unalt_sites)) %>% 
  map(map, ~t(.)) %>% map(map, clean.empty) %>% purrr::map(setNames, c("altered", "unaltered"))
# observed
obs <- purrr::map(tables, ~purrr::map(., function(x) simpairs(x) %>% dist2edgelist(x))) %>%
  purrr::map(bind_rows, .id = 'status') %>% bind_rows(.id = "taxon") %>% mutate(Taxon_status = paste(taxon, status, sep = "_"))

# expected
ss1 <- map(tables, map_int, ncol)  
rps <- 1000 

nbat.a <- resamp(PAn[[1]], reps = rps, sites = ss1$bat["altered"]) 
bat.u <- resamp(PAn[[1]], reps = rps, sites = ss1$bat["unaltered"])
bird.a <- resamp(PAn[[2]], reps = rps, sites = ss1$bird["altered"])
bird.u <- resamp(PAn[[2]], reps = rps, sites = ss1$bird["unaltered"])

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
     
 ### Collate results ###
     ## For the output of this analysis 
  # with reps = 100 (bird results are randomly subsampled to reduce data file size), 
  # load results/output tables with the following code: 
  load("./Results/out_bat_100r_PAn3.RData")
  load("./Results/out_bird_randsamp_100r_PAn3_1.7m.RData")
  
  ## ***NOTE*** the output tables take up roughly 6 gigabytes of memory. 
  
#### Add diet and shared/unique categories #####
  out$diet.Sp1 <- spp[out$Sp1,"guild"]
  out$diet.Sp2 <- spp[out$Sp2,"guild"]
  out$diet.pair <- map2(out$diet.Sp1, out$diet.Sp2, function(x, y) c(x,y)) %>% map(sort) %>% map(paste, collapse = "-") %>% unlist()
  related <- c("C-CI", "CI-IG", "NF-NI", "I-NI", "CI-NI", "CI-I", "CI-FI", "CI-IN","F-FN", "FI-NF", "FI-NI", 
              "FN-IN", "FI-I", "FI-FN", "FI-IN", "FN-N", "FG-FI", "FG-FN", "FG-IG", "FI-IG", "FG-G", "I-IN", "IN-N", 
              "N-NI", "N-NF", "F-FI", "FG-NF", "F-FG", "G-IG", "I-IG","IG-IN")  # pairs that share one diet source.
  out$diet.match <- as.numeric(out$diet.Sp1 == out$diet.Sp2)*2
  out$diet.match[out$diet.pair %in% related] <- 1
  out$diet.match[out$diet.match == 0] <- "Different"
  out$diet.match[out$diet.match == 1] <- "Related"
  out$diet.match[out$diet.match == 2] <- "Same"
  
  out$cat.Sp1[out$Sp1 %in% uniquesp] <- "Unique"
  out$cat.Sp1[out$Sp1 %in% sharedsp] <- "Shared"
  out$cat.Sp2[out$Sp2 %in% uniquesp] <- "Unique"
  out$cat.Sp2[out$Sp2 %in% sharedsp] <- "Shared"
  
  out$cat.pair <- paste(out$cat.Sp1, out$cat.Sp2, sep = "-")
  out$cat.pair[out$cat.pair == "Shared-Unique"] <- "Unique-Shared"
  
  out$cat.group[out$cat.pair == "Unique-Unique"] <- "Unique"
  out$cat.group[out$cat.pair == "Unique-Shared"] <- "Unique"
  out$cat.group[out$cat.pair == "Shared-Shared"] <- "Shared"

#### Results collated as summary tables from output ####  
  
  # d1: overall [NOT PRESENTED IN MANUSCRIPT]
  d1.prop.all <- out %>% group_by(subsample, Taxon_status, type) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), sd = sd(Z.Score), count = length(Z.Score))
  d1.prop.cat <- out %>% group_by(subsample, Taxon_status, cat.group, type) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), sd = sd(Z.Score), count = length(Z.Score))
  d1.mag.all <- out %>% group_by(subsample, Taxon_status, posnegzero(Z.Score), type) %>% 
    summarise(`Mean magnitude` = mean(Z.Score), count = length(Z.Score)) %>% filter(!`posnegzero(Z.Score)` == "ZERO")
  d1.mag.cat <- out %>% group_by(subsample, Taxon_status, cat.group, posnegzero(Z.Score), type) %>% 
    summarise(`Mean magnitude` = mean(Z.Score), count = length(Z.Score)) %>% filter(!`posnegzero(Z.Score)` == "ZERO")
  
  # d2: Diet.match (same, related, and different pairs) [PRESENTED IN MAIN TEXT]
    # proportion of agg:seg, overall by diet group
    d2.prop.all <- out %>% group_by(subsample, Taxon_status, diet.match, type) %>% 
      summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
    
    # proportion of agg:seg, by diet group and shared/unique
    d2.prop.cat <- out %>% group_by(subsample, Taxon_status, diet.match, cat.group, type) %>% 
      summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
    
    # magnitude of agg & seg, overall by diet group 
    d2.mag.all <- out %>% group_by(subsample, Taxon_status, diet.match, posnegzero(Z.Score), type) %>% 
      summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!`posnegzero(Z.Score)` == "ZERO")
    
    # magnitude of agg & seg, by diet group and shared/unique    
    d2.mag.cat <- out %>% group_by(subsample, Taxon_status, diet.match, cat.group, posnegzero(Z.Score), type) %>% 
      summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!`posnegzero(Z.Score)` == "ZERO")
  
 
       
  # d3: Within-guild analysis: Diet.pair == "Same" [PRESENTED IN SUPPLEMENT]
  d3.prop.all <- out[out$diet.match == "Same",] %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
  d3.prop.cat <- out[out$diet.match == "Same",] %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type, cat.group) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
  d3.mag.all <- out[out$diet.match == "Same",]  %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type, posnegzero(Z.Score)) %>% 
    summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!`posnegzero(Z.Score)` == "ZERO")
  d3.mag.cat <- out[out$diet.match == "Same",]  %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type, posnegzero(Z.Score), cat.group) %>% 
    summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!`posnegzero(Z.Score)` == "ZERO")

  # d4: Related pairs analysis: Diet.pair == "Related"  [NOT PRESENTED IN MANUSCRIPT]
  d4.prop.all <- out[out$diet.match == "Related",] %>% na.omit() %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
  d4.prop.cat <- out[out$diet.match == "Related",] %>% na.omit() %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type, cat.group) %>% 
    summarise(seg = percneg(Z.Score), agg = percpos(Z.Score), count = length(Z.Score))
  d4.mag.all <- out[out$diet.match == "Related",] %>% na.omit() %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type, posnegzero(Z.Score)) %>% 
    summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!`posnegzero(Z.Score)` == "ZERO")
  d4.mag.cat <- out[out$diet.match == "Related",] %>% na.omit() %>% group_by(subsample, Taxon_status, diet.pair, diet.match, type, posnegzero(Z.Score), cat.group) %>% 
    summarise(avmag = mean(Z.Score), count = length(Z.Score)) %>% filter(!`posnegzero(Z.Score)` == "ZERO")
  
#

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
  
  
  