## Code to prepare raw data in ./Raw_Data folder 
## Outputs will be used in Analysis_Script.R

## Adjust relative paths to compile environment.
if(grepl(pattern = "Manuscript", getwd())){
  rt <- ".."
}else{rt <- "."}

#file.path(rt, 'Code/HelperFunctions.R') %>% source()  # also sourced from Analysis_script.R
l <- file.path(rt, "Raw_Data") %>% list.files(".txt", full.names = T)
dat <- lapply(l, read.delim) %>% setNames(file.path(rt, "Raw_Data") %>% list.files(".txt"))
message("Formatting raw data")
### update species/genera that can't be resolved ####
# Genus misspelled
dat[[6]]$genus <- dat[[6]]$genus %>% str_replace_all("Cychlaris", replacement = "Cyclarhis")
dat[[8]] <- dat[[8]] %>% filter(species != "Cychlaris gujanensis")

dat[[6]]$genus <- dat[[6]]$genus %>% str_replace_all("Phylidor", replacement = "Philydor")
dat[[8]] <- dat[[8]] %>% filter(!species %in% c("Phylidor pyrrhodes", "Phylidor erythrocercus", "Phylidor ruficaudatus"))

# Genus reclassified
dat[[6]]$genus <- dat[[6]]$genus %>% str_replace_all("Pipromorpha", replacement = "Mionectes")
dat[[6]]$species <- dat[[6]]$species %>% str_replace_all("oleaginea", replacement = "oleagineus")
dat[[8]] <- dat[[8]] %>% filter(species != "Pipromorpha oleaginea")

dat[[6]]$genus <- dat[[6]]$genus %>% str_replace_all("Teleonema", replacement = "Pipra")
dat[[8]] <- dat[[8]] %>% filter(species != "Teleonema filicauda")

dat[[6]]$genus <- dat[[6]]$genus %>% str_replace_all("Platypsaris", replacement = "Pachyramphus")
dat[[8]] <- dat[[8]] %>% mutate(species = ifelse(species == "Platypsaris minor", "Pachyramphus minor", species))

dat[[6]]$genus[which(dat[[6]]$genus == "Saltator" & dat[[6]]$species == "grossus")] <- "Pitylus"
dat[[8]] <- dat[[8]] %>% mutate(species = ifelse(species == "Orzyoborus angolensis", "Oryzoborus angolensis", species))

dat[[6]]$genus[which(dat[[6]]$genus == "Idioptilon")] <- "Hemitriccus"
dat[[6]]$species[which(dat[[6]]$genus == "Idioptilon")] <- "margaritaceiventer" 
dat[[8]] <- dat[[8]] %>% filter(!species %in% c("Idioptilon margaritaceiventer", "Idioptilon margaritaceiventris"))

dat[[6]]$genus[which(dat[[6]]$genus == "Acrochordopus")] <- "Phyllomyias"
dat[[8]] <- dat[[8]] %>% mutate(species = ifelse(species == "Acrochordopus burmeisteri", "Phyllomyias burmeisteri", species))

dat[[6]]$genus[which(dat[[6]]$genus == "Tanagra")] <- "Euphonia"
dat[[8]] <- dat[[8]] %>% filter(species != "Tanagra gouldi")

dat[[6]]$genus[which(dat[[6]]$genus == "Lanio" & dat[[6]]$species == "penicillata")] <- "Eucometis"
dat[[8]] <- dat[[8]] %>% filter(species != "Lanio penicillata")

dat[[6]]$genus[which(dat[[6]]$genus == "Hydropsalis" & dat[[6]]$species == "albicollis")] <- "Nyctidromus"
dat[[8]] <- dat[[8]] %>% filter(species != "Hydropsalis albicollis")

dat[[6]]$genus[which(dat[[6]]$genus == "Sciaphylax" & dat[[6]]$species == "hemimelaena")] <- "Myrmeciza"
dat[[8]] <- dat[[8]] %>% filter(species != "Sciaphylax hemimelaena")

dat[[6]]$genus[which(dat[[6]]$genus == "Pogonotriccus")] <- "Phylloscartes"
dat[[8]] <- dat[[8]] %>% mutate(species = ifelse(species == "Pogonotriccus ophthalmicus", "Phylloscartes ophthalmicus", species))

dat[[6]]$genus[which(dat[[6]]$genus == "Myrmelastes" & dat[[6]]$species == "hyperythrus")] <- "Myrmeciza"
dat[[8]] <- dat[[8]] %>% mutate(species = ifelse(species == "Myrmelastes hyperythrus", "Myrmeciza hyperythrus", species))

dat[[6]]$genus <- dat[[6]]$genus %>% str_replace_all("Ceratopipra", replacement = "Pipra")
dat[[8]] <- dat[[8]] %>% filter(!species %in% c("Ceratopipra erythrocephala", "Ceratopipra rubrocapilla"))

## species name updated
dat[[6]]$species[which(dat[[6]]$genus == "Xiphorhynchus" & dat[[6]]$species == "atlanticus")] <- "fuscus"
dat[[8]] <- dat[[8]] %>% filter(species != "Xiphorhynchus atlanticus")

dat[[6]]$species[which(dat[[6]]$genus == "Saltator" & dat[[6]]$species == "azarae")] <- "coerulescens"
dat[[8]] <- dat[[8]] %>% filter(species != "Saltator azarae")

#Diclidurus virgo ==> Diclidurus albus
dat[[2]]$species[grep("virgo", dat[[2]]$species)] <- "albus"
dat[[4]] <- dat[[4]] %>% filter(species != "Diclidurus virgo")

#Noctilio labialis ==> Noctilio albiventris
dat[[2]]$species[grep("labialis", dat[[2]]$species)] <- "albiventris"
dat[[4]] <- dat[[4]] %>% filter(species != "Noctilio labialis")

#Sturnira occidentalis ==> subspecies of Sturnira ludovici
dat[[2]]$species[grep("occidentalis", dat[[2]]$species)] <- "ludovici"
dat[[4]] <- dat[[4]] %>% filter(species != "Sturnira occidentalis")

#Myotis pilosatibialis is a species complex; not even monophyletic


##### Site metadata ####
sitedat <- dat[grep("samples", names(dat))] %>% setNames(c("bat", "bird")) %>% bind_rows(.id = "taxon") %>% 
  select(taxon, reference.no, sample.no, sample.name, country, ecozone, latitude, longitude, habitat, altered.habitat, MAT, MAP, richness, fragment.size)
sitedat$status[!sitedat$altered.habitat == ""] <- "Altered" 
sitedat$status[sitedat$altered.habitat == ""] <- "Unaltered" 
sitedat$siteid <- paste0("X", sitedat$sample.no)

sitedat <- sitedat %>% unite("p.sample", taxon, reference.no, latitude, longitude, habitat, altered.habitat, sep = "_", remove = F) %>% 
  mutate(p.sample = as.factor(p.sample) %>% as.numeric())


#### Species metadata ####
spp <- dat[grep("species", names(dat))] %>% bind_rows() %>% 
  select(species, life.form, individuals, samples, diet.1, diet.2, mass_g = mass..g.)

rownames(spp) <-  gsub(" ", "_", spp$species)
spp <- spp %>% mutate_if(is.character, list(~na_if(.,""))) 

### Get taxonomy and pairwise distances from OTL ####
message("Generating phylogenetic trees")
tax_otl <- get_otl_taxonomy(spp, lifeform = c("bat", "bird"), contxt = c("Mammals", "Birds"))
tax_dist <- map2(tax_otl, c("Mammals", "Birds"), get_pairwise_dist)
# drop subspp.
tax_otl <- map(tax_otl, ~.x %>% mutate(species_name = word(unique_name, 1,2, sep = " ")))

spp <- left_join(spp, bind_rows(tax_otl) %>% select(1:2), by = c("species" = "search_string"))

# Manually resolve diet conflicts in merged species -- info from Birds of the World

    # guild_conflicts <- spp %>% group_by(unique_name, life.form) %>%
    #  summarise(n = n(), diet.1 = str_unique(diet.1) %>% str_flatten_comma(na.rm = T),
    #            diet.2 = str_unique(diet.2) %>% str_flatten_comma(na.rm = T)) %>% filter(grepl(",", diet.1) | grepl(",", diet.2)) %>%
    #  na.omit()

guild_conflicts <- read_csv("./Data/guild_conflicts.csv")
spp <- spp %>% filter(!unique_name %in% guild_conflicts$unique_name) %>% rbind(guild_conflicts)  # replace conflict rows with reconciled data

## All hummingbirds to NI
hummingbirds <- read_csv("Data/hummingbirds.csv")
hummgen <- pull(hummingbirds, genus) %>% unique()
spp$diet.1[which(word(spp$unique_name, 1,1) %in% hummgen)] <- "nectarivore"
spp$diet.2[which(word(spp$unique_name, 1,1) %in% hummgen)] <- "invertivore"

## Fill missing data
missing <- read_csv("Data/missing_diet.csv")
spp <- full_join(spp, missing, by = c("unique_name")) %>%
  mutate(diet.1 = coalesce(diet.1, BotWdiet.1), 
         diet.2 = coalesce(diet.2, BotWdiet.2)) %>% 
  select(-BotWdiet.1, -BotWdiet.2)


## update waterbird diets with aquatic invertivore guild
waterbirds <- read_csv("./Data/waterbirds.csv")
spp <- spp %>% filter(!unique_name %in% waterbirds$unique_name) %>% rbind(waterbirds) 

# create guild categories ####
message("Creating dietary guild categories")
spp$diet.2[is.na(spp$diet.2)] <- spp$diet.1[is.na(spp$diet.2)]
guild_table <- spp %>% select(diet.1, diet.2) %>% unique() %>% arrange(diet.1, diet.2) %>% mutate(guild = c("A", "IA", "AO", "PA", "C", "CI", "F", "FG", "FI", "FO", "GA", "FG", "G", "IG", "GO", "IA", "CI", "FI", "IG", "I", "NI", "IO", "NI", "N", "FO", "O", "PA", "PC", "P", "S", NA)) 
spp <- spp %>% group_by(unique_name, life.form) %>% 
  summarise(individuals = sum(individuals), samples = sum(samples), mass_g = mean(mass_g, na.rm = T), 
            diet.1 = str_unique(diet.1) %>% str_flatten_comma(na.rm = T), 
            diet.2 = str_unique(diet.2) %>% str_flatten_comma(na.rm = T)) %>%
  full_join(guild_table)

#### Occurrences ####
message("Creating species-by-site matrix")
PAn <- dat[grep("register", names(dat))] %>% 
  purrr::map2(tax_otl, 
              ~filter(.x,!species %in% c("sp.", "spp.", "sp. l", "indet", "indet.")) %>% 
                mutate(name = paste(genus, species, sep = " ")) %>% 
                left_join(sitedat, by = c("sample.no", "sample.name", "country", "ecozone")) %>%
                left_join(.y, by = c("name"="search_string")) %>%
                filter(!is.na(species_name)) %>%
                pivot_wider(id_cols = "species_name", names_from = "p.sample", 
                            values_from = "count", values_fill = 0, values_fn = sum) %>% 
                data.frame(check.names = FALSE) %>% namerows()
  ) %>%
  setNames(c("bat", "bird"))


#### Calculate habitat preference ####
unalt_sites <- sitedat$p.sample[sitedat$status == "Unaltered"]

nsamp <- table(sitedat$taxon, sitedat$status)

u = purrr::map(PAn, ~apply(., 1, function(x) names(which(x > 0))) %>% 
                 sapply(function(x) length(which(x %in% (unalt_sites))))) 
a = purrr::map(PAn, ~apply(., 1, function(x) names(which(x > 0))) %>% 
                 sapply(function(x) length(which(!x %in% (unalt_sites)))))
# cosmo <- map2(u, a, data.frame) %>% 
#   purrr::map(~mutate(.,name = rownames(.))) %>% bind_rows(.id = 'taxon') %>% 
#   setNames(c("taxon", "unaltered", "altered", "name")) %>% 
#   mutate(
#     occurrences = unaltered+altered, 
#     Talt = nsamp[taxon,"Altered"],
#     Tunalt = nsamp[taxon, "Unaltered"],
#     cosmopolitan = FETmP_(Talt, Tunalt, altered, unaltered)
#   )
# 
# cosmo$group[qnorm(cosmo$cosmopolitan) < -1] <- "restricted"
# cosmo$group[qnorm(cosmo$cosmopolitan) > 1] <- "synanthropic"
# cosmo$group[is.na(cosmo$group)] <- "cosmopolitan"
# 
# cosmo$group_abbr <- substr(cosmo$group, start = 1, stop = 5)
# rownames(cosmo) <- cosmo$name

#### Match Biogeography ####
message("Matching altered and intact biogeography")
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


################
#clear unneeded objects
rm(a, u, dat, l, nsamp, rt, guild_conflicts, hummingbirds, hummgen, waterbirds, missing, keep)

