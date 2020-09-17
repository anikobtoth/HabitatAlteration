##    Habitat alteration reduces food competition and local    ###
##      functional diversity in Neotropical bats and birds.    ##### 
##                     Anikó B. Tóth                           ##
# Submitted 1. Apr 2019 
# Figures script
# All figure scripts will run if Analysis_Script.R is run first;
# Please note warnings about data file sizes. 

library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
##  MAIN TEXT 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Figure 1: PCoA ####
#ordinations, colored by alteration type
p1 <- ggplot(ord[[1]], aes(x = X1, y = X2, col = status)) + geom_point() + labs(x = "Comp1", y = "Comp2", col = "alteration") + theme(axis.text = element_text(size = 12)) + labs(col = "status") + 
  scale_colour_manual(values = c("red", "dodgerblue"))
p2 <- ggplot(ord[[2]], aes(x = X1, y = X2, col = status)) + geom_point() + labs(x = "Comp1", y = "Comp2", col = "alteration") + theme(axis.text = element_text(size = 12)) + labs(col = "status") + 
  scale_colour_manual(values = c("red", "dodgerblue")) 
p <- plot_grid(p1+theme(legend.position = "none"), 
               p2+theme(legend.position = "none", axis.title.y = element_blank()), ncol = 2, labels = c("A", "B"), 
               align = "lr")
plot_grid(p, get_legend(p1), rel_widths = c(1, .2))

## specs #################
colors <- c("#FFCB5C", "#FF791F", "#2F99DC", "#084887")
specs <- list(geom_density(lwd = 1.4, alpha = .3),
              facet_wrap(parameter~., scales = "free_x"), 
              scale_color_manual(values = colors), 
              scale_fill_manual(values = colors))

# FIGURE 2: bat results##############
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bat_medianFALSE_100reps_interaction_baggingFALSE.RData")
a1 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + specs
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bat_medianTRUE_100reps_interaction_baggingFALSE.RData")
a2 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bat_medianFALSE_100reps_interaction_baggingTRUE.RData")
a3 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bat_medianTRUE_100reps_interaction_baggingTRUE.RData")
a4 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")

legend <- get_legend(a1+ theme(legend.box.margin = margin(0, 0, 0, 1)))

p1 <- plot_grid(a1 + theme(legend.position = "none"), a2, a3, a4, align = "hv", nrow = 2)
plot_grid(p1, legend, nrow = 1, rel_widths = c(2, 0.5))

# FIGURE 3: bird results################
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bird_medianTRUE_100reps_interaction_baggingFALSE.RData")
r1 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bird_medianFALSE_100reps_interaction_baggingFALSE.RData")
r2 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bird_medianTRUE_100reps_interaction_baggingTRUE.RData")
r3 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bird_medianFALSE_100reps_interaction_baggingTRUE.RData")
r4 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")

legend <- get_legend(r1+ theme(legend.box.margin = margin(0, 0, 0, 7)))

p1 <- plot_grid(r1 + theme(legend.position = "none"), r2, r3, r4, align = "hv", nrow = 2)
plot_grid(p1, legend, nrow = 1, rel_widths = c(2, 0.5))



# Figs 2 and 3 with extra labels: ################
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bat_medianFALSE_100reps_interaction_baggingFALSE.RData")
bat_FF <- out
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bat_medianFALSE_100reps_interaction_baggingTRUE.RData")
bat_FT <- out
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bat_medianTRUE_100reps_interaction_baggingFALSE.RData")
bat_TF <- out
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bat_medianTRUE_100reps_interaction_baggingTRUE.RData")
bat_TT <- out

resbat <- bind_rows(list("common_only-bagged" = bat_TT, "common_only-not_bagged" = bat_TF, "all-bagged" = bat_FT, "all-not_bagged" = bat_FF), .id = "model") %>% 
  separate(model, into = c("species", "bagging"), sep = "-") %>% 
  mutate(species = factor(species, levels = c("common_only", "all")), 
         bagging = factor(bagging, levels = c("not_bagged", "bagged")))
ggplot(resbat %>% filter(parameter != "deviance"), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  geom_density(alpha = 0.25, lwd = 1) + 
  facet_wrap(bagging~species+parameter, scales = "free", nrow = 2, ncol = 4) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + ggtitle("Bats")


load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bird_medianFALSE_100reps_interaction_baggingFALSE.RData")
bird_FF <- out
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bird_medianFALSE_100reps_interaction_baggingTRUE.RData")
bird_FT <- out
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bird_medianTRUE_100reps_interaction_baggingFALSE.RData")
bird_TF <- out
load("E:/HabitatAlteration/Results/Omega/interaction/out_jagsfit_bird_medianTRUE_100reps_interaction_baggingTRUE.RData")
bird_TT <- out

resbird <- bind_rows(list("common_only-bagged" = bird_TT, "common_only-not_bagged" = bird_TF, "all-bagged" = bird_FT, "all-not_bagged" = bird_FF), .id = "model") %>% 
  separate(model, into = c("species", "bagging"), sep = "-") %>% 
  mutate(species = factor(species, levels = c("common_only", "all")), 
         bagging = factor(bagging, levels = c("not_bagged", "bagged")))

ggplot(resbird %>% filter(parameter != "deviance"), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  geom_density(alpha = 0.25, lwd = 1) + 
  facet_wrap(bagging~species+parameter, scales = "free", nrow = 2, ncol = 4) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + ggtitle("Birds")





######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
##  SUPPLEMENT 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### FIGURE S1: Maps ####
par(mfrow = c(2,2), mar = c(0,0,0,0), oma = c(1,1,1,1))

maps::map("world", xlim = c(-125, -30), ylim = c(-40, 25))
points(latitude~longitude, coords1$bat, pch = 16, cex=.7, col = "dodgerblue")
points(latitude~longitude, coords2$bat, pch = 16, cex=.7, col = "red")
mtext("Bats-original sites", side = 3, line = 1)
box()

maps::map("world", xlim = c(-125, -30), ylim = c(-40, 25))
points(latitude~longitude, coords1$bird, pch = 16, cex=.7, col = "dodgerblue")
points(latitude~longitude, coords2$bird, pch = 16, cex=.7, col = "red")
mtext("Birds-original sites", side = 3, line = 1)
box()

maps::map("world", xlim = c(-125, -30), ylim = c(-40, 25))
points(latitude~longitude, coords1.keep$bat, pch = 16, cex=.7, col = "dodgerblue")
points(latitude~longitude, coords2.keep$bat, pch = 16, cex=.7, col = "red")
mtext("Bats-included sites", side = 3, line = 1)
box()

maps::map("world", xlim = c(-125, -30), ylim = c(-40, 25))
points(latitude~longitude, coords1.keep$bird, pch = 16, cex=.7, col = "dodgerblue")
points(latitude~longitude, coords2.keep$bird, pch = 16, cex=.7, col = "red")
mtext("Birds-included sites", side = 3, line = 1)
box()


## FIGURE S6: Altered habitat types ########

pbat <- ggplot(pr.bat, aes(x = altered_habitat, y = value*100, fill= altered_habitat)) + geom_boxplot(outlier.size = 0.3) +
  facet_wrap(variable~., nrow = 4, ncol = 1, scales = "free_y") + geom_jitter(width = 0.15, size = 0.3) +
  scale_fill_manual(values = c(gg_color_hue(9), "gray50")) +
  theme(axis.text.x = element_blank()) + labs(y = "Percent of species co-occurring", x = "Habitat type")


pbird <- ggplot(pr.bird, aes(x = altered_habitat, y = value*100, fill= altered_habitat)) + geom_boxplot(outlier.size = 0.3) +
  facet_wrap(variable~., nrow = 4, ncol = 2, scales = "free_y") + geom_jitter(width = 0.15, size = 0.3) +
  theme(axis.text.x = element_blank()) + labs(y = "Percent of species co-occurring", x = "Habitat type", fill = "Type of alteration") + 
  scale_fill_manual(values = c(gg_color_hue(9), "gray50"), 
                    limits = c("combined", "cropland", "disturbed forest", "fragment", "inhabited area", "pasture", "plantation", "secondary forest", "unaltered"),
                    labels = c("combined", "cropland", "disturbed forest", "fragment", "inhabited area", "pasture", "plantation", "secondary forest", "unaltered"))

plot_grid(pbird+ theme(legend.position = c(1, 0), legend.justification = c(1, 0), legend.text = element_text(size = 10)) + guides(fill=guide_legend(ncol=2)), 
          pbat+theme(legend.position="none"), ncol = 2, rel_widths = c(1.9, 1), labels = c("Birds", "Bats"))




### GUILD-BY-GUILD [NOT PRESENTED IN MANUSCRIPT] #####
# prop agg
d <- d3.prop.cat
d <- d %>% filter(!is.na(diet.pair))
 keep <- d %>% group_by(diet.pair, Taxon_status) %>% summarise(avg = mean(count), mean(agg)) %>% filter(avg>15) 
 keep <- names(table(keep$diet.pair)[table(keep$diet.pair) >1])
 d <- d[d$diet.pair %in% keep,]

o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.pair", "cat.group"), all = T)
pAgg <- m %>% select(Taxon_status, diet.pair, cat.group, subsample = subsample.x, expected = agg.x, observed = agg.y) %>% 
  mutate(obsDexp = observed/expected)
pAgg$taxon <- word(pAgg$Taxon_status, 1, 1, sep = "_")
pAgg$status <- word(pAgg$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.pair+cat.group~status, value.var = 'obsDexp', data = pAgg)
obsexp <- obsexp %>% filter(is.finite(altered), is.finite(unaltered)) %>% na.omit()
ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.pair, fill = diet.pair)) + 
  geom_abline(slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
  geom_hline(yintercept = 1, col = "gray") +
  geom_vline(xintercept = 1, col = "gray") +
  geom_point(size = 0.5) + 
  stat_ellipse(geom= "polygon", alpha = 0.1) +
  facet_grid(taxon~cat.group, scales = "fixed") + 
  scale_colour_hue(h = c(270, 20, 330, 60, 300, 110, 216), c = 100, l = c(50, 60, 90, 50, 80, 85, 70)) +
  scale_fill_hue(h = c(270, 20, 330, 60, 300, 110, 216), c = 100, l = c(50, 60, 90, 50, 80, 85, 70))

## magnitude 
d <- d3.mag.cat
d <- d %>% filter(!is.na(diet.pair))
keep <- d %>% group_by(diet.pair, Taxon_status) %>% summarise(avg = mean(count), mean(avmag)) %>% filter(avg>15) 
keep <- names(table(keep$diet.pair)[table(keep$diet.pair) >1])
d <- d[d$diet.pair %in% keep,]

o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.pair", "cat.group", "pnz"), all = T)

avmag <- m %>% select(Taxon_status, diet.pair, cat.group, pnz, subsample = subsample.x, expected = avmag.x, observed = avmag.y) %>% 
  mutate(obsDexp = observed/expected)
avmag$taxon <- word(avmag$Taxon_status, 1, 1, sep = "_")
avmag$status <- word(avmag$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.pair+cat.group+pnz~status, value.var = 'obsDexp', data = avmag, fun.aggregate = median)
obsexp <- obsexp %>% filter(is.finite(altered), is.finite(unaltered)) %>% na.omit()

ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.pair, fill = diet.pair)) + 
  geom_abline(slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
  geom_hline(yintercept = 1, col = "gray") +
  geom_vline(xintercept = 1, col = "gray") +
  geom_point(size = 0.5) + 
  stat_ellipse(geom= "polygon", alpha = 0.1) +
  facet_grid(taxon+pnz~cat.group, scales = "free") + 
  scale_colour_hue(h = c(270, 20, 330, 60, 300, 110, 216), c = 100, l = c(50, 60, 90, 50, 80, 85, 70)) +
  scale_fill_hue(h = c(270, 20, 330, 60, 300, 110, 216), c = 100, l = c(50, 60, 90, 50, 80, 85, 70))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
## OLD STUFF 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### DIET OVERLAP ####

colors2 <- c("#34008C", "#E8AC00")
specs <- list(geom_hline(yintercept = 1, col = "gray"),
              geom_vline(xintercept = 1, col = "gray"),
              stat_ellipse(geom= "polygon", alpha = 0.3, aes(col = diet.match, fill = diet.match)),
              geom_point(aes(col = diet.match, fill = diet.match), size = 0.5), 
              geom_rug(aes(col = diet.match), alpha = 0.5),
              scale_colour_manual(values = colors2),
              scale_fill_manual(values = colors2),
              theme(legend.position = "none")
)
# FIGURE 2: Mag.all ####
d <- d2.mag.all
obsexp <- obsDexp(d, split.var = "type", data.var = "avmag", Taxon_status, diet.match, pnz)
#b <- bayesPairedTtest(obsexp, split.var = "diet.match", taxon, pnz)

anno <- obsexp %>% group_by(taxon, pnz) %>% 
  summarise() %>% ungroup() %>% 
  mutate(label = paste0("  ", LETTERS[1:4]), x = -Inf, y = Inf)

ggplot(obsexp, aes(x = unaltered, y = altered)) +
  specs + facet_wrap(pnz~taxon, scales = "free") +
  geom_text(data = anno, aes(x = x, y = y, label = label), vjust = "inward", hjust = "inward", fontface = 2)

# FIGURE 3: Mag shared/unique cat.pair ####

d <- d2.mag.catp
obsexp <- obsDexp(d, split.var = "type", data.var = "avmag", Taxon_status, diet.match, cat.pair, pnz)

# bayesian paired t-test on competing and non-competing pairs
#b <- bayesPairedTtest(obsexp, split.var = "diet.match", taxon, pnz, cat.pair)

anno <- obsexp %>% group_by(cat.pair, taxon, pnz) %>% 
  summarise() %>% ungroup() %>% 
  mutate(label = paste0("  ", LETTERS[1:12]), x = -Inf, y = Inf)

ggplot(obsexp, aes(x = unaltered, y = altered)) + 
  specs+ facet_wrap(cat.pair+taxon~pnz, scales = "free") +
  geom_text(data = anno, aes(x = x, y = y, label = label), vjust = "inward", hjust = "inward", fontface = 2)

# FIGURES 4 and 5: Mag syn/cosmo/restr cosmo.pair ####

d <- d5.mag.catp
obsexp <- obsDexp(d, split.var = "type", data.var = "avmag", Taxon_status, diet.match, cosmo.pair, pnz)
# bayesian paired t-test on competing and non-competing pairs
#b <- bayesPairedTtest(obsexp, split.var = "diet.match", taxon, pnz, cosmo.pair)

obsexp$cosmo.pair <- factor(obsexp$ cosmo.pair, levels = c("synan-synan", "cosmo-synan", "cosmo-cosmo", "restr-synan", "restr-cosmo", "restr-restr"))
anno <- obsexp %>% filter(taxon == "bat") %>% group_by(pnz, cosmo.pair) %>% 
  summarise() %>% ungroup() %>% 
  mutate(label = paste0("  ", LETTERS[1:12]), x = -Inf, y = Inf)

ggplot(obsexp %>% filter(taxon == "bat"), aes(x = unaltered, y = altered)) + 
  specs+ facet_wrap(cosmo.pair~pnz, scales = "free", ncol = 2, nrow = 6) +
  geom_text(data = anno, aes(x = x, y = y, label = label), vjust = "inward", hjust = "inward", fontface = 2)

anno <- obsexp %>% filter(taxon == "bird") %>% group_by(pnz, cosmo.pair) %>% 
  summarise() %>% ungroup() %>% 
  mutate(label = paste0("  ", LETTERS[1:12]), x = -Inf, y = Inf)

ggplot(obsexp %>% filter(taxon == "bird"), aes(x = unaltered, y = altered)) + 
  specs+ facet_wrap(cosmo.pair~pnz, scales = "free", ncol = 2, nrow = 6)+
  geom_text(data = anno, aes(x = x, y = y, label = label), vjust = "inward", hjust = "inward", fontface = 2)


### Supplement old stuff
### FIGURE S2: Explanation for plots #####
meansx = c(0.8, 1, 1.2, 1.2, 1.2, 1, 0.8, 0.8)
meansy = c(0.8, 0.8, 0.8, 1, 1.2, 1.2, 1.2, 1 )
d <- map2(meansx, meansy, function(x, y) data.frame(unaltered=rnorm(n= 100, mean = x, sd = 0.2), altered = rnorm(n = 100, mean = y, sd = 0.2))) %>% 
  setNames(c("bl", "bc", "br", "mr", "tr", "tc", "tl", "ml")) %>% bind_rows(.id = "cluster")
d$group[d$cluster %in% c("bl", "tr")] <- "Both-Unchanged"
d$group[d$cluster %in% c("br", "tl")] <- "Both-Reversed"
d$group[d$cluster %in% c("bc", "tc")] <- "Altered"
d$group[d$cluster %in% c("mr", "ml")] <- "Unaltered"

d$diet[d$cluster %in% c("bc", "br", "bl", "ml")] <- "Different"
d$diet[d$cluster %in% c("mr", "tc", "tl", "tr")] <- "Same"

ggplot(d, aes(x = unaltered, y = altered)) + 
  geom_point(aes(col = diet, fill = diet), size = 0.8) + stat_ellipse(geom= "polygon", alpha = 0.3, aes(col = diet, fill = diet)) + 
  facet_wrap(~group, scales = "free") +
  geom_rug(aes(col = diet), alpha = 0.5) + scale_color_manual(values = colors2) +
  scale_fill_manual(values = colors2)

## FIGURE S3: Prop.all ####
d <- d2.prop.all
obsexp <- obsDexp(d, split.var = "type",data.var = "agg", Taxon_status, diet.match)
#b <- bayesPairedTtest(obsexp, split.var = "diet.match", taxon)

ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.match, fill = diet.match)) + 
  specs+ 
  facet_grid(taxon~., scales = "fixed") + geom_abline(slope = 1, intercept = 0, lty = 2, col = "gray")

## FIGURE S4: Prop shared/unique cat.pair #### 
d <- d2.prop.catp
obsexp <- obsDexp(d, split.var = "type",data.var = "agg", Taxon_status, diet.match, cat.pair) %>% 
  mutate(pnz = "Proportion")
#b <- bayesPairedTtest(obsexp, split.var = "diet.match", taxon, pnz, cat.pair)

ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.match, fill = diet.match)) + 
  specs+  
  facet_wrap(taxon~cat.pair, scales = "free") 

## FIGURE S5: Prop syn/cosmo/rest cosmo.pair #### 
d <- d5.prop.catp
obsexp <- obsDexp(d, split.var = "type",data.var = "agg", Taxon_status, diet.match, cosmo.pair) %>% 
  mutate(pnz = "Proportion")
b <- bayesPairedTtest(obsexp, split.var = "diet.match", taxon, pnz, cosmo.pair)

obsexp$cosmo.pair <- factor(obsexp$ cosmo.pair, levels = c("synan-synan", "cosmo-synan", "cosmo-cosmo", "restr-synan", "restr-cosmo", "restr-restr"))
ggplot(obsexp, aes(x = unaltered, y = altered)) + 
  specs+  
  facet_wrap(taxon~cosmo.pair, scales = "free", ncol = 6, nrow = 2) 


# TABLE S1: Sig tests ####
l <- list.files("./Results/Bayes_Paired_Ttests", full.names = T)
s <- list()
t <- list()
u <- list()
for(i in seq_along(l)){
  load(l[[i]])
  s[[i]] <-  b %>% purrr::map(~.$stats %>% data.frame) %>% purrr::map(~.["mu_diff", c("mean", "sd", "X..comp")]) %>% bind_rows(.id = "id")
  t[[i]] <-  b %>% purrr::map(~.$stats %>% data.frame) %>% purrr::map(~.["sigma_diff", c("mean", "sd", "X..comp")]) %>% bind_rows(.id = "id")
  u[[i]] <-  b %>% purrr::map(~.$stats %>% data.frame) %>% purrr::map(~.["eff_size", c("mean", "sd", "X..comp")]) %>% bind_rows(.id = "id")
}

sigtests <- map2(s, t, merge, by="id", all = TRUE)
sigtests <- map2(sigtests, u, merge, by = "id", all = TRUE)
sigtests <- bind_rows(sigtests) %>% setNames(c("id", "mean_mu_diff", "sd_mu_diff", "p_mu_diff", "mean_sd_diff", "sd_sd_diff", "p_sd_diff", "mean_eff_size", "sd_eff_size", "p_eff_size"))
sigtests$id <- substr(sigtests$id, start = 1, stop = nchar(sigtests$id)-1)
sigtests <- sigtests %>% separate(id, c("status", "taxon", "type", "group"), sep = "_")
sigtests$type[sigtests$type == ""] <- "Proportion"
sigtests$group[sigtests$group %>% is.na] <- "all"
sigtests$group[sigtests$group == ""] <- "all"

#####

#####