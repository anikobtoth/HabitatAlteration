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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

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




### DIET OVERLAP ####

colors2 <- c("#34008C", "#E8AC00")
specs <- list(geom_hline(yintercept = 1, col = "gray"),
              geom_vline(xintercept = 1, col = "gray"),
              stat_ellipse(geom= "polygon", alpha = 0.3, aes(col = diet.match, fill = diet.match)),
              geom_point(aes(col = diet.match, fill = diet.match), size = 0.5), 
              geom_rug(aes(col = diet.match), alpha = 0.5),
              #geom_point(aes(x = 1, y = 1), col = "black", size = 1.2),
              scale_colour_manual(values = colors2),
              scale_fill_manual(values = colors2),
              theme(legend.position = "none"),
              geom_text(data = anno, aes(x = x, y = y, label = label), vjust = "inward", hjust = "inward", fontface = 2)
)
# FIGURE 2: Mag.all ####
d <- d2.mag.all

d <- d %>% filter(!is.na(diet.match))
o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.match", "pnz"), all = T)
avmag <- m %>% select(Taxon_status, diet.match, pnz, subsample = subsample.x, expected = avmag.x, observed = avmag.y) %>% 
  mutate(obsDexp = observed/expected)
avmag$taxon <- word(avmag$Taxon_status, 1, 1, sep = "_")
avmag$status <- word(avmag$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.match+pnz~status, value.var = 'obsDexp', data = avmag, fun.aggregate = median)

anno <- obsexp %>% group_by(taxon, pnz) %>% 
  summarise() %>% ungroup() %>% 
  mutate(label = paste0("  ", LETTERS[1:4]), x = -Inf, y = Inf)

ggplot(obsexp, aes(x = unaltered, y = altered)) +
  specs +
  facet_wrap(pnz~taxon, scales = "free")

# FIGURE 3: Mag shared/unique cat.pair ####

d <- d2.mag.catp

d <- d %>% filter(!is.na(diet.match))
o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.match", "cat.pair", "pnz"), all = T)
avmag <- m %>% select(Taxon_status, diet.match, cat.pair, pnz, subsample = subsample.x, expected = avmag.x, observed = avmag.y) %>% 
  mutate(obsDexp = observed/expected)
avmag$taxon <- word(avmag$Taxon_status, 1, 1, sep = "_")
avmag$status <- word(avmag$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.match+cat.pair+pnz~status, value.var = 'obsDexp', data = avmag, fun.aggregate = median)

anno <- obsexp %>% group_by(taxon, pnz, cat.pair) %>% 
  summarise() %>% ungroup() %>% 
  mutate(label = paste0("  ", LETTERS[1:12]), x = -Inf, y = Inf)

ggplot(obsexp, aes(x = unaltered, y = altered)) + 
  specs+ 
  facet_wrap(cat.pair+taxon~pnz, scales = "free") 


# FIGURE 4: Mag syn/cosmo/restr cosmo.pair ####

d <- d5.mag.catp

d <- d %>% filter(!is.na(diet.match))
o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.match", "cosmo.pair", "pnz"), all = T)
avmag <- m %>% select(Taxon_status, diet.match, cosmo.pair, pnz, subsample = subsample.x, expected = avmag.x, observed = avmag.y) %>% 
  mutate(obsDexp = observed/expected)
avmag$taxon <- word(avmag$Taxon_status, 1, 1, sep = "_")
avmag$status <- word(avmag$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.match+cosmo.pair+pnz~status, value.var = 'obsDexp', data = avmag, fun.aggregate = median)
n <- expand.grid(unique(obsexp$pnz), unique(obsexp$cosmo.pair)) %>% mutate(n = paste(.$Var2, .$Var1, sep = "_")) %>% pull(n) %>% expand.grid(unique(obsexp$taxon)) %>% mutate(n = paste(.$Var2, .$Var1, sep = "_")) %>% pull(n)
b <- besttest(obsexp, split.var = "diet.match.2", taxon, cosmo.pair, pnz) %>% purrr::map(setNames, n)
b <- purrr::map(b, ~.[!sapply(., is_null)])
s <- b %>% purrr::map(~purrr::map_dbl(., ~summary(.) %>% as.data.frame() %>% pull(`%>compVal`) %>% `[`(3)))

anno <- obsexp %>% filter(taxon == "bat") %>% group_by(pnz, cosmo.pair) %>% 
  summarise() %>% ungroup() %>% 
  mutate(label = paste0("  ", LETTERS[1:12]), x = -Inf, y = Inf)

ggplot(obsexp %>% filter(taxon == "bat"), aes(x = unaltered, y = altered)) + 
  specs+ 
  facet_wrap(pnz~cosmo.pair, scales = "free", ncol = 6, nrow = 2) 


anno <- obsexp %>% filter(taxon == "bird") %>% group_by(pnz, cosmo.pair) %>% 
  summarise() %>% ungroup() %>% 
  mutate(label = paste0("  ", LETTERS[1:12]), x = -Inf, y = Inf)

ggplot(obsexp %>% filter(taxon == "bird"), aes(x = unaltered, y = altered)) + 
  specs+ 
  facet_wrap(pnz~cosmo.pair, scales = "free", ncol = 6, nrow = 2) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
##  SUPPLEMENT 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

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


### FIGURE S2: Explanation for plots #####
## TODO
## FIGURE S3: Prop.all ####
d <- d2.prop.all
d <- d %>% filter(!is.na(diet.match))
o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.match"), all = T)
pAgg <- m %>% select(Taxon_status, diet.match, subsample = subsample.x, expected = agg.x, observed = agg.y) %>% 
  mutate(obsDexp = observed/expected)
pAgg$taxon <- word(pAgg$Taxon_status, 1, 1, sep = "_")
pAgg$status <- word(pAgg$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.match~status, value.var = 'obsDexp', data = pAgg)

ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.match, fill = diet.match)) + 
  specs+ 
  facet_grid(taxon~., scales = "fixed") 

## FIGURE S4: Prop shared/unique cat.pair #### 
 d <- d2.prop.catp
 
 d <- d %>% filter(!is.na(diet.match))
 o <- d[d$type == "observed",]
 e <- d[d$type == "expected",]
 m <- merge(e, o, by = c("Taxon_status", "diet.match", "cat.pair"), all = T)
 pAgg <- m %>% select(Taxon_status, diet.match, cat.pair, subsample = subsample.x, expected = agg.x, observed = agg.y) %>% 
   mutate(obsDexp = observed/expected)
 pAgg$taxon <- word(pAgg$Taxon_status, 1, 1, sep = "_")
 pAgg$status <- word(pAgg$Taxon_status, 2, 2, sep = "_")
 
 obsexp <- dcast(subsample+taxon+diet.match+cat.pair~status, value.var = 'obsDexp', data = pAgg)
 
ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.match, fill = diet.match)) + 
  specs+  
  facet_wrap(taxon~cat.pair, scales = "free") 
 
## FIGURE S5: Prop syn/cosmo/rest cosmo.pair #### 
d <- d5.prop.catp

d <- d %>% filter(!is.na(diet.match.2))
o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.match.2", "cosmo.pair"), all = T)
pAgg <- m %>% select(Taxon_status, diet.match.2, cosmo.pair, subsample = subsample.x, expected = agg.x, observed = agg.y) %>% 
  mutate(obsDexp = observed/expected)
pAgg$taxon <- word(pAgg$Taxon_status, 1, 1, sep = "_")
pAgg$status <- word(pAgg$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.match.2+cosmo.pair~status, value.var = 'obsDexp', data = pAgg)

n <- expand.grid(unique(obsexp$cosmo.pair), unique(obsexp$taxon)) %>% mutate(n = paste(.$Var2, .$Var1, sep = "_")) %>% pull(n)
b <- besttest(obsexp, split.var = "diet.match.2", taxon, cosmo.pair) %>% purrr::map(setNames, n)
b <- purrr::map(b, ~.[!sapply(., is_null)])
s <- b %>% purrr::map(~purrr::map_dbl(., ~summary(.) %>% as.matrix() %>% as.data.frame() %>% pull(`%>compVal`) %>% `[`(3)))

ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.match.2, fill = diet.match.2)) + 
  specs+  
  facet_wrap(taxon~cosmo.pair, scales = "free") 


## FIGURE S6: Altered habitat types ########

pbat <- ggplot(pr.bat, aes(x = altered_habitat, y = value*100, fill= altered_habitat)) + geom_boxplot(outlier.size = 0.3) +
  facet_wrap(variable~., nrow = 4, ncol = 1, scales = "free_y") + geom_jitter(width = 0.15, size = 0.3) +
  panel_border(remove = F, colour = "black") + scale_fill_manual(values = c(gg_color_hue(9), "gray50")) +
  theme(axis.text.x = element_blank()) + labs(y = "Percent of species co-occurring", x = "Habitat type")


pbird <- ggplot(pr.bird, aes(x = altered_habitat, y = value*100, fill= altered_habitat)) + geom_boxplot(outlier.size = 0.3) +
  facet_wrap(variable~., nrow = 4, ncol = 2, scales = "free_y") + geom_jitter(width = 0.15, size = 0.3) +
  panel_border(remove = F, colour = "black") + 
  theme(axis.text.x = element_blank()) + labs(y = "Percent of species co-occurring", x = "Habitat type", fill = "Type of alteration") + 
  scale_fill_manual(values = c(gg_color_hue(9), "gray50"), 
                    limits = c("cropland", "disturbed forest", "fragment", "pasture", "plantation", "rural", "secondary forest", "suburban", "urban", "unaltered"),
                    labels = c("cropland", "disturbed forest", "fragment", "pasture", "plantation", "rural", "secondary forest", "suburban", "urban", "unaltered"))

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

#####
#####
#####