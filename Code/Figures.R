##    Habitat alteration reduces food competition and local    ###
##      functional diversity in Neotropical bats and birds.    ##### 
##                     Anikó B. Tóth                           ##

# Submitted 1. Apr 2019 
# Figures script
# All figure scripts will run if Analysis_Script.R is run first;
# Please note warnings about data file sizes. 

library(ggplot2)
library(cowplot)

### Maps ####
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


### Figure 1: PCoA ####
#ordinations, colored by alteration type
p1 <- ggplot(ord[[1]], aes(x = X1, y = X2, col = status)) + geom_point() + labs(x = "Comp1", y = "Comp2", col = "alteration") + theme(axis.text = element_text(size = 12)) + labs(col = "status") + 
  scale_colour_manual(values = c("red", "dodgerblue")) + panel_border(remove = FALSE, col = "black")
p2 <- ggplot(ord[[2]], aes(x = X1, y = X2, col = status)) + geom_point() + labs(x = "Comp1", y = "Comp2", col = "alteration") + theme(axis.text = element_text(size = 12)) + labs(col = "status") + 
  scale_colour_manual(values = c("red", "dodgerblue")) + panel_border(remove = FALSE, col = "black")
p <- plot_grid(p1+theme(legend.position = "none"), 
               p2+theme(legend.position = "none", axis.title.y = element_blank()), ncol = 2, labels = c("A", "B"), 
               align = "lr")
plot_grid(p, get_legend(p1), rel_widths = c(1, .2))



#### Figures 2-5: scatterplots #####

specs <- list(scale_fill_discrete(h = c(250, 130), c = c(20, 100), l = c(40, 50)),  
              theme(legend.position = "none", axis.text = element_text(size = 12)))

#### Diet overlap [PRESENTED IN MAIN TEXT] ####
t <- "bird"  # choose 'bat' (Figs. 2 and 4) or 'bird' (Figs. 3 and 5)
colors3 <- c("#34008C", "#E8004D", "#E8AC00")

  # FIGS 2 and 4 proportion of aggregations and segregations ####
  d <- data.frame(d2.prop.all)
  d <- d[d$taxon == t,]
  d2 <- dcast(subsample+taxon+diet.match~status, value.var = 'agg', data = d)
  
  ggplot(d2, aes(x = Unaltered, y = Altered, col = diet.match)) + geom_point() + 
    panel_border(remove = F, col = "black") + 
    geom_abline(mapping = NULL, data = NULL, slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
    scale_colour_manual(values = colors3) + labs(col = "Diet Overlap")
  
  # FIGS 3 and 5: magnitude of aggregations and segregations.####
  
  d <- d2.mag.cat
  d <- d[d$taxon == t,]
  d$`posnegzero(Z.Score)` <- factor(d$`posnegzero(Z.Score)`, levels = c("Segregation", "Aggregation"))
  
  d2 <- dcast(subsample+taxon+diet.match+`posnegzero(Z.Score)`+cat.group~status, value.var = 'avmag', data = d)
  
  p2 <- ggplot(d2, aes(x = Unaltered, y = Altered, col = diet.match)) + geom_point(size = 0.8) + 
    facet_wrap(cat.group~`posnegzero(Z.Score)`, scales = "free") + panel_border(remove = F, col = "gray30") + 
    geom_abline(mapping = NULL, data = NULL, slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray35") +
    scale_colour_manual(values = colors3)+ labs(col = "Diet Overlap")
  
  d <- d2.mag.all
  d <- d[d$taxon == t,]
  d$`posnegzero(Z.Score)` <- factor(d$`posnegzero(Z.Score)`, levels = c("Segregation", "Aggregation"))
  
  d2 <- dcast(subsample+taxon+diet.match+`posnegzero(Z.Score)`~status, value.var = 'avmag', data = d)
  
  p1 <- ggplot(d2, aes(x = Unaltered, y = Altered, col = diet.match)) + geom_point(size = 0.8) + 
    facet_wrap(`posnegzero(Z.Score)`~., scales = "free") + panel_border(remove = F, col = "gray30") + 
    geom_abline(mapping = NULL, data = NULL, slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray35") +
    scale_colour_manual(values = colors3)+ labs(col = "Diet Overlap")
  
  plot_grid(p1 + theme(legend.position = "none", axis.title.x = element_blank()) ,
            p2 + theme(legend.position = "none"), rel_heights = c(1, 1.85), ncol = 1)
  

### Guild by guild [PRESENTED IN SUPPLEMENT] ####

t <- "bird"  # choose 'bat' or 'bird'

# Proportion aggregated
d <- d3.prop.all
d <- d[d$taxon == t,]
  # keep guilds with 15 pairs or more per subsample on average (this equates to roughly 15 species in the guild, showing up in both site types)
keep <- d %>% group_by(diet.pair, status) %>% summarise(avg = mean(count), mean(agg)) %>% filter(avg>15) 
keep <- names(table(keep$diet.pair)[table(keep$diet.pair) >1])
d <- d[d$diet.pair %in% keep,]
d2 <- dcast(subsample+taxon+diet.pair+diet.match~status, value.var = 'agg', data = d)

ggplot(d2, aes(x = Unaltered, y = Altered, col = diet.pair)) + geom_point() + 
  panel_border(remove = F) + 
  geom_abline(mapping = NULL, data = NULL, slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
  scale_colour_hue(h = c(270, 20, 330, 60, 300, 110, 216), c = 100, l = c(50, 60, 90, 50, 80, 85, 70))

# magnitude of aggregations and segregations.
d <- d3.mag.cat
d <- d[d$taxon == t & d$diet.pair %in% keep,]
d2 <- dcast(subsample+taxon+diet.pair+diet.match+`posnegzero(Z.Score)`+cat.group~status, value.var = 'avmag', data = d)

ggplot(d2, aes(x = Unaltered, y = Altered, col = diet.pair)) + geom_point(size = 0.5) + 
  facet_wrap(`posnegzero(Z.Score)`~cat.group, scales = "free") + panel_border(remove = F) + 
  geom_abline(mapping = NULL, data = NULL, slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
  scale_colour_hue(h = c(270, 20, 330, 60, 300, 110, 216), c = 100, l = c(50, 60, 90, 50, 80, 85, 70))


# Table S1
s <- spp %>% filter(rownames(.) %in% unlist(map(tables, map, rownames))) %>% group_by(life.form, guild) %>% dplyr::summarise(count = length(guild)) %>% na.omit()


## Altered habitat types figure ########

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
#### SUPPLEMENTARY Boxplots (Figs 2-5 presented as boxplots instead of scatter plots) ######

load("~/Desktop/MacQuarie_PhD/Thesis/Chapter_3/R_Files/Results/d2_4xsummarytables.RData")
specs <- list(scale_fill_manual(values = c("red", "dodgerblue")),  
              theme(legend.position = "none", axis.text = element_text(size = 12)))

#####
d <- data.frame(d2.prop.all)
t <- "bat" # choose "bird" (Figs 3 and 5) or "bat" (Figs. 2 and 4)
d$all <- "All"

p1 <- ggplot(d[d$taxon == t,], aes(y = agg, x = status, fill = status)) + geom_boxplot(notch = T) + 
  facet_grid(diet.match~all, scales = "free") + specs + coord_flip() + labs(y = "proportion aggregated") + 
  panel_border(remove = FALSE)

d <- data.frame(d2.prop.cat)
d$cat.group <- factor(d$cat.group, levels = c("Unique", "Shared"))
sp <- d %>% group_by(taxon, status, diet.match, cat.group) %>% summarise(N = mean(count)) %>% filter(taxon == t, status == "Unaltered")
sa <- d %>% group_by(taxon, status, diet.match, cat.group) %>% summarise(N = mean(count)) %>% filter(taxon == t, status == "Altered")
p2 <- ggplot(d[d$taxon == t,], aes(y = agg, x = status, fill = status)) + 
  geom_boxplot(notch = T, width = 0.5, outlier.size = 0.4) + facet_grid(diet.match~cat.group, scales = "free") + 
  specs + coord_flip() + labs(y = "proportion aggregated")  +
  annotate("text", label = expression(paste(bar(N), " =")), size = 3, x = Inf, y = Inf, hjust = 2.7, vjust = 1.2) +
  annotate("text", label = round(sp$N), size = 3, x = Inf, y = Inf, hjust = 1, vjust = 1.4 ) + 
  annotate("text", label = expression(paste(bar(N), " =")), size = 3, x = -Inf, y = Inf, hjust = 2.7, vjust = -0.8) +
  annotate("text", label = round(sa$N), size = 3, x = -Inf, y = Inf, hjust = 1, vjust = -1) + 
  panel_border(remove = FALSE)

plot_grid(p1, p2+theme(axis.text.y = element_blank(), axis.title.y = element_blank()), rel_widths = c(2,3),  labels = c("A", "B"), label_x = c(0.94, 0.96))

###
d <- data.frame(d2.mag.all)
d$posnegzero.Z.Score. <- factor(d$posnegzero.Z.Score., levels = c("Segregation", "Aggregation"))
p1 <- ggplot(d[d$taxon == t,], aes(y = avmag, x = status, fill = status)) + geom_boxplot(notch = T) + 
  facet_grid(diet.match~posnegzero.Z.Score., scales = "free_x") + 
  ggtitle("All") + coord_flip() + specs + panel_border(remove = FALSE) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title.x = element_blank()) 


d <- data.frame(d2.mag.cat)
d$posnegzero.Z.Score. <- factor(d$posnegzero.Z.Score., levels = c("Segregation", "Aggregation"))
sp <- d %>% group_by(taxon, status, diet.match, posnegzero.Z.Score., cat.group) %>% summarise(N = mean(count)) %>% filter(taxon == t, status == "Unaltered")
sa <- d %>% group_by(taxon, status, diet.match, posnegzero.Z.Score., cat.group) %>% summarise(N = mean(count)) %>% filter(taxon == t, status == "Altered")

p2 <- ggplot(d[d$taxon == t,], aes(y = avmag, x = status, fill = status)) + 
  geom_boxplot(notch = T, width = 0.5, outlier.size = 0.4) + facet_grid(diet.match~posnegzero.Z.Score.+cat.group, scales = "free") + ggtitle("Components") + specs + coord_flip() + labs(y = "Mean strength") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + panel_border(remove = FALSE) +
  annotate("text", label = expression(paste(bar(N), " =")), size = 3, x = Inf, y = Inf, hjust = 2.7, vjust = 1.2) +
  annotate("text", label = round(sp$N), size = 3, x = Inf, y = Inf, hjust = 1, vjust = 1.4 ) + 
  annotate("text", label = expression(paste(bar(N), " =")), size = 3, x = -Inf, y = Inf, hjust = 2.7, vjust = -0.8) +
  annotate("text", label = round(sa$N), size = 3, x = -Inf, y = Inf, hjust = 1, vjust = -1)

plot_grid(p1, p2, ncol = 1, rel_heights = c(1,1.5), labels = c("A", "B"))  

##
#




### New scatterplot figures: obs/exp with unalt plotted against alt #####
### DIET OVERLAP ####
#ggplot(d2.prop.all %>% filter(!is.na(diet.match)), aes(y = agg, fill = type, x = diet.match)) + geom_boxplot(notch = T) + facet_wrap(.~Taxon_status, scales = "free") 

colors3 <- c("#34008C", "#E8004D", "#E8AC00")
## Proportion of aggregations
ggplot(d2.prop.cat %>% filter(!is.na(diet.match)), aes(y = agg, fill = type, x = diet.match)) + geom_boxplot(notch = T) + facet_grid(cat.group~Taxon_status, scales = "free") 
## Cat.group ####
d <- d2.prop.cat

d <- d %>% filter(!is.na(diet.match))
o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.match", "cat.group"), all = T)
pAgg <- m %>% select(Taxon_status, diet.match, cat.group, subsample = subsample.x, expected = agg.x, observed = agg.y) %>% 
  mutate(obsDexp = observed/expected)
pAgg$taxon <- word(pAgg$Taxon_status, 1, 1, sep = "_")
pAgg$status <- word(pAgg$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.match+cat.group~status, value.var = 'obsDexp', data = pAgg)

 ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.match, fill = diet.match)) + 
  geom_abline(mapping = NULL, data = NULL, slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
  geom_hline(yintercept = 1, col = "gray") +
  geom_vline(xintercept = 1, col = "gray") +
  #geom_point(size = 0.5) + 
  stat_ellipse(geom= "polygon", alpha = 0.3) +
  facet_grid(taxon~cat.group, scales = "fixed") + 
  scale_colour_manual(values = colors3) +
  scale_fill_manual(values = colors3)  

## cat.pair #### 
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
   geom_abline(mapping = NULL, data = NULL, slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
   geom_hline(yintercept = 1, col = "gray") +
   geom_vline(xintercept = 1, col = "gray") +
   #geom_point(size = 0.5) + 
   stat_ellipse(geom= "polygon", alpha = 0.3) +
   facet_grid(taxon~cat.pair, scales = "free_y") + 
   scale_colour_manual(values = colors3) +
   scale_fill_manual(values = colors3)  
 
# strength of associations
#ggplot(d2.mag.all %>% filter(!is.na(diet.match)), aes(y = abs(avmag), fill = type, x = diet.match)) + geom_boxplot(notch = T) + facet_grid(`posnegzero(Z.Score)`~Taxon_status, scales = "free") 
ggplot(d2.mag.cat %>% filter(!is.na(diet.match)), aes(y = abs(avmag), fill = type, x = diet.match)) + geom_boxplot(notch = T) + facet_grid(cat.group+`posnegzero(Z.Score)`~Taxon_status, scales = "free") 

# cat.group ####
d <- d2.mag.cat

d <- d %>% filter(!is.na(diet.match))
o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.match", "cat.group", "posnegzero(Z.Score)"), all = T)
avmag <- m %>% select(Taxon_status, diet.match, cat.group, `posnegzero(Z.Score)`, subsample = subsample.x, expected = avmag.x, observed = avmag.y) %>% 
  mutate(obsDexp = observed/expected)
avmag$taxon <- word(avmag$Taxon_status, 1, 1, sep = "_")
avmag$status <- word(avmag$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.match+cat.group+`posnegzero(Z.Score)`~status, value.var = 'obsDexp', data = avmag, fun.aggregate = median)

ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.match, fill = diet.match)) + 
  geom_abline(mapping = NULL, data = NULL, slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
  geom_hline(yintercept = 1, col = "gray") +
  geom_vline(xintercept = 1, col = "gray") +
  #geom_point(size = 0.5) + 
  stat_ellipse(geom= "polygon", alpha = 0.1) +
  facet_grid(cat.group+taxon~`posnegzero(Z.Score)`, scales = "free") + 
  scale_colour_manual(values = colors3) +
  scale_fill_manual(values = colors3)  

# cat.pair ####

d <- d2.mag.catp

d <- d %>% filter(!is.na(diet.match))
o <- d[d$type == "observed",]
e <- d[d$type == "expected",]
m <- merge(e, o, by = c("Taxon_status", "diet.match", "cat.pair", "posnegzero(Z.Score)"), all = T)
avmag <- m %>% select(Taxon_status, diet.match, cat.pair, `posnegzero(Z.Score)`, subsample = subsample.x, expected = avmag.x, observed = avmag.y) %>% 
  mutate(obsDexp = observed/expected)
avmag$taxon <- word(avmag$Taxon_status, 1, 1, sep = "_")
avmag$status <- word(avmag$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.match+cat.pair+`posnegzero(Z.Score)`~status, value.var = 'obsDexp', data = avmag, fun.aggregate = median)

ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.match, fill = diet.match)) + 
  geom_abline(mapping = NULL, data = NULL, slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
  geom_hline(yintercept = 1, col = "gray") +
  geom_vline(xintercept = 1, col = "gray") +
  #geom_point(size = 0.5) + 
  stat_ellipse(geom= "polygon", alpha = 0.1) +
  facet_wrap(cat.pair+taxon~`posnegzero(Z.Score)`, scales = "free") + 
  scale_colour_manual(values = colors3) +
  scale_fill_manual(values = colors3)  

### GUILD-BY-GUILD #####
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
m <- merge(e, o, by = c("Taxon_status", "diet.pair", "cat.group", "posnegzero(Z.Score)"), all = T)

avmag <- m %>% select(Taxon_status, diet.pair, cat.group, `posnegzero(Z.Score)`, subsample = subsample.x, expected = avmag.x, observed = avmag.y) %>% 
  mutate(obsDexp = observed/expected)
avmag$taxon <- word(avmag$Taxon_status, 1, 1, sep = "_")
avmag$status <- word(avmag$Taxon_status, 2, 2, sep = "_")

obsexp <- dcast(subsample+taxon+diet.pair+cat.group+`posnegzero(Z.Score)`~status, value.var = 'obsDexp', data = avmag, fun.aggregate = median)
obsexp <- obsexp %>% filter(is.finite(altered), is.finite(unaltered)) %>% na.omit()

ggplot(obsexp, aes(x = unaltered, y = altered, col = diet.pair, fill = diet.pair)) + 
  geom_abline(slope = 1, intercept = 0, show.legend = NA, lty = 2, col = "gray30") +
  geom_hline(yintercept = 1, col = "gray") +
  geom_vline(xintercept = 1, col = "gray") +
  geom_point(size = 0.5) + 
  stat_ellipse(geom= "polygon", alpha = 0.1) +
  facet_grid(taxon+`posnegzero(Z.Score)`~cat.group, scales = "free") + 
  scale_colour_hue(h = c(270, 20, 330, 60, 300, 110, 216), c = 100, l = c(50, 60, 90, 50, 80, 85, 70)) +
  scale_fill_hue(h = c(270, 20, 330, 60, 300, 110, 216), c = 100, l = c(50, 60, 90, 50, 80, 85, 70))

