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
tax = "bat"
data_plot <- stan_data_fun(filter(contables, taxon == tax), medn = FALSE)[[2]]
bat_data_plot <- as.data.frame(bat_winner) %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = 8) %>% 
  pivot_wider() %>% cbind(bat_data) %>% full_join(bat_data_plot) %>% 
  mutate(diet.match = recode(diet.match, "different" = "control", "same" = "intraguild"), 
         status = recode(status, "unaltered" = "intact"))

colors <- c("#ADADAD", "#525252", "#2F99DC", "#00487C")
# bat plot
bat_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(status, diet.match))) + 
  geom_jitter(alpha = 0.3, width = 0.1) + 
  #geom_hex(aes(fill = stat(count))) +
  geom_ribbon(data = filter(CI_df_bat, !D95q), aes(x = D, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = D, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(slope = -0.12, intercept = 0.41, col = "black", lwd = 0.4, lty = 2) + 
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = D, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(status~diet.match) + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() + theme(legend.position = "none")

## Bird plot ####

tax = "bird"

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

### Figure S5: PCoA ####
#ordinations, colored by alteration type
altered <- map(PAna, as.matrix) %>% map(cosine) %>% map(as.dist, upper = F) %>% map(~return(1-.)) %>% 
  map(cmdscale) %>% map(data.frame) %>% bind_rows(.id = "taxon")
intact <- map(PAnu, as.matrix) %>% map(cosine) %>% map(as.dist, upper = F) %>% map(~return(1-.)) %>% 
  map(cmdscale) %>% map(data.frame)  %>% bind_rows(.id = "taxon")
ord <- bind_rows(altered = altered, intact = intact, .id = "status")  

ord <- ord %>% ggplot(aes(x = X1, y = X2, col = status)) + geom_point() +
  labs(x = "Comp1", y = "Comp2", col = "alteration") + 
  theme(axis.text = element_text(size = 12)) + labs(col = "status") + 
  facet_grid(~taxon)+
  scale_colour_manual(values = c("limegreen", "darkgreen"))




######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## OLD STUFF ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## specs 
colors <- c("#FFCB5C", "#FF791F", "#2F99DC", "#084887")
specs <- list(geom_density(lwd = .8, alpha = .3),
              facet_wrap(parameter~., scales = "free_x"), 
              scale_color_manual(values = colors), 
              scale_fill_manual(values = colors))

# FIGURE 2: bat results#
load("./Results/Omega/interaction/out_jagsfit_bat_medianTRUE_100reps_interaction_baggingFALSE.RData")
a1 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs 
load("./Results/Omega/interaction/out_jagsfit_bat_medianFALSE_100reps_interaction_baggingFALSE.RData")
a2 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")
load("./Results/Omega/interaction/out_jagsfit_bat_medianTRUE_100reps_interaction_baggingTRUE.RData")
a3 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")
load("./Results/Omega/interaction/out_jagsfit_bat_medianFALSE_100reps_interaction_baggingTRUE.RData")
a4 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")

legend <- get_legend(a1+ theme(legend.box.margin = margin(0, 0, 0, 1)))

p1 <- plot_grid(a1 + theme(legend.position = "none"), a2, a3, a4, align = "hv", nrow = 2)
plot_grid(p1, legend, nrow = 1, rel_widths = c(2, 0.5))

# FIGURE 3: bird results #
load("./Results/Omega/interaction/out_jagsfit_bird_medianTRUE_100reps_interaction_baggingFALSE.RData")
r1 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs
load("./Results/Omega/interaction/out_jagsfit_bird_medianFALSE_100reps_interaction_baggingFALSE.RData")
r2 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")
load("./Results/Omega/interaction/out_jagsfit_bird_medianTRUE_100reps_interaction_baggingTRUE.RData")
r3 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")
load("./Results/Omega/interaction/out_jagsfit_bird_medianFALSE_100reps_interaction_baggingTRUE.RData")
r4 <- ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  specs + theme(legend.position = "none")

legend <- get_legend(r1+ theme(legend.box.margin = margin(0, 0, 0, 7)))

p1 <- plot_grid(r1 + theme(legend.position = "none"), r2, r3, r4, align = "hv", nrow = 2)
plot_grid(p1, legend, nrow = 1, rel_widths = c(2, 0.5))



# Figs 2 and 3 with extra labels: #
load("./Results/Omega/interaction/out_jagsfit_bat_medianFALSE_100reps_interaction_baggingFALSE.RData")
bat_FF <- out
load("./Results/Omega/interaction/out_jagsfit_bat_medianFALSE_100reps_interaction_baggingTRUE.RData")
bat_FT <- out
load("./Results/Omega/interaction/out_jagsfit_bat_medianTRUE_100reps_interaction_baggingFALSE.RData")
bat_TF <- out
load("./Results/Omega/interaction/out_jagsfit_bat_medianTRUE_100reps_interaction_baggingTRUE.RData")
bat_TT <- out

resbat <- bind_rows(list("common_only-bagged" = bat_TT, "common_only-not_bagged" = bat_TF, "all-bagged" = bat_FT, "all-not_bagged" = bat_FF), .id = "model") %>% 
  separate(model, into = c("species", "bagging"), sep = "-") %>% 
  mutate(species = factor(species, levels = c("common_only", "all")), 
         bagging = factor(bagging, levels = c("not_bagged", "bagged")))
ggplot(resbat %>% filter(parameter != "deviance"), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  geom_density(alpha = 0.25, lwd = 1) + 
  facet_wrap(bagging~species+parameter, scales = "free", nrow = 2, ncol = 4) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + ggtitle("Bats")


load("./Results/Omega/interaction/out_jagsfit_bird_medianFALSE_100reps_interaction_baggingFALSE.RData")
bird_FF <- out
load("./Results/Omega/interaction/out_jagsfit_bird_medianFALSE_100reps_interaction_baggingTRUE.RData")
bird_FT <- out
load("./Results/Omega/interaction/out_jagsfit_bird_medianTRUE_100reps_interaction_baggingFALSE.RData")
bird_TF <- out
load("./Results/Omega/interaction/out_jagsfit_bird_medianTRUE_100reps_interaction_baggingTRUE.RData")
bird_TT <- out

resbird <- bind_rows(list("common_only-bagged" = bird_TT, "common_only-not_bagged" = bird_TF, "all-bagged" = bird_FT, "all-not_bagged" = bird_FF), .id = "model") %>% 
  separate(model, into = c("species", "bagging"), sep = "-") %>% 
  mutate(species = factor(species, levels = c("common_only", "all")), 
         bagging = factor(bagging, levels = c("not_bagged", "bagged")))

ggplot(resbird %>% filter(parameter != "deviance"), aes(x = value, col = status_dietmatch, fill = status_dietmatch)) + 
  geom_density(alpha = 0.25, lwd = 1) + 
  facet_wrap(bagging~species+parameter, scales = "free", nrow = 2, ncol = 4) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + ggtitle("Birds")


## Figures for CPEG talk 
out %>% filter(parameter == "omega") %>% 
  separate(status_dietmatch, c("status", "diet")) %>% 
  mutate(parameter = recode(parameter, omega = "Co-occurrence"), 
         status = recode(status, Altered = "Disturbed", Unaltered = "Intact")) %>% 
  mutate(status = factor(status, levels = c("Intact", "Disturbed")), 
         diet = factor(diet, levels = c("Same", "Different"))) %>% 
  ggplot(aes(x = value, col = diet, fill = diet)) + geom_density(alpha = 0.4, lwd = 1) + 
  facet_grid(status~parameter, scales= "free_x")



