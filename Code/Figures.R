### Effects of phylogenetic distance, niche overlap and habitat ###
### alteration on spatial co-occurrence patterns in Neotropical ###
### bats and birds.  

### By Aniko B. Toth, John Alroy, S Kathleen Lyons, and Andrew P. Allen

# Submitted 15 July 2024 

# Figures script
# All figure scripts will run if Analysis_Script.R is run first;
# Please note warnings about computing requirements 

library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(hexbin)
library(viridis)
library(gganimate)
library(latex2exp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
##  MAIN TEXT 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Figure 1 #########

generate_data <- function(n, slope, intercept, sd) {
  x <- seq(-2, 2, length.out = n)
  y <- intercept + slope * x + rnorm(n, sd = sd)
  data.frame(x, y)
}

# Generate data for each panel
set.seed(42)
data_list <- list(
  A = bind_rows(list(
    low = generate_data(500, -0.6, 1, 0.2),
    high = generate_data(500, -0.4, 0.5, 0.2)
  ),.id = "Diet overlap"),
  B = bind_rows(list(
    low = generate_data(500, 0.02, 0.5, 0.2),
    high = generate_data(500, 0.022, 1, 0.2)
  ),.id = "Diet overlap"),
  C = bind_rows(list(
    low = generate_data(500, -0.50, 0.81, 0.2),
    high = generate_data(500, -0.53, 0.8, 0.2)
  ),.id = "Diet overlap"),
  D =bind_rows(list(
    low = generate_data(500, -0.31, 0.81, 0.2),
    high = generate_data(500, -0.3, 0.8, 0.5)
  ), .id = "Diet overlap")
) %>% bind_rows(.id= "panel")

data_list %>% ggplot(aes(x = x, y = y, col = `Diet overlap`)) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1) + facet_wrap(~panel) +
  geom_text(data = data_list %>% 
              group_by(panel) %>%
              slice(1) %>%
              mutate(label = panel),
            aes(label = label, x = -1.9, y = 2.6), color = "black", fontface = "bold") +
  scale_color_manual(values = c("#0070C0", "gray25")) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), strip.background = element_blank(),
        strip.text.x = element_blank(), text = element_text(family = "Times New Roman")) +
  labs(x = expression("standardized phylogenetic distance ("*italic(D)*")"), y = TeX("$\\theta$"))

## Figure 2 #####
# fnchd lpmf helper functions used
data <- expand.grid(Na = c(2,8), Nb =c(2,8), N = 50, theta = -2:2, Nab = 0:8) %>% 
  filter(Na == Nb) %>% arrange(Na, theta) %>%
  mutate(Nab_theta = exp(fnchypergeo_lpmf(Nab, Na, Nb, N, theta)), 
         label = ifelse(Na == 2, "italic(N)[A]~'='~italic(N)[B]~'='~2", "italic(N)[A]~'='~italic(N)[B]~'='~8")) %>% 
  na.omit()

custom_breaks <- c(10e-12, 10e-8, 10e-5, 10e-2, 1)
ggplot(data, aes(x = Nab, y = Nab_theta, group = theta, col = theta))+ 
  geom_point(size = 2) + geom_line(lwd = 1) + 
  scale_y_log10(
     breaks = custom_breaks) +
  facet_wrap(~label, labeller = label_parsed, scales= "free") +
  labs(y = expression('probability,'~italic(f)*'('~italic(N)[AB]~'|'~italic(theta)~')'), 
       x= expression('number of co-occurrences,'~italic(N)[AB]), 
       col = expression(theta)) + 
  theme(panel.grid.major = element_line(color = "gray80"))+
  theme_light() + scale_color_gradientn(colors = c("black", "#480355", "#8C2981FF", "#DE4968FF", "#FE9F6DFF"))

## Figure 3 ####
### created in Ppt ###

## Figure 4 ####
colors6 <- c("#ADADAD", "#525252", "#C595CE", "#480355", "#2F99DC", "#00487C")

#### BAT theta posterior mean and CI calculations ----
tax = "bat"
bat_data <- stan_data_fun(filter(contables_full, taxon == tax))[[1]]
data_plot <- stan_data_fun(filter(contables_full, taxon == tax))[[2]]
bat_post <- as.data.frame(bat_FULL_winner)  

bat_data_plot <- bat_post %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bat_data) %>% full_join(data_plot) %>% 
  mutate(DietOvlp =  fct_relevel(DietOvlp, "low", "medium", "high"),
         Habitat = recode(Habitat, "unaltered" = "intact"))

CI_df_bat <- bat_post %>% calc_parameters() %>% CIcalc(data_plot) 

linefit <- expand.grid(DietOvlp = c("high", "medium", "low")) %>% 
  mutate(intercept = c(mean(bat_post$b_Intercept), 
                       mean(bat_post$b_Intercept+bat_post$b_DietOvlpmedium),
                       mean(bat_post$b_Intercept+bat_post$b_DietOvlphigh)), 
         slope = c(mean(bat_post$b_PhyloD))) 

#### Bat plot FIG 4 ####
# FIG. 4A hexbin, no habitat dimension
F4A <- bat_data_plot %>% #mutate(DietOvlp = recode(DietOvlp, "different" = "low", "intersecting" = "medium", "intraguild" = "high")) %>% 
  ggplot(aes(x = D, y = theta_vl)) + 
  geom_hex() + facet_grid(.~DietOvlp) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bat, !D95q, Dbin < 0), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "orange", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope),lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group") + 
  labs(col = "group", y = expression(theta), x = expression("standardized phylogenetic distance ("*italic(D)*")")) + 
  scale_color_manual(values = colors) +
  theme_light()
#### BIRD theta posterior mean and CI calculations ----
# calculate posterior slope (a) and intercept (b) for each of 4 groups
tax = "bird"
bird_data <- stan_data_fun(filter(contables_full, taxon == tax))[[1]]
data_plot <- stan_data_fun(filter(contables_full, taxon == tax))[[2]]
bird_post <- as.data.frame(bird_FULL_winner)
CI_df_bird_full <- bird_post %>% calc_parameters() %>% CIcalc(data_plot)

bird_data_plot <- bird_post %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bird_data) %>% full_join(data_plot) %>% 
  mutate(DietOvlp =  fct_relevel(DietOvlp, "low", "medium", "high"))

#### Bird plot FIG 4 ####
linefit_full <- expand.grid(DietOvlp = c("high","medium", "low")) %>%
  mutate(intercept = mean(bird_post$b_Intercept),
         slope  = mean(bird_post$b_PhyloD))

F4B <-  bird_data_plot %>% ggplot(aes(x = jitter(D), y = theta_vl)) + 
  geom_hex() + facet_grid(.~DietOvlp) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bird_full, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.5) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "orange", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit_full, aes(intercept = intercept, slope = slope),lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group", y = expression(theta), x = expression("standardized phylogenetic distance ("*italic(D)*")")) + 
  theme_light() 

plot_grid(F4A + theme(plot.margin = unit(c(1,.2,0,.2), "cm")), 
          F4B + theme(plot.margin = unit(c(1,.2 ,1,.2), "cm")), 
          ncol = 1, align = "v", labels = c("A Bats", "B Birds"))
## Figure 5 SD plot ----
sd <- list(bat_FULL_winner, bird_RP_winner, bird_FULL_winner, bird_RP_winner) %>% 
  map(~data.frame(summary(.x)$random) %>% select(1:4) %>% 
        setNames(c("Estimate", "SE", "l-95CI", "u-95CI")) %>% rownames_to_column()) %>% 
  setNames(c("bat_Full", "bat_Restricted-pool", "bird_Full", "bird_Restricted-pool")) %>% 
  bind_rows(.id = "model") %>% 
  separate(model, into = c("taxon", "model"), sep = "_") %>% 
  mutate(Variance = ifelse(grepl("high", rowname), "high-overlap", 
                           ifelse(grepl("medium", rowname), "medium-overlap", 
                                  ifelse(grepl("low", rowname), "low-overlap", 
                                         ifelse(grepl("altered", rowname), "altered", 
                                                ifelse(grepl("intact", rowname), "intact", ""))))) %>% 
           fct_relevel("low-overlap", "medium-overlap", "high-overlap", "intact", "altered")) 

colors <- c("#595959","#B382BE","#0070C0", "darkgreen", "limegreen")

sd %>% ggplot(aes(x = Variance, y = Estimate, fill = Variance)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = `l-95CI`, ymax = `u-95CI`), width = 0.1, lwd = 0.5) +
  facet_wrap(taxon~model, scales = "free") + 
  theme_light() + labs(fill = "group") +
  scale_fill_manual(values = colors) + 
  theme(legend.position = "none") + 
  labs(x = "", y = expression('standard deviation of ' * theta))
# End MAIN TEXT figures ----

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



### Figure S2: PCoA ####
#ordinations, colored by alteration type
PAnu <- PAn %>% map(~return(.[,colnames(.) %in% unalt_sites])) # separate unaltered
PAna <- PAn %>% map(~return(.[,!colnames(.) %in% unalt_sites])) # and altered sites

altered <- map(PAna, as.matrix) %>% map(cosine) %>% map(as.dist, upper = F) %>% map(~return(1-.)) %>% 
  map(cmdscale) %>% map(data.frame) %>% bind_rows(.id = "taxon")
intact <- map(PAnu, as.matrix) %>% map(cosine) %>% map(as.dist, upper = F) %>% map(~return(1-.)) %>% 
  map(cmdscale) %>% map(data.frame)  %>% bind_rows(.id = "taxon")
ord1 <- bind_rows(altered = altered, intact = intact, .id = "status")  

altered <- map(tobinary(PAna), as.matrix) %>% map(cosine) %>% map(as.dist, upper = F) %>% map(~return(1-.)) %>% 
  map(cmdscale) %>% map(data.frame) %>% bind_rows(.id = "taxon")
intact <- map(tobinary(PAnu), as.matrix) %>% map(cosine) %>% map(as.dist, upper = F) %>% map(~return(1-.)) %>% 
  map(cmdscale) %>% map(data.frame)  %>% bind_rows(.id = "taxon")
ord2 <- bind_rows(altered = altered, intact = intact, .id = "status")  

ord <- list(abundance = ord1, binary = ord2) %>% bind_rows(.id = "datatype")

ordplot <- ord %>% ggplot(aes(x = X1, y = X2, col = status)) + geom_point() +
  labs(x = "Comp1", y = "Comp2", col = "alteration") + 
  theme(axis.text = element_text(size = 12)) + 
  labs(col = "Habitat") + 
  facet_grid(taxon~datatype)+
  scale_colour_manual(values = c("limegreen", "darkgreen")) + 
  theme_light()

richdist <- ggplot(rich, aes(x = richness, colour = status, fill = status)) + 
  geom_density(lwd = 1, alpha = 0.4, adjust = 1) + facet_grid(taxon~metric) +
  scale_colour_manual(values = c("limegreen", "darkgreen")) + 
  scale_fill_manual(values = c("limegreen", "darkgreen")) +
  scale_x_log10() + labs(col = "Habitat", fill = "Habitat") +
  theme_light()

plot_grid(richdist + theme(legend.position = "none"), ordplot, ncol = 2, labels = c("A", "B")) 

######

## Model checks Figs S4-S6 -----
library(cowplot)
## normality check
p1 <- bird_data_plot %>% ggplot(aes(x = theta_vl)) + geom_histogram(binwidth = .25) + facet_grid(status~DietOvlp) + labs(x = expression(theta~"estimate"))
p2 <- bat_data_plot %>% ggplot(aes(x = theta_vl)) + geom_histogram(binwidth = .25) + facet_grid(status~DietOvlp) + labs(x = expression(theta~"estimate"))
plot_grid(p2, p1, labels = c("A", "B"), ncol = 2)

## bat bias check
p1 <- bat_data_plot %>% ggplot(aes(x = log(N_sb), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.3) + geom_smooth(method = "loess") + 
  facet_grid(status~DietOvlp, scales = "free_y") + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = expression(ln(N[set])))

p2 <- bat_data_plot %>% ggplot(aes(x = log(occ1), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess", se=F) + 
  facet_grid(status~DietOvlp, scales = "free") + 
  ylim(c(0, 1.5))+
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of rarer species in pair")

p3 <-  bat_data_plot %>% ggplot(aes(x = log(occ2), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess", se=F) + 
  facet_grid(status~DietOvlp, scales = "free") + 
  ylim(c(0, 1.5))+
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of more common species in pair")

plot_grid(p2, p3, p1, labels = c("A", "B", "C"), ncol = 1)

## bird bias check
p1 <- bird_data_plot %>% slice_sample(n = 100000) %>% ggplot(aes(x = log(N_sb), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.3) + geom_smooth(method = "loess") + facet_grid(status~DietOvlp) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = expression(ln(N[set])))

p2 <- bird_data_plot %>% slice_sample(n = 50000) %>% ggplot(aes(x = jitter(log(occ1), amount = 0.01), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess", se = FALSE) + facet_grid(status~DietOvlp) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of rarer species in pair")

p3 <- bird_data_plot %>% slice_sample(n = 50000) %>% ggplot(aes(x = jitter(log(occ2), amount = 0.01), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess", se = FALSE) + facet_grid(status~DietOvlp) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of more common species in pair")

plot_grid(p2, p3, p1, labels = c("A", "B", "C"), ncol = 1)

## Rand distributions 
plot_grid(
ggplot(rand_out, aes(x = Estimate)) + geom_histogram(bins = 10) + facet_wrap(taxon~rowname, scales = "free"),
ggplot(rand_out_RP, aes(x = Estimate)) + geom_histogram(bins = 10) + facet_wrap(taxon~rowname, scales = "free"), 
ncol = 1, labels = c("A", "B"))

# end model checks----

# model tables for supplement ----
list(`bat FULL` = bat_FULL_winner, 
     `bat RP` = bat_RP_winner, 
     `bird FULL` = bird_FULL_winner, 
     `bird RP` = bird_RP_winner) %>% map(summary) %>% map(summary_effects) %>% 
  bind_rows(.id = "model") %>% as_tibble() %>% separate(model, into = c("taxon", "model")) %>%
  full_join(rand_out %>% group_by(taxon, model, Variable) %>% 
              summarise(`median Rand` = median(Estimate), 
                        `L-95% CI Rand` = quantile(Estimate, 0.025), 
                        `U-95% CI Rand` = quantile(Estimate, 0.975))) %>% 
  select(-Rhat, -Bulk_ESS, -Tail_ESS) %>%
  write.csv("./Results/brms_model_coeffs.csv")

# bat randeff
readRDS("./stan/bat_loo_randeff.rds") %>% 
  make_loo_table(FALSE, re_steps, "./Results/bat_lc_rand_table.csv")

readRDS("./stan/bat_rp_loo_randeff.rds") %>% 
  make_loo_table(FALSE, re_steps, "./Results/bat_rp_lc_rand_table.csv")

# bird randeff
readRDS("./stan/bird_loo_randeff.rds") %>% 
  make_loo_table(FALSE, re_steps, "./Results/bird_lc_rand_table.csv")

readRDS("./stan/bird_rp_loo_randeff.rds") %>% 
  make_loo_table(FALSE, re_steps, "./Results/bird_rp_lc_rand_table.csv")

# bat fixedeff

readRDS("./stan/bat_loo_fixedeff.rds") %>% 
  make_loo_table(TRUE, fx_steps, "./Results/bat_lc_fixed_table.csv")

readRDS("./stan/bat_rp_loo_fixedeff.rds") %>% 
  make_loo_table(TRUE, fx_steps, "./Results/bat_rp_lc_fixed_table.csv")

# bird fixedeff
readRDS("./stan/bird_loo_fixedeff.rds") %>% 
  make_loo_table(TRUE , fx_steps, "./Results/bird_lc_fixed_table.csv")

readRDS("./stan/bird_rp_loo_fixedeff.rds") %>% 
  make_loo_table(TRUE , fx_steps, "./Results/bird_rp_lc_fixed_table.csv")

## End model tables for supplement ----

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
##   Alternative plotting formats, UNPUBLISHED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Bats scatterplot DietOvl by Habitat ####
bat_data_plot %>% ggplot(aeps(x = D, y = theta_vl, col = interaction(Habitat, DietOvlp))) +
  geom_jitter(alpha = 0.3, width = 0.15, size = 0.5) +
  geom_ribbon(data = filter(CI_df_bat, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(Habitat~DietOvlp) +
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") +
  scale_color_manual(values = colors6) +
  theme_light() + theme(legend.position = "none")

#### Bats hexbins DietOvlp by Habitat ####
bat_data_plot %>% ggplot(aes(x = jitter(D), y = theta_vl)) + 
  geom_hex() + facet_grid(.~DietOvlp) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bat, !D95q, Dbin < 0), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "orange", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() 

#### Birds scatterplot DietOvlp by Habitat ####
bird_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(Habitat, DietOvlp))) + 
  geom_jitter(alpha = 0.2, width = 0.15, size = 0.5) + 
  geom_ribbon(data = filter(CI_df_bird_full, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit_full, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(Habitat~DietOvlp) + labs(col = "group") + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors6) +
  theme_light() + theme(legend.position = "none")

#### Birds hexbins DietOvlp by Habitat ####
birdfull <- bird_data_plot %>% ggplot(aes(x = jitter(D), y = theta_vl)) + 
  geom_hex() + facet_grid(Habitat~DietOvlp) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bird_full, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.5) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "orange", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit_full, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  theme_light() 

#######

## Restricted-Pool plots, not in paper #####
#### BAT Restricted-Pool PLOT  ####
tax = "bat"
bat_data <- stan_data_fun(filter(contables, taxon == tax))[[1]]
data_plot <- stan_data_fun(filter(contables, taxon == tax))[[2]]
bat_post <- as.data.frame(bird_RP_winner)  
CI_df_bat <- bat_post %>% calc_parameters() %>% CIcalc(data_plot)

bat_data_plot <- as.data.frame(bird_RP_winner) %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bat_data) %>% full_join(data_plot) %>% 
  mutate(DietOvlp = recode(DietOvlp, "different" = "low", "same" = "high"),
         Habitat = recode(Habitat, "unaltered" = "intact"))

# bat plot
linefit <- expand.grid(Habitat = c("intact", "altered"), DietOvlp = c("high", "medium", "low")) %>% 
  mutate(intercept = c(0.41, 0.36, 0.77, 0.72, 0.42, 0.37), 
         slope = c(-0.11, -0.11, -0.11, -0.11, -0.11, -0.11))  

bat_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(Habitat, DietOvlp))) + 
  geom_jitter(alpha = 0.2, width = 0.15, size = 0.5) + 
  geom_ribbon(data = filter(CI_df_bat, !D95q, Dbin < 0), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(Habitat~DietOvlp) + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() + theme(legend.position = "none")



#### BIRD Restricted-Pool PLOT ####
tax = "bird"
bird_data <- stan_data_fun(filter(contables, taxon == tax))[[1]]
data_plot <- stan_data_fun(filter(contables, taxon == tax))[[2]]
bird_post <- as.data.frame(bird_RP_winner)
CI_df_bird_full <- bird_post %>% calc_parameters() %>% CIcalc(data_plot)

bird_data_plot_nt <- bird_post %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bird_data) %>% full_join(data_plot) %>% 
  mutate(DietOvlp =  fct_relevel(DietOvlp, "low", "medium", "high"))

#Restricted-pool models linefit
linefit_nt <- expand.grid(Habitat = c("intact", "altered"), DietOvlp = c("high","medium", "low")) %>%
  mutate(intercept = c(0.37, 0.31, 0.37, 0.31, 0.37, 0.31),
         slope  = c(-0.09, -0.1, -0.09, -0.1, -0.09, -0.1))
# restricted-pool hexbins
birdnt <- bird_data_plot_nt %>% ggplot(aes(x = D, y = theta_vl)) + 
  geom_hex() + facet_grid(Habitat~DietOvlp) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bird_nt, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bird_nt, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit_nt, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bird_nt, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  theme_light() 
#### Birds overlay of full and restricted-pool#####
bird_data_plot_diff <- anti_join(bird_data_plot_full, bird_data_plot_nt)
bird_data_plot <- list(retained = bird_data_plot_nt, turnover = bird_data_plot_diff) %>% 
  bind_rows(.id = "type")

bird_data_plot %>%
  ggplot(aes(x = theta_vl, y = as.factor(D%/%3), col = status)) + 
  geom_density_ridges(scale = 1, alpha = 0.4, bandwidth = .15, lwd = 1) + 
  facet_grid(type~DietOvlp) +
  scale_color_manual(values = colors) +
  xlim(c(-1.1, 2.2)) +
  labs(y = "binned standardised phylogenetic distance", x = expression(theta))

#combined bird data NT and FULL 
colours2 <- c("#FF4E33", "#000000")
bird_data_plot <- list(FULL = bird_data_plot_full, NT = bird_data_plot_nt) %% bind_rows(.id = model)
bird_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = model)) +
  geom_jitter(alpha = 0.2, width = 0.15, size = 0.5) +
  facet_grid(Habitat~DietOvlp) +
  geom_abline(data = linefit_bird, aes(intercept = intercept, slope = slope, col = model), lty = 2,lwd = 0.4) +
  scale_color_manual(values = colours2) +
  labs(y = expression(theta), x = "standardised phylogenetic distance")
