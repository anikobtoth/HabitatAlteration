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

# bird and bat theta posterior mean and CI calculations ----

# calculate posterior slope (a) and intercept (b) for each of 4 groups
bird_post <- as.data.frame(bird_winner)
# slopes 
x <- data.frame(iteration = 1:nrow(bird_post), 
                slope_control_altered = bird_post$`b_D_cut`, 
                slope_control_intact = bird_post$`b_D_cut`,
                slope_intraguild_altered = bird_post$`b_D_cut` + bird_post$`b_D_cut:diet.matchsame`, 
                slope_intraguild_intact = bird_post$`b_D_cut` + bird_post$`b_D_cut:diet.matchsame`, 
                # intercepts 
                intercept_control_altered = bird_post$`b_Intercept`,
                intercept_control_intact = bird_post$`b_Intercept` + bird_post$`b_statusunaltered`,
                intercept_intraguild_altered = bird_post$`b_Intercept` + bird_post$`b_diet.matchsame`,
                intercept_intraguild_intact = bird_post$`b_Intercept` + bird_post$`b_diet.matchsame` + bird_post$`b_statusunaltered`) %>%
  pivot_longer(contains("_")) %>% 
  separate(name, into = c("variable", "diet.match", "status"), sep = "_") %>% 
  pivot_wider(names_from = "variable", values_from = "value")

CI_df_bird <-seq(range(bird_data_plot$D)[1], range(bird_data_plot$D)[2], length.out = 35) %>%
  map(~x %>% mutate(D = .x, theta = slope*D+intercept)) %>%
  bind_rows() %>%
  group_by(diet.match, status, D) %>%
  summarise(theta_025 = quantile(theta, 0.025),
            theta_975 = quantile(theta, 0.975),
            theta_mu  = mean(theta)) %>%
  full_join(bird_data_plot %>% group_by(diet.match, status) %>% 
              summarise(d025 = quantile(D, 0.025), d975 = quantile(D, 0.975))) %>%
  mutate(D95q = D > d025 & D < d975)


bat_post <- as.data.frame(bat_winner)  
x <- data.frame(iteration = 1:nrow(bat_post), 
                slope= bat_post$`b_D_cut`, 
                intercept = bat_post$`b_Intercept`)

CI_df_bat <-seq(range(bat_data_plot$D)[1], range(bat_data_plot$D)[2], length.out = 35) %>%
  map(~x %>% mutate(D = .x, theta = slope*D+intercept)) %>%
  bind_rows() %>%
  group_by(D) %>%
  summarise(theta_025 = quantile(theta, 0.025),
            theta_975 = quantile(theta, 0.975),
            theta_mu  = mean(theta)) %>%
  cross_join(bat_data_plot %>% group_by(diet.match, status) %>% 
               summarise(d025 = quantile(D, 0.025), d975 = quantile(D, 0.975))) %>%
  mutate(D95q = D > d025 & D < d975)
# end mean and CI calculations ----

colors <- c("#ADADAD", "#525252", "#2F99DC", "#00487C")

## BAT PLOT FIG 2 ####
tax = "bat"
data_plot <- stan_data_fun(filter(contables, taxon == tax), medn = FALSE)[[2]]
bat_data_plot <- as.data.frame(bat_winner) %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = 8) %>% 
  pivot_wider() %>% cbind(bat_data) %>% full_join(data_plot) %>% 
  mutate(diet.match = recode(diet.match, "different" = "control", "same" = "intraguild"), 
         status = recode(status, "unaltered" = "intact"))

# bat plot
bat_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(status, diet.match))) + 
  geom_jitter(alpha = 0.3, width = 0.1) + 
  geom_ribbon(data = filter(CI_df_bat, !D95q), aes(x = D, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = D, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(slope = -0.12, intercept = 0.41, col = "black", lwd = 0.4, lty = 2) + 
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = D, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(status~diet.match) + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() + theme(legend.position = "none")
######


## Bird plot FIG 3 ####

tax = "bird"
data_plot <- stan_data_fun(filter(contables, taxon == tax), medn = FALSE)[[2]]
bat_data_plot <- as.data.frame(bird_winner) %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = 8) %>% 
  pivot_wider() %>% cbind(bird_data) %>% full_join(data_plot) %>% 
  mutate(diet.match = recode(diet.match, "different" = "control", "same" = "intraguild"), 
         status = recode(status, "unaltered" = "intact"))

# bird plot
linefit <- expand.grid(status = c("intact", "altered"), diet.match = c("intraguild", "control")) %>% 
  mutate(intercept = c(0.24, 0.41, 0.24, 0.41), 
         slope  = c(-0.18, -0.18, -0.02, -0.02))

bird_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(status, diet.match))) + 
  geom_jitter(alpha = 0.3, width = 0.1) + 
  geom_ribbon(data = CI_df_bird, aes(x = D, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bird, D95q), aes(x = D, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bird, D95q), aes(x = D, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(status~diet.match) + labs(col = "group") + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() + theme(legend.position = "none")


######

## SD plot FIG 4 ----
sd <- list(bat_winner, bat_winner_NT, bird_winner, bird_winner_NT) %>% 
  map(~data.frame(summary(.x)$random) %>% select(1:4) %>% 
        setNames(c("Estimate", "SE", "l-95CI", "u-95CI")) %>% rownames_to_column()) %>% 
  setNames(c("bat_full", "bat_NT", "bird_full", "bird_NT")) %>% 
  bind_rows(.id = "model") %>% 
  separate(model, into = c("taxon", "model"), sep = "_") %>% 
  mutate(status = ifelse(grepl("unaltered", rowname), "intact", "altered"), 
         diet.match = ifelse(grepl("same", rowname), "intraguild", "control")) %>%
  rbind(filter(.,taxon == "bird", diet.match == "control") %>% 
          mutate(diet.match = "intraguild")) %>% select(-rowname)

sd %>% ggplot(aes(x = interaction(status, diet.match), y = Estimate, fill = interaction(status, diet.match))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = `l-95CI`, ymax = `u-95CI`), width = 0.1, lwd = 0.5) +
  facet_grid(model~taxon) + 
  theme_light() + labs(fill = "group") +
  scale_fill_manual(values = colors) + 
  coord_flip(ylim = c(0.3, 0.81)) + 
  theme(legend.position = "none") + 
  labs(x = "", y = paste('Standard deviation of', expression(theta)))
# End SD plot ----


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

## Model checks Figs S6-S8 -----
## normality check
p1 <- bird_data_plot %>% ggplot(aes(x = theta_vl)) + geom_histogram() + facet_grid(status~diet.match)
p2 <- bat_data_plot %>% ggplot(aes(x = theta_vl)) + geom_histogram() + facet_grid(status~diet.match)
plot_grid(p2, p1, labels = c("A", "B"), ncol = 2)

## bat bias check
p1 <- bat_data_plot %>% ggplot(aes(x = log(N_sb), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.3) + geom_smooth(method = "loess") + facet_grid(status~diet.match) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = expression(ln(N[set])))

p2 <- bat_data_plot %>% ggplot(aes(x = log(occ1), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess") + facet_grid(status~diet.match) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of rarer species in pair")
plot_grid(p2, p1, labels = c("A", "B"), ncol = 2)

## bird bias check
p1 <- bird_data_plot %>% slice_sample(n = 100000) %>% ggplot(aes(x = log(N_sb), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.3) + geom_smooth(method = "loess") + facet_grid(status~diet.match) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = expression(ln(N[set])))

p2 <- bird_data_plot %>% slice_sample(n = 50000) %>% ggplot(aes(x = log(occ1), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess") + facet_grid(status~diet.match) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of rarer species in pair")
plot_grid(p2, p1, labels = c("A", "B"), ncol = 2)

# end model checks----



