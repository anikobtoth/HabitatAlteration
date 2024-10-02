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

colors <- c("#ADADAD", "#525252", "#2F99DC", "#00487C", "#C595CE", "#480355")

# calculate posterior slopes and intercepts
calc_parameters <- function(post){
  params <- list()
  c <- 0
  coeffs <- grep(names(post), pattern = "^b_", value = TRUE)
  for(i in c("control", "intersecting", "same")){
    for(j in c("altered", "unaltered")){
      c <- c+1
      slope <- data.frame(post$`b_PhyloD`,
                          post[,grep(coeffs, pattern = paste0("PhyloD:DietGroup", i), value = TRUE)],
                          post[,grep(coeffs, pattern = paste0("PhyloD:Habitat", j), value = TRUE)]) %>%
        rowSums()
      
      intcoeffs <- coeffs[!grepl(coeffs, pattern = "PhyloD")] # remove slope coefficients
      intercept <- data.frame(post$`b_Intercept`, 
                              post[,grep(intcoeffs, pattern = paste0("DietGroup", i, "|Habitat", j), value = TRUE)]) %>%
        rowSums()
      
      params[[c]] <- tibble(iteration = 1:nrow(post), DietGroup = i, Habitat = j, slope, intercept)
    }
  }
  params <- bind_rows(params) %>% mutate(Habitat = recode(Habitat, "unaltered" = "intact"), 
                                         DietGroup = recode(DietGroup, "same" = "intraguild"))
  return(params)
}

# Calculate credible intervals
CIcalc <- function(params, plotdata){
    seq(range(plotdata$D)[1], range(plotdata$D)[2], length.out = 100) %>%
    map(~params %>% mutate(Dbin = .x, theta = slope*Dbin+intercept)) %>%
    bind_rows() %>%
    group_by(DietGroup, Habitat, Dbin) %>%
    summarise(theta_025 = quantile(theta, 0.025),
              theta_975 = quantile(theta, 0.975),
              theta_mu  = mean(theta)) %>%
    full_join(plotdata %>% group_by(DietGroup, Habitat) %>% 
                summarise(d025 = quantile(D, 0.025), d975 = quantile(D, 0.975))) %>%
    mutate(D95q = Dbin > d025 & Dbin < d975)
}


# BAT theta posterior mean and CI calculations ----
tax = "bat"
data_plot <- stan_data_fun(filter(contables, taxon == tax))[[2]]
bat_post <- as.data.frame(bat_winner)  

CI_df_bat <- bat_post %>% calc_parameters() %>% CIcalc(data_plot)
# CI_df_bat <- data.frame(iteration = 1:nrow(bat_post), 
#                 #slope
#                 slope = bat_post$`b_PhyloD`, 
#                 # intercepts
#                 intercept_control_intact      =  bat_post$`b_Intercept`,
#                 intercept_control_altered     =  bat_post$`b_Intercept`,
#                 intercept_intersecting_intact =  bat_post$`b_Intercept` + bat_post$`b_DietGroupintersecting`,
#                 intercept_intersecting_altered = bat_post$`b_Intercept` + bat_post$`b_DietGroupintersecting`,
#                 intercept_intraguild_intact   = bat_post$`b_Intercept` + bat_post$`b_DietGroupsame`,
#                 intercept_intraguild_altered  = bat_post$`b_Intercept` + bat_post$`b_DietGroupsame`) %>% 
#   pivot_longer(contains("_")) %>% 
#   separate(name, into = c("variable", "DietGroup", "Habitat"), sep = "_") %>% 
#   pivot_wider(names_from = "variable", values_from = "value") %>%
#   CIcalc(data_plot)

# end mean and CI calculations ----

## BAT PLOT FIG 2 ####
bat_data_plot <- as.data.frame(bat_winner) %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bat_data) %>% full_join(data_plot) %>% 
  mutate(Habitat = recode(Habitat, "unaltered" = "intact"))

# bat plot
linefit <- expand.grid(Habitat = c("intact", "altered"), DietGroup = c("intraguild", "intersecting", "control")) %>% 
  mutate(intercept = c(0.44, 0.44, 0.73, 0.73, 0.40, 0.40), 
         slope = c(-0.09, -0.09, -0.09, -0.09, -0.09, -0.09))  

bat_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(Habitat, DietGroup))) + 
  geom_jitter(alpha = 0.3, width = 0.1, size = 0.5) + 
  geom_ribbon(data = filter(CI_df_bat, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(Habitat~DietGroup) + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() + theme(legend.position = "none")
######

# BIRD theta posterior mean and CI calculations ----
# calculate posterior slope (a) and intercept (b) for each of 4 groups
tax = "bird"
data_plot <- stan_data_fun(filter(contables, taxon == tax))[[2]]
bird_post <- as.data.frame(bird_winner)
CI_df_bird <- bird_post %>% calc_parameters() %>% CIcalc(data_plot)
# CI_df_bird <- data.frame(iteration = 1:nrow(bird_post), 
#                      #slope = bird_post$`b_PhyloD`,
#                      slope_control_altered = bird_post$`b_PhyloD`, 
#                      slope_control_intact = bird_post$`b_PhyloD`,
#                      slope_intersecting_altered = bird_post$`b_PhyloD`+ bird_post$`b_PhyloD:DietGroupintersecting`, 
#                      slope_intersecting_intact = bird_post$`b_PhyloD`+ bird_post$`b_PhyloD:DietGroupintersecting`,
#                      slope_intraguild_altered = bird_post$`b_PhyloD` + bird_post$`b_PhyloD:DietGroupsame`, 
#                      slope_intraguild_intact = bird_post$`b_PhyloD` + bird_post$`b_PhyloD:DietGroupsame`, 
#                      # intercepts 
#                      intercept_control_altered = bird_post$`b_Intercept`,
#                      intercept_control_intact = bird_post$`b_Intercept`+ bird_post$`b_Habitatunaltered`,
#                      intercept_intersecting_altered = bird_post$`b_Intercept` + bird_post$`b_DietGroupintersecting`,
#                      intercept_intersecting_intact = bird_post$`b_Intercept`+ bird_post$`b_DietGroupintersecting`+ bird_post$`b_Habitatunaltered`, 
#                      intercept_intraguild_altered = bird_post$`b_Intercept`+ bird_post$`b_DietGroupsame`, 
#                      intercept_intraguild_intact = bird_post$`b_Intercept`+ bird_post$`b_DietGroupsame` + bird_post$`b_Habitatunaltered`
# ) %>%
#   pivot_longer(contains("_")) %>% 
#   separate(name, into = c("variable", "DietGroup", "Habitat"), sep = "_") %>% 
#   pivot_wider(names_from = "variable", values_from = "value") %>%
#   CIcalc(data_plot)

## Bird plot FIG 3 ####

bird_data_plot <- bird_post %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bird_data) %>% full_join(data_plot) %>% 
  mutate(DietGroup = recode(DietGroup, "different" = "control", "same" = "intraguild"), 
         Habitat = recode(Habitat, "unaltered" = "intact"))

# bird plot
linefit <- expand.grid(Habitat = c("intact", "altered"), DietGroup = c("intraguild","intersecting", "control")) %>% 
  mutate(intercept = c(0.44, 0.38, 0.36, 0.30, 0.38, 0.32), 
         slope  = c(-0.11, -0.11, -0.1, -0.1, -0.04, -0.04))

bird_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(Habitat, DietGroup))) + 
  geom_jitter(alpha = 0.2, width = 0.1) + 
  geom_ribbon(data = filter(CI_df_bird, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bird, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bird, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(Habitat~DietGroup) + labs(col = "group") + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() + theme(legend.position = "none")



#######
##### No-Turnover plots #####

## BAT NO TURNOVER PLOT  ####
tax = "bat"
data_plot <- stan_data_fun(filter(contables, taxon == tax))[[2]]
bat_post <- as.data.frame(bat_nt_winner)  
params <- data.frame(iteration = 1:nrow(bat_post), 
                     #slope
                     slope= bat_post$`b_PhyloD`, 
                     # intercepts
                     intercept_control_intact      =  bat_post$`b_Intercept`,
                     intercept_control_altered     =  bat_post$`b_Intercept`,
                     intercept_intersecting_intact =  bat_post$`b_Intercept` + bat_post$`b_DietGroupintersecting`,
                     intercept_intersecting_altered = bat_post$`b_Intercept` + bat_post$`b_DietGroupintersecting`,
                     intercept_intraguild_intact   = bat_post$`b_Intercept` + bat_post$`b_DietGroupintraguild`,
                     intercept_intraguild_altered  = bat_post$`b_Intercept` + bat_post$`b_DietGroupintraguild`) %>% 
  pivot_longer(contains("_")) %>% 
  separate(name, into = c("variable", "DietGroup", "Habitat"), sep = "_") %>% 
  pivot_wider(names_from = "variable", values_from = "value") %>%
  CIcalc(data_plot)

bat_data_plot <- as.data.frame(bat_nt_winner) %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bat_data) %>% full_join(data_plot) %>% 
  mutate(Habitat = recode(Habitat, "unaltered" = "intact"))

# bat plot
linefit <- expand.grid(Habitat = c("intact", "altered"), DietGroup = c("intraguild", "intersecting", "control")) %>% 
  mutate(intercept = c(0.41, 0.41, 0.74, 0.74, 0.39, 0.39), 
         slope = c(-0.11, -0.11, -0.11, -0.11, -0.11, -0.11))  

bat_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(Habitat, DietGroup))) + 
  geom_jitter(alpha = 0.2, width = 0.1) + 
  geom_ribbon(data = filter(CI_df_bat, !D95q, Dbin < 0), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(Habitat~DietGroup) + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() + theme(legend.position = "none")



######

## SD plot FIG 4 ----
sd <- list(bat_winner, bat_nt_winner, bird_winner, bird_nt_winner) %>% 
  map(~data.frame(summary(.x)$random) %>% select(1:4) %>% 
        setNames(c("Estimate", "SE", "l-95CI", "u-95CI")) %>% rownames_to_column()) %>% 
  setNames(c("bat_full", "bat_NT", "bird_full", "bird_NT")) %>% 
  bind_rows(.id = "model") %>% 
  separate(model, into = c("taxon", "model"), sep = "_") %>% 
  mutate(Habitat = ifelse(grepl("unaltered", rowname), "intact", "altered"), 
         DietGroup = ifelse(grepl("same", rowname), "intraguild", "control")) %>%
  rbind(filter(.,taxon == "bird", DietGroup == "control") %>% 
          mutate(DietGroup = "intraguild")) %>% select(-rowname)

sd %>% ggplot(aes(x = interaction(Habitat, DietGroup), y = Estimate, fill = interaction(Habitat, DietGroup))) + 
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
p1 <- bird_data_plot %>% ggplot(aes(x = theta_vl)) + geom_histogram() + facet_grid(status~DietGroup)
p2 <- bat_data_plot %>% ggplot(aes(x = theta_vl)) + geom_histogram() + facet_grid(status~DietGroup)
plot_grid(p2, p1, labels = c("A", "B"), ncol = 2)

## bat bias check
p1 <- bat_data_plot %>% ggplot(aes(x = log(N_sb), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.3) + geom_smooth(method = "loess") + facet_grid(status~DietGroup) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = expression(ln(N[set])))

p2 <- bat_data_plot %>% ggplot(aes(x = log(occ1), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess") + facet_grid(status~DietGroup) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of rarer species in pair")
plot_grid(p2, p1, labels = c("A", "B"), ncol = 2)

## bird bias check
p1 <- bird_data_plot %>% slice_sample(n = 100000) %>% ggplot(aes(x = log(N_sb), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.3) + geom_smooth(method = "loess") + facet_grid(status~DietGroup) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = expression(ln(N[set])))

p2 <- bird_data_plot %>% slice_sample(n = 50000) %>% ggplot(aes(x = log(occ1), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess") + facet_grid(status~DietGroup) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of rarer species in pair")
plot_grid(p2, p1, labels = c("A", "B"), ncol = 2)

# end model checks----



