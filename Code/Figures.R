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
library(hexbin)
library(viridis)
library(gganimate)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
##  MAIN TEXT 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

colors <- c("#ADADAD", "#525252", "#2F99DC", "#00487C", "#C595CE", "#480355")

# calculate posterior slopes and intercepts
calc_parameters <- function(post){
  params <- list()
  c <- 0
  coeffs <- grep(names(post), pattern = "^b_", value = TRUE)
  for(i in c("low", "medium", "high")){
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
    group_by(DietGroup, Dbin) %>%
    summarise(theta_025 = quantile(theta, 0.025),
              theta_975 = quantile(theta, 0.975),
              theta_mu  = mean(theta)) %>%
    full_join(plotdata %>% group_by(DietGroup) %>% 
                summarise(d025 = quantile(D, 0.025), d975 = quantile(D, 0.975))) %>%
    mutate(D95q = Dbin > d025 & Dbin < d975)
}


# BAT theta posterior mean and CI calculations ----
tax = "bat"
bat_data <- stan_data_fun(filter(contables, taxon == tax))[[1]]
data_plot <- stan_data_fun(filter(contables, taxon == tax))[[2]]
bat_post <- as.data.frame(bat_winner)  

bat_data_plot <- bat_post %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bat_data) %>% full_join(data_plot) %>% 
  mutate(DietGroup =  fct_relevel(DietGroup, "low", "medium", "high"),
         Habitat = recode(Habitat, "unaltered" = "intact"))

CI_df_bat <- bat_post %>% calc_parameters() %>% CIcalc(data_plot) %>% 
  mutate(DietGroup = recode(DietGroup, "low" = "control", "medium" = "intersecting", "high" = "intraguild") %>% fct_relevel("low", "medium", "high"))

linefit <- expand.grid(DietGroup = c("high", "medium", "low")) %>% 
  mutate(intercept = c(mean(bat_post$b_Intercept+bat_post$b_DietGroupsame), 
                       mean(bat_post$b_Intercept+bat_post$b_DietGroupintersecting),
                       mean(bat_post$b_Intercept)), 
         slope = c(mean(bat_post$b_PhyloD))) 

## BAT PLOT FIG 2 ####
# points
bat_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(Habitat, DietGroup))) + 
  geom_jitter(alpha = 0.3, width = 0.15, size = 0.5) + 
  geom_ribbon(data = filter(CI_df_bat, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(Habitat~DietGroup) + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() + theme(legend.position = "none")

# hexbins
bat_data_plot %>% ggplot(aes(x = jitter(D), y = theta_vl)) + 
  geom_hex() + facet_grid(.~DietGroup) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bat, !D95q, Dbin < 0), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "orange", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() 

# hexbin, no habitat dimension
F2A <- bat_data_plot %>% mutate(DietGroup = recode(DietGroup, "different" = "low", "intersecting" = "medium", "intraguild" = "high")) %>% 
  ggplot(aes(x = D, y = theta_vl)) + 
  geom_hex() + facet_grid(.~DietGroup) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bat, !D95q, Dbin < 0), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.4) +
  geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "orange", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit, aes(intercept = intercept, slope = slope),lwd = 0.4) +
  #geom_ribbon(data = filter(CI_df_bat, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group") + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light()
######

# BIRD theta posterior mean and CI calculations ----
# calculate posterior slope (a) and intercept (b) for each of 4 groups
tax = "bird"
bird_data <- stan_data_fun(filter(contables, taxon == tax))[[1]]
data_plot <- stan_data_fun(filter(contables, taxon == tax))[[2]]
bird_post <- as.data.frame(bird_winner)
CI_df_bird <- bird_post %>% calc_parameters() %>% CIcalc(data_plot) %>%
  mutate(DietGroup = recode(DietGroup, "low" = "control", "medium" = "intersecting", "high" = "intraguild") %>% fct_relevel("low", "medium", "high"))

bird_data_plot <- bird_post %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bird_data) %>% full_join(data_plot) %>% 
  mutate(DietGroup = recode(DietGroup, "low" = "control", "medium" = "intersecting", "high" = "intraguild") %>% fct_relevel("low", "medium", "high"), 
         Habitat = recode(Habitat, "unaltered" = "intact"))

## Bird plots FIG 3 ####

linefit_full <- expand.grid(DietGroup = c("high","medium", "low")) %>%
  mutate(intercept = mean(bird_post$b_Intercept),
         slope  = mean(bird_post$b_PhyloD))

#NT models linefit
 linefit_nt <- expand.grid(Habitat = c("intact", "altered"), DietGroup = c("high","medium", "low")) %>%
   mutate(intercept = c(0.37, 0.31, 0.37, 0.31, 0.37, 0.31),
          slope  = c(-0.09, -0.1, -0.09, -0.1, -0.09, -0.1))

bird_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(Habitat, DietGroup))) + 
  geom_jitter(alpha = 0.2, width = 0.15, size = 0.5) + 
  geom_ribbon(data = filter(CI_df_bird_full, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit_full, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  facet_grid(Habitat~DietGroup) + labs(col = "group") + 
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  scale_color_manual(values = colors) +
  theme_light() + theme(legend.position = "none")

birdfull <- bird_data_plot_full %>% ggplot(aes(x = jitter(D), y = theta_vl)) + 
  geom_hex() + facet_grid(Habitat~DietGroup) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bird_full, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.5) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "orange", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit_full, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  theme_light() 

birdnt <- bird_data_plot_nt %>% ggplot(aes(x = D, y = theta_vl)) + 
  geom_hex() + facet_grid(Habitat~DietGroup) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bird_nt, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.25) +
  geom_ribbon(data = filter(CI_df_bird_nt, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "yellow", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit_nt, aes(intercept = intercept, slope = slope), lty = 2,lwd = 0.4) +
  geom_ribbon(data = filter(CI_df_bird_nt, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  theme_light() 

bird_data_plot_diff <- anti_join(bird_data_plot_full, bird_data_plot_nt)
bird_data_plot <- list(retained = bird_data_plot_nt, turnover = bird_data_plot_diff) %>% 
  bind_rows(.id = "type")

bird_data_plot %>%
  ggplot(aes(x = theta_vl, y = as.factor(D%/%3), col = status)) + 
  geom_density_ridges(scale = 1, alpha = 0.4, bandwidth = .15, lwd = 1) + 
  facet_grid(type~DietGroup) +
  scale_color_manual(values = colors) +
  xlim(c(-1.1, 2.2)) +
  labs(y = "binned standardised phylogenetic distance", x = expression(theta))

#combined bird data NT and FULL
colours2 <- c("#FF4E33", "#000000")
bird_data_plot <- list(FULL = bird_data_plot_full, NT = bird_data_plot_nt) %% bind_rows(.id = model)
bird_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = model)) + 
     geom_jitter(alpha = 0.2, width = 0.15, size = 0.5) + 
     facet_grid(Habitat~DietGroup) +
     geom_abline(data = linefit_bird, aes(intercept = intercept, slope = slope, col = model), lty = 2,lwd = 0.4) +
     scale_color_manual(values = colours2) +
     labs(y = expression(theta), x = "standardised phylogenetic distance")

F2B <-  bird_data_plot_full %>% ggplot(aes(x = jitter(D), y = theta_vl)) + 
  geom_hex() + facet_grid(.~DietGroup) + scale_fill_viridis_c(trans = "log", labels = label_number(accuracy = 1)) + 
  geom_ribbon(data = filter(CI_df_bird_full, !D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), inherit.aes = FALSE, fill = "gray50", alpha = 0.5) +
  geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_025, ymax = theta_975), fill = "orange", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(data = linefit_full, aes(intercept = intercept, slope = slope),lwd = 0.4) +
  #geom_ribbon(data = filter(CI_df_bird_full, D95q), aes(x = Dbin, ymin = theta_mu, ymax = theta_mu), col = "black", inherit.aes = FALSE) +
  labs(col = "group", y = expression(theta), x = "standardised phylogenetic distance") + 
  theme_light() 

plot_grid(F2A, F2B, ncol = 1, align = "v", labels = c("A", "B"))
#######
##### No-Turnover plots #####

## BAT NO TURNOVER PLOT  ####
tax = "bat"
bat_data <- stan_data_fun(filter(contables, taxon == tax))[[1]]
data_plot <- stan_data_fun(filter(contables, taxon == tax))[[2]]
bat_post <- as.data.frame(bat_nt_winner)  
CI_df_bat <- bat_post %>% calc_parameters() %>% CIcalc(data_plot)

bat_data_plot <- as.data.frame(bat_nt_winner) %>% colMeans() %>% 
  data.frame(names(.)) %>% setNames(c("value", "col")) %>% 
  filter(grepl('theta', col)) %>% separate(col, into=c("name", "index"), sep = "\\[") %>% 
  pivot_wider() %>% cbind(bat_data) %>% full_join(data_plot) %>% 
  mutate(DietGroup = recode(DietGroup, "different" = "low", "same" = "high"),
         Habitat = recode(Habitat, "unaltered" = "intact"))

# bat plot
linefit <- expand.grid(Habitat = c("intact", "altered"), DietGroup = c("high", "medium", "low")) %>% 
  mutate(intercept = c(0.41, 0.36, 0.77, 0.72, 0.42, 0.37), 
         slope = c(-0.11, -0.11, -0.11, -0.11, -0.11, -0.11))  

bat_data_plot %>% ggplot(aes(x = D, y = theta_vl, col = interaction(Habitat, DietGroup))) + 
  geom_jitter(alpha = 0.2, width = 0.15, size = 0.5) + 
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
sd <- list(bat_winner, bat_nt_winner, bird_fx18_re20, bird_nt_winner) %>% 
  map(~data.frame(summary(.x)$random) %>% select(1:4) %>% 
        setNames(c("Estimate", "SE", "l-95CI", "u-95CI")) %>% rownames_to_column()) %>% 
  setNames(c("bat_Full", "bat_No-turnover", "bird_FULL", "bird_No-turnover")) %>% 
  bind_rows(.id = "model") %>% 
  separate(model, into = c("taxon", "model"), sep = "_") %>% 
  mutate(Habitat = ifelse(grepl("unaltered", rowname), "intact", ifelse(grepl("altered", rowname), "altered", "")), 
         DietGroup = ifelse(grepl("same", rowname), "within-guild", ifelse(grepl("intersecting", rowname), "intermediate", ifelse(grepl("different", rowname), "control", "" )))) 

colors <- c("#ADADAD", "#2F99DC", "#C595CE")

sd %>% ggplot(aes(x = interaction(Habitat, DietGroup, sep = ""), y = Estimate, fill = interaction(Habitat, DietGroup, sep = ""))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = `l-95CI`, ymax = `u-95CI`), width = 0.1, lwd = 0.5) +
  facet_wrap(taxon~model, scales = "free") + 
  theme_light() + labs(fill = "group") +
  scale_fill_manual(values = colors) + 
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

### Figure S4: PCoA ####
#ordinations, colored by alteration type
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

plot_grid(richdist + theme(legend.position = "none"), ordplot, ncol = 2) 

######

## Model checks Figs S5-S7 -----
library(cowplot)
## normality check
p1 <- bird_data_plot %>% ggplot(aes(x = theta_vl)) + geom_histogram(binwidth = .25) + facet_grid(status~DietGroup) + labs(x = expression(theta~"estimate"))
p2 <- bat_data_plot %>% ggplot(aes(x = theta_vl)) + geom_histogram(binwidth = .25) + facet_grid(status~DietGroup) + labs(x = expression(theta~"estimate"))
plot_grid(p2, p1, labels = c("A", "B"), ncol = 2)

## bat bias check
p1 <- bat_data_plot %>% ggplot(aes(x = log(N_sb), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.3) + geom_smooth(method = "loess") + 
  facet_grid(status~DietGroup, scales = "free_y") + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = expression(ln(N[set])))

p2 <- bat_data_plot %>% ggplot(aes(x = log(occ1), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess", se=F) + 
  facet_grid(status~DietGroup, scales = "free") + 
  ylim(c(0, 1.5))+
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of rarer species in pair")
plot_grid(p2, p1, labels = c("A", "B"), ncol = 1)

## bird bias check
p1 <- bird_data_plot %>% slice_sample(n = 100000) %>% ggplot(aes(x = log(N_sb), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.3) + geom_smooth(method = "loess") + facet_grid(status~DietGroup) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = expression(ln(N[set])))

p2 <- bird_data_plot %>% slice_sample(n = 50000) %>% ggplot(aes(x = log(occ1), y = sqrt(abs(theta_vl-theta_mu)))) + 
  geom_point(size = 0.8, alpha = 0.4) + geom_smooth(method = "loess") + facet_grid(status~DietGroup) + 
  labs(y = expression(paste(sqrt(abs(theta-bar(theta))))), 
       x = "logged occupancy of rarer species in pair")
plot_grid(p2, p1, labels = c("A", "B"), ncol = 1)

# end model checks----

# model tables for supplement ----
list(batfull = bat_winner, 
     batNT = bat_nt_winner, 
     birdfull = bird_winner, 
     birdNT = bird_nt_winner) %>% map(~summary(.x)$fixed %>% rownames_to_column("Variable")) %>% 
  bind_rows(.id = "model") %>% as_tibble() %>% write.csv("./Results/brms_model_coeffs.csv")

bloo <- readRDS("./stan/bat_randeff/bat_loo_randeff.rds") # done
bat_lc_rand <- bloo %>% #setNames(word(names(.), 3, 3, sep = "_")) %>% 
  loo_compare() %>% data.frame() %>% rownames_to_column() %>%
  left_join(data.frame(re_steps) %>% rownames_to_column(), by = "rowname") %>%
  select(random_id = rowname, `random effects structure` = re_steps,elpd_diff, se_diff) %>%
  mutate(z = elpd_diff/se_diff,
         p = pnorm(z))
write_csv(bat_lc_rand, "./Results/bat_lc_rand_table.csv")

bloo <- readRDS("./stan/bird_randeff/bird_loo_randeff.rds")
bird_lc_rand <- bloo %>% loo_compare() %>% data.frame() %>% rownames_to_column() %>%
  left_join(data.frame(re_steps) %>% rownames_to_column(), by = "rowname") %>%
  select(random_id = rowname, `random effects structure` = re_steps,elpd_diff, se_diff) %>%
  mutate(z = elpd_diff/se_diff,
         p = pnorm(z))
write_csv(bird_lc_rand, "./Results/bird_lc_rand_table.csv")


bloo <- readRDS("./stan/bat_fixedeff/bat_loo_fixdeff.rds") #done
bat_lc_fixd <- bloo %>% setNames(word(names(.), 2, 2, sep = "_")) %>% loo_compare() %>% data.frame() %>% rownames_to_column() %>%
  left_join(data.frame(fx_steps) %>% rownames_to_column(), by = "rowname") %>%
  select(fixed_id = rowname, `fixed effects structure` = fx_steps,elpd_diff, se_diff) %>%
  mutate(z = elpd_diff/se_diff,
         p = pnorm(z))
write_csv(bat_lc_fixd, "./Results/bat_lc_fixed_table.csv")


bloo <- readRDS("./stan/bird_fixedeff/bird_loo_fixedeff.rds") #done
bird_lc_fixd <- bloo %>% setNames(word(names(.), 2, 2, sep = "_")) %>% loo_compare() %>% data.frame() %>% rownames_to_column() %>%
  left_join(data.frame(fx_steps) %>% rownames_to_column(), by = "rowname") %>%
  select(fixed_id = rowname, `fixed effects structure` = fx_steps,elpd_diff, se_diff) %>%
  mutate(z = elpd_diff/se_diff,
         p = pnorm(z))
write_csv(bird_lc_fixd, "./Results/bird_lc_fixed_table.csv")


## End model tables for supplement ----

