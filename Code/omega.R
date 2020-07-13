# JAGS Script by Drew Allen

library(rjags)
library(R2jags)

data2 <- map(tables, map, cont_table) %>% map(bind_rows, .id = "status") %>% bind_rows(.id = "taxon")
data <- data2 %>% filter(taxon == "bird")

n <- paste(data$Sp1, data$Sp2, data$status, sep = "-")
jags.params <- c("av.ln.omega","sigma",paste0('ln.omega[',1:nrow(data),']'))

jags.inits <- function() {
	list(av.ln.omega=rnorm(1),sigma = runif(1))
}

jags.data <- list(N = nrow(data),           # number of pairs
                  presSp1 = data$presSp1,   # number of sites with species 1
                  absSp1  = data$absSp1,    # number of sites without species 1
                  presSp2 = data$presSp2,   # number of sites with species 2
                  presBoth = data$presBoth) # number of sites with both species

jags.model <- function() {

	# Priors
	av.ln.omega ~ dnorm(0, 0.01)      # average for log(omega)
	sigma       ~ dunif(0, 100)       # hyperparameter sigma for random effect
  tau        <- 1 / (sigma * sigma) # need to convert from sd to precision for dnorm

	# Likelihood
	for (i in 1:N) {
		presBoth[i]  ~ dhyper(presSp1[i],absSp1[i],presSp2[i],exp(ln.omega[i]))
	  ln.omega[i] <- av.ln.omega + rand[i]
	  rand[i]      ~ dnorm(0,tau)
	}

}

set.seed(123)
jags.fit <- jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.params, model.file = jags.model,
			   n.chains = 3, n.iter = 50000, n.burnin = 10000, n.thin = 100) 

# credible interval for ln.omega includes 0
# lots of variability among species pairs, as indicated by the size of sigma versus ln.omega
jags.fit

# trace plots look good
traceplot(jags.fit, mfrow = c(1, 2))

# Extract posterior samples for parameters
jags.fit.mcmc <- as.mcmc(jags.fit) %>% as.matrix() %>% as.data.frame() %>%
  setNames(c("av.ln.omega", "deviance", 
             n[c(1:length(n)) %>% as.character() %>% order()], "sigma")) %>% 
  mutate(rep = 1:1200)
jags.out <-jags.fit.mcmc %>% select(-av.ln.omega, -deviance, -sigma) %>% melt(id.vars = "rep") %>% 
  separate(variable, into = c("Sp1", "Sp2", "status"), sep = "-")
jags.hyper <- jags.fit.mcmc %>% select(av.ln.omega, deviance, sigma, rep) %>% melt(id.vars = "rep")
hist(jags.fit.mcmc[,'av.ln.omega'])
hist(jags.fit.mcmc[,'sigma'])
hist(jags.fit.mcmc[,'ln.omega[1]'])

test <- jags.out %>% group_by(variable) %>% summarise(omega = mean(value, na.rm = TRUE))
test1 <- out %>% filter(type == "observed" & Taxon_status == "bat_altered") %>% group_by(id) %>% summarise(fetmp = mean(Z.Score, na.rm = TRUE))

g <- merge(test, test1, by.x = "variable", by.y = "id", all.y = TRUE)
g <- g %>% separate(variable, into = c("Sp1", "Sp2"), sep = "-")
g <- merge(g, data) 
ggplot(g, aes(x = presSp1, y = presSp2, col = x)) + geom_point()

plot3d(g[,c("presSp1", "presSp2", "presBoth")], col = g$x)



## out_omega
out_omega <- jags.out %>% filter(!variable %in% c("av.ln.omega", "sigma", "deviance")) %>% 
  separate(variable, into = c("Sp1", "Sp2"), sep = "-")
out_omega$diet.Sp1 <- spp[out_omega$Sp1,"guild"]
out_omega$diet.Sp2 <- spp[out_omega$Sp2,"guild"]

related <- c("C-CI", "CI-IG", "NF-NI", "I-NI", "CI-NI", "CI-I", "CI-FI", "CI-IN","F-FN", "FI-NF", "FI-NI", 
             "FN-IN", "FI-I", "FI-FN", "FI-IN", "FN-N", "FG-FI", "FG-FN", "FG-IG", "FI-IG", "FG-G", "I-IN", "IN-N", 
             "N-NI", "N-NF", "F-FI", "FG-NF", "F-FG", "G-IG", "I-IG","IG-IN") 
out_omega <- out_omega[-which(paste(out_omega$diet.Sp1, out_omega$diet.Sp2, sep = "-") %in% related),]

#out_omega$diet.pair <- map2(out_omega$diet.Sp1, out_omega$diet.Sp2, function(x, y) c(x,y)) %>% map(sort) %>% map(paste, collapse = "-") %>% unlist()
out_omega$diet.match <- as.numeric(out_omega$diet.Sp1 == out_omega$diet.Sp2)
out_omega$diet.match[out_omega$diet.match == 0] <- "Different"
out_omega$diet.match[out_omega$diet.match == 1] <- "Same"

out_omega$cat.Sp1[out_omega$Sp1 %in% uniquesp] <- "Unique"
out_omega$cat.Sp1[out_omega$Sp1 %in% sharedsp] <- "Shared"
out_omega$cat.Sp2[out_omega$Sp2 %in% uniquesp] <- "Unique"
out_omega$cat.Sp2[out_omega$Sp2 %in% sharedsp] <- "Shared"

out_omega$cat.pair <- paste(out_omega$cat.Sp1, out_omega$cat.Sp2, sep = "-")
out_omega$cat.pair[out_omega$cat.pair == "Shared-Unique"] <- "Unique-Shared"

out_omega$cat.group[out_omega$cat.pair == "Unique-Unique"] <- "Unique"
out_omega$cat.group[out_omega$cat.pair == "Unique-Shared"] <- "Unique"
out_omega$cat.group[out_omega$cat.pair == "Shared-Shared"] <- "Shared"

### Cosmopolitan/restricted groupings #  
out_omega$cosmo.Sp1 <- cosmo[out_omega$Sp1,"group_abbr"]
out_omega$cosmo.Sp2 <- cosmo[out_omega$Sp2,"group_abbr"]
out_omega$cosmo.pair <- paste(out_omega$cosmo.Sp1, out_omega$cosmo.Sp2, sep = "-")

out_omega$cosmo.pair[out_omega$cosmo.pair == "cosmo-restr"] <- "restr-cosmo"
out_omega$cosmo.pair[out_omega$cosmo.pair == "synan-cosmo"] <- "cosmo-synan"
out_omega$cosmo.pair[out_omega$cosmo.pair == "synan-restr"] <- "restr-synan"

glm_omega <- out_omega %>% na.omit() %>% 
  group_by(rep, cosmo.pair, diet.match) %>% 
  summarise(mean = mean(value), median = median(value)) %>% 
  pivot_wider(id_cols = c("rep", "cosmo.pair"), names_from = "diet.match", values_from = "mean") 
