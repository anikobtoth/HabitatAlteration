# JAGS Script by Drew Allen

library(rjags)
library(R2jags)
library(MCMCvis)
library(bayestestR)

data <- contables %>% filter(taxon == "bird")

data <- data[order(data$status,data$Sp1,data$Sp2),]
#data$diet.match[data$diet.match=='Related'] <- 'Same'
data$diet.match[data$diet.match=='Related'] <- 'Similar'
data$xx <- as.numeric(interaction(factor(data$status),factor(data$diet.match))) ## categories for analysis

#n <- paste(data$Sp1, data$Sp2, data$status, sep = "-")

# pair adjacent species in permuted vector of species and find row numbers
idx_rng <- function(sp) {
  x <- sample(1:sp)
  nidx <- floor(sp/2)
  y <- 0
  z <- 0
  for (i in 2:(sp - 1)) {
    z[i] = z[i - 1] + sp - i + 1
  }
  for (i in 0:(nidx-1)) {
    y[i+1] = z[min(x[i*2+1],x[i*2+2])] + abs(x[i*2+1] - x[i*2+2])
  }
  return(sort(y))
}

######## MODEL 1: ESTIMATE OMEGA USING FISHERS NONCENTRAL HYPERGEOMETRIC #######
jags.params <- c(paste0('av.ln.omega[',1:4,']'),"sigma")

jags.inits <- function() {
  list(av.ln.omega=rnorm(4),sigma = runif(1))
}

jags.model <- function() {
  
  # Priors
  for (i in 1:4) {
    av.ln.omega[i] ~ dnorm(0,0.001)
  }	
  sigma       ~ dunif(0, 100)       # hyperparameter sigma for random effect
  tau        <- 1 / (sigma * sigma) # convert from sd to precision for dnorm
  
  # Likelihood
  for (i in 1:N) {
    presBoth[i]  ~ dhyper(presSp1[i],absSp1[i],presSp2[i],exp(ln.omega[i]))
    ln.omega[i] <- av.ln.omega[xx[i]] + rand[i]
    rand[i]      ~ dnorm(0,tau)
  }
  
}

######### MODEL 2: ALLOW SIGMA TO DIFFER ##########
jags.params <- c(paste0('av.ln.omega[',1:4,']'),paste0('sigma[',1:4,']'))

jags.inits <- function() {
  list(av.ln.omega=rnorm(4),sigma = runif(4))
}

jags.model <- function() {
  
  # Priors
  for (i in 1:4) {
    av.ln.omega[i] ~ dnorm(0,0.01)
    sigma[i]       ~ dunif(0, 8)       # hyperparameter sigma for random effect
    tau[i]        <- 1 / (sigma[i] * sigma[i]) # convert from sd to precision for dnorm
  }
  
  # Likelihood
  for (i in 1:N) {
    presBoth[i]  ~ dhyper(presSp1[i],absSp1[i],presSp2[i],exp(ln.omega[i]))
    ln.omega[i] <- av.ln.omega[xx[i]] + rand[i]
    rand[i]      ~ dnorm(0,tau[xx[i]])
  }
  
}

###### FIT MODEL ############## 

temp_altered   <- data[data$status=='altered',]
temp_unaltered <- data[data$status=='unaltered',]
sp_altered <- length(unique(c(temp_altered$Sp1,temp_altered$Sp2)))
sp_unaltered <- length(unique(c(temp_unaltered$Sp1,temp_unaltered$Sp2)))


jags.fit <- list() # container for fitted models

for(i in 1:10) {

  ncat <- 2
  while(ncat != 4) {  ## Ensure all categories are represented in the random draw
    temp <- rbind(temp_altered[idx_rng(sp_altered),],
                  temp_unaltered[idx_rng(sp_unaltered),])
    temp <- temp %>% filter(diet.match != "Similar")
    ncat <- length(unique(temp$xx))
  }
  
  jags.data <- list(N = nrow(temp),           # number of pairs
                    presSp1  = temp$presSp1,  # number of sites with species 1
                    absSp1   = temp$absSp1,   # number of sites without species 1
                    presSp2  = temp$presSp2,  # number of sites with species 2
                    presBoth = temp$presBoth, # number of sites with both species
                    xx = temp$xx)             # index
  
  set.seed(123)
  jags.fit[[i]] <- jags(data = jags.data, inits = jags.inits, 
                        parameters.to.save = jags.params, model.file = jags.model,
                        n.chains = 3, n.iter = 50000, n.burnin = 10000, n.thin = 100)
  
}

#### EXAMINE Posteriors ########
out <- jags.fit %>% map(~as.mcmc(.) %>% as.matrix() %>% as.data.frame()) %>% 
  setNames(1:length(.)) %>% bind_rows(.id = "rep") %>%
  setNames(c("rep", "omega-Altered_Different", "omega-Unaltered_Different", 
             "omega-Altered_Same", "omega-Unaltered_Same", "deviance-all", 
             "sigma-Altered_Different", "sigma-Unaltered_Different", 
             "sigma-Altered_Same", "sigma-Unaltered_Same"))

out <- out %>% pivot_longer(cols = 2:10, names_to = "Group", values_to= "avg.ln.omega") %>% 
  separate(Group, into = c("parameter", "status_dietmatch"), sep = "-")

ggplot(out %>% filter(!parameter %in% c("deviance")), aes(x = avg.ln.omega, col = status_dietmatch)) + 
  geom_density() + facet_grid(parameter~., scales = "free_x")
  
# credible interval for ln.omega includes 0
# lots of variability among species pairs, as indicated by the size of sigma versus ln.omega
#jags.fit

# trace plots look good
#traceplot(jags.fit, mfrow = c(1, 2))

# # Extract posterior samples for parameters
# jags.fit.mcmc <- as.mcmc(jags.fit) %>% as.matrix() %>% as.data.frame() %>%
#   setNames(c("av.ln.omega", "deviance", 
#              n[c(1:length(n)) %>% as.character() %>% order()], "sigma")) %>% 
#   mutate(rep = 1:1200)
# jags.out <-jags.fit.mcmc %>% select(-av.ln.omega, -deviance, -sigma) %>% melt(id.vars = "rep") %>% 
#   separate(variable, into = c("Sp1", "Sp2", "status"), sep = "-")
# jags.hyper <- jags.fit.mcmc %>% select(av.ln.omega, deviance, sigma, rep) %>% melt(id.vars = "rep")
# hist(jags.fit.mcmc[,'av.ln.omega'])
# hist(jags.fit.mcmc[,'sigma'])
# hist(jags.fit.mcmc[,'ln.omega[1]'])
# 
# test <- jags.out %>% group_by(variable) %>% summarise(omega = mean(value, na.rm = TRUE))
# test1 <- out %>% filter(type == "observed" & Taxon_status == "bat_altered") %>% group_by(id) %>% summarise(fetmp = mean(Z.Score, na.rm = TRUE))
# 
# g <- merge(test, test1, by.x = "variable", by.y = "id", all.y = TRUE)
# g <- g %>% separate(variable, into = c("Sp1", "Sp2"), sep = "-")
# g <- merge(g, data) 
# ggplot(g, aes(x = presSp1, y = presSp2, col = x)) + geom_point()
# 
# plot3d(g[,c("presSp1", "presSp2", "presBoth")], col = g$x)
# 
# 
# 
# ## out_omega
# out_omega <- jags.out %>% filter(!variable %in% c("av.ln.omega", "sigma", "deviance")) %>% 
#   separate(variable, into = c("Sp1", "Sp2"), sep = "-")
# 
# out_omega <- diet_cat(out_omega, spp, related = FALSE)
# 

##### other categorical variables ######
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

######
######
######
