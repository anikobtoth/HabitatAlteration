# JAGS Script by Drew Allen

library(rjags)
library(R2jags)
library(MCMCvis)
library(bayestestR)

# indep draw functions
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

recursive.indepdraw <- function(rows, shuf, data, pairs) {
  if(length(shuf) == 0) {return(rows)
  }else if(length(rows) > pairs) {return(rows)
  }else{
    row <- numeric()
    rsp <- character()
    for(i in 1:floor(length(shuf)/2)){
      trial <- which((data$Sp1 == shuf[(i*2)-1] & data$Sp2 == shuf[i*2]) |
                       (data$Sp2 == shuf[(i*2)-1] & data$Sp1 == shuf[i*2]))
      
      if(length(trial) > 0){
        row[i] <- trial
      }else{
        rsp <- c(rsp, shuf[(i*2)-1], shuf[i*2])
      }
    }
    row <- na.omit(row)
    print(row)
    print(rsp)
    if(length(row) == 0 && length(nsp) < 5) {return(rows)
    }else{rows <- c(rows, row)
    return(recursive.indepdraw(rows, shuf = sample(rsp), data = data[-rows,], pairs))
    } 
  }
}

indep_drw <- function(data, reps = 10, pairs){
  out <- list()
  for(j in 1:reps){
    print(j)
    splist <- unique(c(data$Sp1, data$Sp2))
    shuf <- sample(splist)
    rows <- numeric()
    rows <- recursive.indepdraw(rows, shuf, data, pairs)
    out[[j]] <- data[rows,]
    data <- data[-rows,]
  }
  
  return(out)
}

# omega
omega <- function(data, tax, related = FALSE, interaction = FALSE, bagging = FALSE, reps = 100, median = TRUE) {
  ## Prep data ##
  message("Preparing data")
  data <- contables %>% filter(taxon == tax)
  
  if(median){
    occs <- rbind(data %>% select(status, Sp1, presSp1) %>% setNames(c("status", "sp", "pres")), 
                  data %>% select(status, Sp2, presSp2) %>% setNames(c("status", "sp", "pres"))) %>% 
      unique()  
    meds <- occs %>% group_by(status) %>% summarise(med = median(pres)) %>% pull(med)
    
    alt <- occs %>% filter(status == "altered" & pres > meds[1]) %>% pull("sp")
    unalt <- occs %>% filter(status == "unaltered" & pres > meds[2]) %>% pull("sp")
    
    data <- rbind(data %>% filter(status == "altered" & Sp1 %in% alt & Sp2 %in% alt),
                  data %>% filter(status == "unaltered" & Sp1 %in% unalt & Sp2 %in% unalt))
    }
  
  if(related){
    data$diet.match[data$diet.match=='Related'] <- 'Same'   #related diets coded as competing
  }else{
    data$diet.match[data$diet.match=='Related'] <- 'Similar' # related diets excluded
  }
  data$sa <- ifelse(data$status=='altered',1,0)  # dummy variable for status
  data$ds <- ifelse(data$diet.match=='Same',1,0) # dummy variable for diet match
  data$xx <- as.numeric(interaction(factor(data$status),factor(data$diet.match)))
  
  ##  Define models ##
  message("Defining models")
  if(interaction){
    jags.params <- c(paste0('av.ln.omega[',1:4,']'),paste0('sigma[',1:4,']'))
    jags.inits <- function() {
      list(av.ln.omega=rnorm(4),sigma = runif(4))
    }
    jags.model <- function() {
      
      # Priors
      for (i in 1:4) {
        av.ln.omega[i] ~ dnorm(0,0.01)
        sigma[i]       ~ dunif(0, 6)       # hyperparameter sigma for random effect
        tau[i]        <- 1 / (sigma[i] * sigma[i]) # convert from sd to precision for dnorm
      }
      
      # Likelihood
      for (i in 1:N) {
        presBoth[i]  ~ dhyper(presSp1[i],absSp1[i],presSp2[i],exp(ln.omega[i]))
        ln.omega[i] <- av.ln.omega[xx[i]] + rand[i]
        rand[i]      ~ dnorm(0,tau[xx[i]])
      }
      
    }
    
  } else {
    jags.params <- c('av.ln.omega','altered','same',paste0('sigma[',1:4,']'))
    jags.inits <- function() {
      list(av.ln.omega=rnorm(1),altered=rnorm(1),same=rnorm(1),sigma = runif(4))
    }
    jags.model <- function() {
      # Priors
      av.ln.omega ~ dnorm(0,0.04)
      altered     ~ dnorm(0,0.04)
      same        ~ dnorm(0,0.04)
      for(i in 1:4) {  
        sigma[i]  ~ dunif(0,10)
        tau[i]   <- 1 / (sigma[i] * sigma[i])
      }
      # Likelihood
      for (i in 1:N) {
        presBoth[i]  ~ dhyper(presSp1[i],absSp1[i],presSp2[i],exp(ln.omega[i]))
        ln.omega[i] <- av.ln.omega + sa[i]*altered + ds[i]*same + rand[i]
        rand[i]      ~ dnorm(0,tau[xx[i]])
      }
    }
      
  }
  
  if(bagging){
    data    <- data[order(data$status,data$Sp1,data$Sp2),]
    temp_altered   <- data[data$status=='altered',]
    temp_unaltered <- data[data$status=='unaltered',]
    sp_altered <- length(unique(c(temp_altered$Sp1,temp_altered$Sp2)))
    sp_unaltered <- length(unique(c(temp_unaltered$Sp1,temp_unaltered$Sp2)))
    
    jags.fit <- list() # container for fitted models
    message("Running bagged estimates")
    
    set.seed(123)
     
    for(i in 1:reps) {
      
      ncat <- 2
      while(ncat != 4) {  ## Ensure all categories are represented in the random draw
        temp <- rbind(temp_altered[idx_rng(sp_altered),],
                      temp_unaltered[idx_rng(sp_unaltered),])
        if(!related) temp <- temp %>% filter(diet.match != "Similar")
        ncat <- length(unique(temp$xx))
      }
      
      if(interaction){jags.data <- list(N = nrow(temp),           # number of pairs
                        presSp1  = temp$presSp1,  # number of sites with species 1
                        absSp1   = temp$absSp1,   # number of sites without species 1
                        presSp2  = temp$presSp2,  # number of sites with species 2
                        presBoth = temp$presBoth, # number of sites with both species
                        xx = temp$xx)             # index
      
     }else{ jags.data <- list(N = nrow(temp),           # number of pairs
                        presSp1  = temp$presSp1,  # number of sites with species 1
                        absSp1   = temp$absSp1,   # number of sites without species 1
                        presSp2  = temp$presSp2,  # number of sites with species 2
                        presBoth = temp$presBoth, # number of sites with both species
                        sa = temp$sa,
                        ds = temp$ds,
                        xx = temp$xx)             # index
     }
      jags.fit[[i]] <- jags.parallel(data = jags.data, inits = jags.inits, 
                            parameters.to.save = jags.params, model.file = jags.model,
                            n.chains = 3, n.iter = 50000, n.burnin = 10000, n.thin = 100, n.cluster = detectCores()*.75)
    }
    
    if(interaction){
      out <- jags.fit %>% map(~as.mcmc(.) %>% as.matrix() %>% as.data.frame()) %>% 
        setNames(1:length(.)) %>% bind_rows(.id = "rep") %>%
        setNames(c("rep", "omega-Altered_Different", "omega-Unaltered_Different", 
                   "omega-Altered_Same", "omega-Unaltered_Same", "deviance-all", 
                   "sigma-Altered_Different", "sigma-Unaltered_Different", 
                   "sigma-Altered_Same", "sigma-Unaltered_Same"))
      
    }else{
      out <- jags.fit %>% map(~as.mcmc(.) %>% as.matrix() %>% as.data.frame()) %>% 
        setNames(1:length(.)) %>% bind_rows(.id = "rep")  %>%
        setNames(c("rep", "altered", "avg.ln.omega", "deviance", "same", "sigma-Altered_Different", 
                   "sigma-Unaltered_Different", "sigma-Altered_Same", "sigma-Unaltered_Same")) %>%
        mutate(`omega-Unaltered_Different` = avg.ln.omega,
               `omega-Unaltered_Same` = avg.ln.omega + same, 
               `omega-Altered_Different` = avg.ln.omega + altered, 
               `omega-Altered_Same` = avg.ln.omega + altered + same) %>% 
        select(-altered, -avg.ln.omega, -same)
    }
    out <- out %>% pivot_longer(cols = 2:10, names_to = "Group", values_to= "value") %>%
      separate(Group, into = c("parameter", "status_dietmatch"), sep = "-")
    
  }else{
    
    temp <- data
    if(!related) temp <- temp %>% filter(diet.match != "Similar")
    
    if(interaction){jags.data <- list(N = nrow(temp),           # number of pairs
                                      presSp1  = temp$presSp1,  # number of sites with species 1
                                      absSp1   = temp$absSp1,   # number of sites without species 1
                                      presSp2  = temp$presSp2,  # number of sites with species 2
                                      presBoth = temp$presBoth, # number of sites with both species
                                      xx = temp$xx)             # index
    
    }else{ jags.data <- list(N = nrow(temp),           # number of pairs
                             presSp1  = temp$presSp1,  # number of sites with species 1
                             absSp1   = temp$absSp1,   # number of sites without species 1
                             presSp2  = temp$presSp2,  # number of sites with species 2
                             presBoth = temp$presBoth, # number of sites with both species
                             sa = temp$sa,
                             ds = temp$ds,
                             xx = temp$xx)             # index
    }
    message("Running full estimates")
    set.seed(123)
    jags.fit <- jags.parallel(data = jags.data, inits = jags.inits, 
                                   parameters.to.save = jags.params, model.file = jags.model,
                                   n.chains = 3, n.iter = 50000, n.burnin = 10000, n.thin = 100, 
                                   n.cluster = detectCores()*.75)
  
  message("Formatting posteriors")
  
  out <- as.mcmc(jags.fit) %>% as.matrix() %>% as.data.frame()
  
  if(interaction){
    n <- c("omega-Altered_Different", "omega-Unaltered_Different", 
           "omega-Altered_Same", "omega-Unaltered_Same", "deviance-all", 
           "sigma-Altered_Different", "sigma-Unaltered_Different", 
           "sigma-Altered_Same", "sigma-Unaltered_Same")
    out <- out %>% setNames(n) 
      
  }else{
    n <- c("altered", "avg.ln.omega", "deviance", "same", "sigma-Altered_Different",
           "sigma-Unaltered_Different", "sigma-Altered_Same", "sigma-Unaltered_Same")
    out <- out %>% setNames(n) %>%  
     mutate(`omega-Unaltered_Different` = avg.ln.omega,
            `omega-Unaltered_Same` = avg.ln.omega + same,
            `omega-Altered_Different` = avg.ln.omega + altered,
            `omega-Altered_Same` = avg.ln.omega + altered + same) %>%
     select(-altered, -avg.ln.omega, -same)
 
  }
  out <- out %>% pivot_longer(cols = 1:9, names_to = "Group", values_to= "value") %>%
    separate(Group, into = c("parameter", "status_dietmatch"), sep = "-")
  }
  
  return(out)
}

######
######
######
