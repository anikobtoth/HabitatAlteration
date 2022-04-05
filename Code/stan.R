library(BiasedUrn)
library(dplyr)
library(rstan)

# prepare data from contables 
stan_data_fun <- function(df, medn) {
  if(medn) {
    occs <- rbind(df %>% select(status, Sp1, presSp1) 
                  %>% setNames(c("status", "sp", "pres")), 
                  df %>% select(status, Sp2, presSp2) 
                  %>% setNames(c("status", "sp", "pres"))) %>% 
      unique()  
    meds <- occs %>% group_by(status) %>% 
      summarise(med = median(pres)) %>% pull(med)
    alt <- occs %>% filter(status == "altered" & pres > meds[1]) %>% 
      pull("sp")
    unalt <- occs %>% filter(status == "unaltered" & pres > meds[2]) %>% 
      pull("sp")
    df <- rbind(df %>% filter(status == "altered" 
                                  & Sp1 %in% alt & Sp2 %in% alt),
                  df %>% filter(status == "unaltered" 
                                  & Sp1 %in% unalt & Sp2 %in% unalt))
  }
  df$diet.match[df$diet.match=='Related'] <- 'Similar' # related diets excluded
  df$sa <- ifelse(df$status=='altered',1,0)  # dummy variable for status
  df$ds <- ifelse(df$diet.match=='Same',1,0) # dummy variable for diet match
  df$gg <- as.numeric(interaction(factor(df$sa),factor(df$ds)))
  df <- df %>% filter(diet.match != "Similar")
  
  N_site <- df$presSp1[1] + df$absSp1[1]
  df$sb <- df$presBoth
  df$s1 <- apply(df[c('presSp1','presSp2')],1,min)
  df$s2 <- apply(df[c('presSp1','presSp2')],1,max)
  df$sb_min <- apply(df[c('s1','s2')],1,function(x) max(0,sum(x) - N_site))
  occ <- unique(df[c('s1','s2','sb_min')])
  occ <- occ[order(occ$s1,occ$s2),]
  occ$oo <- 1:nrow(occ)
  df <- merge(occ,df)
  rr <- unique(df[c('oo','gg')])
  rr <- rr[order(rr$oo,rr$gg),]
  rr$rr <- 1:nrow(rr)
  df <- merge(rr,df)
  df2 <- unique(df[c('oo','s1','s2','gg','rr','sb')])
  df2 <- df2[order(df2$rr,df2$sb),]
  df2$nsb <- 0
  for(i in 1:nrow(df2)) {
    df2$nsb[i] <- nrow(df[df$rr==df2$rr[i]&
                            df$sb==df2$sb[i],])
  }
  return(list(N_site = N_site,
              N_dat = nrow(df2),
              dat = df2[c('oo','sb','nsb','rr')],
              N_occ = nrow(occ),
              occ = occ[c('s1','s2')],
              occ_idx = cbind(tapply(1:nrow(df2),df2$oo,min),
                              tapply(1:nrow(df2),df2$oo,max)),
              N_gg = max(df2$gg),
              N_rr = max(df2$rr),
              gg = unique(df2[c('rr','gg')])[,2],
              N_posthoc = max(df2$gg)*(max(df2$gg)-1)/2))
}

# bayesian model
stan_theta <- function(stan_data){

  stan_code <- 
    '
functions {

  // PMF for Fisher noncentral hypergeometric distribution
  real fnchypergeo_lpmf(int[] sb, int[] nsb, vector theta,
                        int[] occ, int N_site) {
    vector[1 + min(occ) - (sum(occ) > N_site ? sum(occ) - N_site : 0)] u;
    vector[1 + min(occ) - (sum(occ) > N_site ? sum(occ) - N_site : 0)] ll;
    vector[num_elements(sb)] ll_summand;
    int sb_min = sum(occ) > N_site ? sum(occ) - N_site : 0;
    int sb_max = min(occ);
    for(i in sb_min:sb_max) {
        u[i - sb_min + 1] = i;
        ll[i - sb_min + 1] = lchoose(occ[1], i) + 
                  lchoose(N_site - occ[1], occ[2] - i);
    }
    for(i in 1:num_elements(nsb)) {
      ll_summand[i] = nsb[i]*(ll[sb[i] - sb_min + 1] + sb[i]*theta[i] - 
                              log_sum_exp(ll + u*theta[i]));
    }    
    return sum(ll_summand);
  }
  
}

data {

  int N_site; // # of sites
  int N_dat; // # of rows in data table (dat)
  int dat[N_dat,4]; // dat columns: oo, sb, nsb, rr
  int N_occ; // # of rows in occupancy table (occ)
  int occ[N_occ,2]; // occ columns: s1, s2
  int occ_idx[N_occ,2]; // row index to occ for dat
  int N_gg; // # of groups
  int N_rr; // # of random effects
  int gg[N_rr]; // group assigments for random effects
  int N_posthoc; // number of posthoc comparisons

}

parameters {

  vector[N_gg] mu; // means
  vector<lower=0>[N_gg] sigma; // standard deviations
  vector[N_rr] theta; // random effects

}

model {

  // priors
  theta ~ normal(mu[gg], sigma[gg]);
  mu ~ normal(0,10);
  for(i in 1:N_gg) sigma[i] ~ normal(0,10) T[0,];

  // loglikelihood
  for (i in 1:N_occ) {
    target += fnchypergeo_lpmf(dat[occ_idx[i,1]:occ_idx[i,2],2] | 
                               dat[occ_idx[i,1]:occ_idx[i,2],3],
                               theta[dat[occ_idx[i,1]:occ_idx[i,2],4]],
                               occ[i], N_site);
  }

}

generated quantities {
  
  vector[N_posthoc] posthoc_mu;
  vector[N_posthoc] posthoc_sigma;
  int pos = 0;
  
  for(i in 1:(N_gg - 1)) {
    for(j in (i+1):N_gg) {
      pos += 1;
      posthoc_mu[pos] = mu[i] - mu[j];
      posthoc_sigma[pos] = sigma[i] - sigma[j];
    }
  }

}
'
# fit model
stan_fit <- stan(model_code=stan_code,
                 data= stan_data, cores = detectCores()*0.5) 
return(stan_fit)
}

# maximum likelihood model
stan_theta_ML <- function(stan_data, stan_fit){
  # Fit ML model
  ml_stan_code <- 
    '
functions {

  // PMF for Fisher noncentral hypergeometric distribution
  real fnchypergeo_lpmf(int[] sb, int[] nsb, vector theta,
                        int[] occ, int N_site) {
    vector[1 + min(occ) - (sum(occ) > N_site ? sum(occ) - N_site : 0)] u;
    vector[1 + min(occ) - (sum(occ) > N_site ? sum(occ) - N_site : 0)] ll;
    vector[num_elements(sb)] ll_summand;
    int sb_min = sum(occ) > N_site ? sum(occ) - N_site : 0;
    int sb_max = min(occ);
    for(i in sb_min:sb_max) {
        u[i - sb_min + 1] = i;
        ll[i - sb_min + 1] = lchoose(occ[1], i) + 
                  lchoose(N_site - occ[1], occ[2] - i);
    }
    for(i in 1:num_elements(nsb)) {
      ll_summand[i] = nsb[i]*(ll[sb[i] - sb_min + 1] + sb[i]*theta[i] - 
                              log_sum_exp(ll + u*theta[i]));
    }    
    return sum(ll_summand);
  }
  
}

data {

  int N_site; // # of sites
  int N_dat; // # of rows in data table (dat)
  int dat[N_dat,4]; // dat columns: oo, sb, nsb, rr
  int N_occ; // # of rows in occupancy table (occ)
  int occ[N_occ,2]; // occ columns: s1, s2
  int occ_idx[N_occ,2]; // row index to occ for dat
  int N_gg; // # of groups (IGNORED)
  int N_rr; // # of FIXED effects
  int gg[N_rr]; // group assigments for fixed effects (IGNORED)
  int N_posthoc; // number of posthoc comparisons (IGNORED)

}

parameters {

  vector[N_rr] omega; // FIXED effects

}

model {

  // loglikelihood
  for (i in 1:N_occ) {
    target += fnchypergeo_lpmf(dat[occ_idx[i,1]:occ_idx[i,2],2] | 
                               dat[occ_idx[i,1]:occ_idx[i,2],3],
                               omega[dat[occ_idx[i,1]:occ_idx[i,2],4]],
                               occ[i], N_site);
  }

}
'
ml_model <- optimizing(stan_model(model_code=ml_stan_code),
                       data = stan_data,hessian=TRUE,
                       init = list(theta=summary(stan_fit,'theta')[[1]][,1]))

return(ml_model)
}

# summary stats
stan_summary <- function(stan_data, stan_fit, ml_model){
  # Summary data for plots
  stan_sum <- cbind.data.frame(unique(stan_data$dat[c('rr','oo')]),
                               gg = stan_data$gg,
                               theta_bayes = summary(stan_fit,'theta')$summary[,1],
                               theta_ml = ml_model$par)
  stan_sum$log_n <- log(tapply(stan_data$dat$nsb,stan_data$dat$rr,sum))
  temp <- stan_data$occ
  temp$oo <- 1:nrow(temp)
  stan_sum <- merge(stan_sum,temp)
  stan_sum$max_sb <- tapply(stan_data$dat$sb,stan_data$dat$rr,max)
  stan_sum$min_sb <- tapply(stan_data$dat$sb,stan_data$dat$rr,min)
  stan_sum$sb_min <- apply(stan_sum[c('s1','s2')],1,
                           function(x) max(0,sum(x)-stan_data$N_site))
  stan_sum$sb_max <- stan_sum$s1
  stan_sum$theta_ml[stan_sum$max_sb==stan_sum$sb_min] <- -Inf
  stan_sum$theta_ml[stan_sum$min_sb==stan_sum$sb_max] <- Inf
  
  # difference between Bayes estimate and maximum likelihood estimate
  stan_sum$diff_ml <- stan_sum$theta_bayes - stan_sum$theta_ml
  stan_sum$diff_ml[stan_sum$diff_ml == -Inf] <- -4
  stan_sum$diff_ml[stan_sum$diff_ml == Inf] <- 4
  stan_sum$diff_mean <- stan_sum$theta_ml-summary(stan_fit,'mu')[[1]][,1][stan_sum$gg]
  return(stan_sum)
}

# visualise shrinkage by comparing bayesian and ml estimates
stan_shrinkage_plots <- function(stan_sum){
  par(mfrow = c(2,1), mar = c(3.5,4,1,1))
  
  # difference between maximum likelihood estimate and group mean
  plot(stan_sum$diff_mean,stan_sum$diff_ml,col=alpha('black',0.1),pch=16,
       xlab='',ylab='',cex=0.75,ylim=c(-2.25,2.25), las = 1)
  mtext(expression("deviation from group mean "~(hat(theta[i])-bar(theta[g(i)]))),
        side=1,line=2.5,cex=1.1)
  mtext(expression("shrinkage "~(theta[i] - hat(theta[i]))),side=2,
        line=2.5,cex=1.1)
  # shrinkage plotted against number of pairs per occupancy group
  plot(stan_sum$log_n ,stan_sum$diff_ml,pch=16,col=alpha('black',0.1),axes=FALSE,
       xlab='',las=3,ylab='',cex=0.75, las = 1)
  abline(h=0,lty=2)
  mtext(expression(number~of~pairs),side=1,line=2.5,cex=1.1)
  mtext(expression("shrinkage "~(theta[i] - hat(theta[i]))),side=2,
        line=2.5,cex=1.1)
  axis(1,at=log(10^(0:3)),labels=parse(text=paste0('10^{',0:3,'}')))
  axis(2,at=c(-4,-2,-1,0,1,2,4),expression(-infinity,-2,-1,0,1,2,infinity), las = 1)
  box()
  # Add break lines along y-axis     
  axis.break(2, -3, style = "zigzag",brw=0.03)
  axis.break(2, 3, style = "zigzag",brw=0.03)
}

# Check for bias in theta estimates
# summary plots
ggplot(stan_sum, aes(x = s1, y = theta_bayes)) + 
  geom_point(size = .5, alpha = .5) + facet_wrap(~gg) + 
  geom_smooth(method = "loess") + 
  geom_line(aes(y = 0), lty = 2) + 
  labs(x = "occupancy of rarer species in pattern")

ggplot(stan_sum, aes(x = log_n, y = theta_bayes)) + 
  geom_point(size = .5, alpha = .5) + facet_wrap(~gg) + 
  geom_smooth(method = "loess") + 
  geom_line(aes(y = 0), lty = 2) + 
  labs(x = "log(number of pairs)")

# Check distribution of estimates (normal distribution is good)
ggplot(stan_sum, aes(x =theta_bayes)) + 
  geom_histogram(bins = 12) + 
  facet_wrap(~gg) + labs(x = "theta")

## SANITY CHECKS ####
# Log-likelihood function for NHD
log_sum_exp <- function(x) log(sum(exp(x - max(x)))) + max(x)
fnchypergeo_lpmf <- function(sb, s1, s2, n_site, omega) {
  sb_min <- max(c(0,s1 + s2 - n_site))
  sb_max <- min(c(s1,s2))
  u <- ll <- ans <- numeric()
  for(i in sb_min:sb_max) {
    u[i - sb_min + 1] = i;
    ll[i - sb_min + 1] = lchoose(s1, i) + lchoose(n_site - s1, s2 - i)
  }
  if(length(sb) > 1) {
    for(i in 1:length(sb)) {
      ans[i] = ll[sb[i] - sb_min + 1] + sb[i]*omega - log_sum_exp(ll + u*omega)
    }
  } else {
    for(i in 1:length(omega)) {
      ans[i] = ll[sb - sb_min + 1] + sb*omega[i] - log_sum_exp(ll + u*omega[i])
    }
  }    
  return(ans)
}

# Check log-likelihood function
# sb <- 2:4; s1 <- 40; s2 <- 4; n_site <- 42; psi <- 1.5
# log(dFNCHypergeo(sb,s1,n_site - s1,s2,psi)) # BiasedUrn
# fnchypergeo_lpmf(sb,s1,s2,n_site,log(psi)) # our function

# Score function
fnchypergeo_score <- function(sb, s1, s2, n_site, omega) {
  sb_min <- max(c(0,s1 + s2 - n_site))
  sb_max <- min(c(s1,s2))
  u <- ll <- ans <- numeric()
  for(i in sb_min:sb_max) {
    u[i - sb_min + 1] = i;
    ll[i - sb_min + 1] = lchoose(s1, i) + lchoose(n_site - s1, s2 - i)
  }
  ans <- numeric()
  if(length(sb) > 1) {
    for(i in 1:length(sb)) {
      ans[i] <- sb[i] - 
        exp(log_sum_exp(log(u) + ll + u*omega) - log_sum_exp(ll + u*omega))
    }
  } else {
    for(i in 1:length(omega)) {
      ans[i] = sb - 
        exp(log_sum_exp(log(u) + ll + u*omega[i]) - log_sum_exp(ll + u*omega[i]))
    }
  }       
  return(ans)
}

# Check score function
# (fnchypergeo_lpmf(1, 5, 20,50,1.11)-fnchypergeo_lpmf(1, 5, 20,50,1.09))/(1.11-1.09)
# fnchypergeo_score(1, 5, 20,50,1.1)

# Fisher Information function
fnchypergeo_FI <- function(s1, s2, n_site, omega) {
  sb_min <- max(c(0,s1 + s2 - n_site))
  sb_max <- min(c(s1,s2))
  ans<- numeric()
  u <- ll <- ans <- numeric()
  for(i in sb_min:sb_max) {
    u[i - sb_min + 1] = i;
    ll[i - sb_min + 1] = lchoose(s1, i) + lchoose(n_site - s1, s2 - i)
  }
  for(i in 1:length(omega)) {
    ans[i] <-
      exp(log_sum_exp(2*log(u) + ll + u*omega[i])-log_sum_exp(ll + u*omega[i])) - 
      exp(log_sum_exp(log(u) + ll + u*omega[i])-log_sum_exp(ll + u*omega[i]))^2
  }
  return(ans)
}

# Check Fisher Information function
fnchypergeo_FIb <- function(s1, s2, n_site, omega) {
  sb_min <- max(c(0,s1 + s2 - n_site))
  sb_max <- min(c(s1,s2))
  ans<- numeric()
  for(i in 1:length(omega)) {
    scores <- fnchypergeo_score(sb_min:sb_max, s1, s2, n_site, omega[i])^2
    probs  <- exp(fnchypergeo_lpmf(sb_min:sb_max, s1, s2, n_site, omega[i]))
    ans[i] <- sum(scores*probs)
  }
  return(ans)
}
# fnchypergeo_FI(1,5,50,-5:5)
# fnchypergeo_FIb(1,5,50,-5:5)


######
#####
####