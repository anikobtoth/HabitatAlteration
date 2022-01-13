library(rstan)
library(BiasedUrn)

# data prep
data <- contables[contables$taxon=='bat',]
medn <- TRUE
if(medn) {
  occs <- rbind(data %>% select(status, Sp1, presSp1) 
                %>% setNames(c("status", "sp", "pres")), 
                data %>% select(status, Sp2, presSp2) 
                %>% setNames(c("status", "sp", "pres"))) %>% 
    unique()  
  meds <- occs %>% group_by(status) %>% 
    summarise(med = median(pres)) %>% pull(med)
  alt <- occs %>% filter(status == "altered" & pres > meds[1]) %>% 
    pull("sp")
  unalt <- occs %>% filter(status == "unaltered" & pres > meds[2]) %>% 
    pull("sp")
  data <- rbind(data %>% filter(status == "altered" 
                                & Sp1 %in% alt & Sp2 %in% alt),
                data %>% filter(status == "unaltered" 
                                & Sp1 %in% unalt & Sp2 %in% unalt))
}
data$diet.match[data$diet.match=='Related'] <- 'Similar' # related diets excluded
data$sa <- ifelse(data$status=='altered',1,0)  # dummy variable for status
data$ds <- ifelse(data$diet.match=='Same',1,0) # dummy variable for diet match
data$gg <- as.numeric(interaction(factor(data$status),factor(data$diet.match)))
temp <- data
temp <- temp %>% filter(diet.match != "Similar")
temp <- temp[temp$presSp1>1&temp$presSp2>1,]

stan_code <- 
'
functions {

  // loglikelihood
  real fnchypergeo_lpmf(int[] sb, int[] nsb, int[] occ, 
                        vector ln_omega, int N_site) {
    vector[num_elements(sb)] ll_num_summand;
    vector[num_elements(sb)] ll_den_summand;
    vector[occ[1] + 1] ll;
    vector[occ[1] + 1] u;
    for(i in 0:occ[1]) {
        u[i+1] = i;
        ll[i+1] = lchoose(occ[1], u[i+1]) + 
                  lchoose(N_site - occ[1], occ[2] - u[i+1]);
    }
    for(i in 1:num_elements(sb)) {
       ll_num_summand[i] = nsb[i]*ll[sb[i] + 1] + nsb[i]*sb[i]*ln_omega[i];
       ll_den_summand[i] = nsb[i]*log_sum_exp(ll + u*ln_omega[i]);
    }    
    return sum(ll_num_summand) - sum(ll_den_summand);
  }
  
}

data {

  int<lower=1> N; // number of combinations for occ/gg/sb
  int<lower=1> N_occ; // number of occupancy patterns for s1/s2
  int<lower=1> N_site; // number of sites
  int<lower=1> N_rr; // number of random effects (combinations of occ/gg)
  int<lower=1> N_gg; // number of groups
  int<lower=1,upper=N_gg> gg[N_rr]; // group assigment for each combination
  int<lower=1,upper=N_occ> oo[N]; // index to occ
  int<lower=0,upper=N_site> sb[N]; // sites with both species
  int<lower=1> nsb[N]; // frequency of each occ/gg/sb combination
  int<lower=1,upper=N_rr> rr[N]; // index to random effects
  int<lower=0,upper=N_site> occ[N_occ,2]; // occupancy table
  int<lower=1,upper=N> oo_l[N_occ]; // pairs index to occ for gg/sb/nsb (lower)
  int<lower=1,upper=N> oo_u[N_occ]; // pairs index to occ for gg/sb/nsb (upper)

}

parameters {

  vector[N_gg] mu; // means
  vector<lower=0>[N_gg] sigma; // standard deviations
  vector[N_rr] ln_omega; // random effects

}

model {

  // priors
  ln_omega ~ normal(mu[gg], sigma[gg]);
  for(i in 1:N_gg) sigma[i] ~ normal(0,5) T[0,];

  // loglikelihood
  for (i in 1:N_occ) {
    target += fnchypergeo_lpmf(sb[oo_l[i]:oo_u[i]] | nsb[oo_l[i]:oo_u[i]],
                               occ[i],ln_omega[rr[oo_l[i]:oo_u[i]]],N_site);
  }

}

generated quantities {

  vector[(N_gg)*(N_gg - 1)/2] posthoc_mu;
  vector[(N_gg)*(N_gg - 1)/2] posthoc_sigma;
  int pos = 0;
  
  for(i in 1:(N_gg - 1)) {
    for(j in (i+1):N_gg) {
      pos += 1;
      posthoc_mu[pos] = ln_omega[i] - ln_omega[j];
      posthoc_sigma[pos] = sigma[i] - sigma[j];
    }
  }

}
'

stan_data_prep <- function(df) {
  temp <- df
  df$sb <- df$presBoth
  df$s1 <- apply(df[c('presSp1','presSp2')],1,min)
  df$s2 <- apply(df[c('presSp1','presSp2')],1,max)
  occ <- unique(df[c('s1','s2')])
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
  return(list(N = nrow(df2),
              N_occ = nrow(occ),
              N_site = df$presSp1[1] + df$absSp1[1],
              N_rr = max(df2$rr),
              N_gg = max(df2$gg),
              gg = unique(df2[c('rr','gg')])[,2],
              oo = df2$oo,
              sb = df2$sb,
              nsb = df2$nsb,
              rr = df2$rr,
              occ = occ[c('s1','s2')],
              oo_l = tapply(1:nrow(df2),df2$oo,min),
              oo_u = tapply(1:nrow(df2),df2$oo,max)))
}
stan_data <- stan_data_prep(temp)

stan_fit <- stan(model_code=stan_code,data= stan_data)

# test of fnchypergeo_lpmf
log_sum_exp <- function(x) log(sum(exp(x - max(x)))) + max(x)
fnchypergeo <- function(sb,occ,ln_omega,N_site) {
  u <- ll <- ll_num_summand <- ll_den_summand <- numeric()
  for(i in 0:occ[1]) {
    u[i+1] = i
    ll[i+1] = lchoose(occ[1], u[i+1]) + 
              lchoose(N_site - occ[1], occ[2] - u[i+1])
  }
  for(i in 1:length(sb)) {
    ll_num_summand[i] = ll[sb[i] + 1] + sb[i]*ln_omega
    ll_den_summand[i] = log_sum_exp(ll + u*ln_omega)
  }    
  return(ll_num_summand -ll_den_summand)
}
sb <- 0:4; occ <- c(4,6); N_site <- 42; omega <- 1.5
log(dFNCHypergeo(sb,occ[1],N_site - occ[1],occ[2],omega)) # biased urn
fnchypergeo(sb,occ,log(omega),N_site) # our function

# Fisher Information
score_fnchypergeo <- function(sb,occ,ln_omega,N_site) {
  u <- ll <- numeric()
  for(i in 0:occ[1]) {
    u[i+1] = i
    ll[i+1] = lchoose(occ[1], u[i+1]) + 
      lchoose(N_site - occ[1], occ[2] - u[i+1])
  }
  ans <- numeric()
  for(i in 1:length(ln_omega)) {
   ans[i] <- sb - sum(u*exp(ll + u*ln_omega[i]))/
             sum(exp(ll + u*ln_omega[i]))
  }
  return(ans)
}
D_fnchypergeo(sb,occ,log(omega),N_site)
FI_fnchypergeo <- function(ln_omega,occ,N_site) {
  ans <- numeric()
  for(i in 1:length(ln_omega)) {
    p <- exp(fnchypergeo(0:occ[1],occ,ln_omega[i],N_site))
    d <- D_fnchypergeo(0:occ[1],occ,ln_omega[i],N_site)
    ans[i] <- sum(p*d^2)
  }
  return(ans)
}
o.sq <- seq(-10,10,.01)
plot(o.sq,FI_fnchypergeo(o.sq,c(2,2),N_site),las=1,
     xlab=expression(omega),type='l',ylim=c(0,2.6),
     ylab=expression(Fisher~Information~I[X](omega)))
lines(o.sq,FI_fnchypergeo(o.sq,c(4,4),N_site),
      lty=2)
lines(o.sq,FI_fnchypergeo(o.sq,c(8,8),N_site),
      lty=3)
lines(o.sq,FI_fnchypergeo(o.sq,c(16,16),N_site),
      lty=4)
legend('topleft',legend=c('2,2','4,4','8,8','16,16'),
       lty=1:4,cex=0.5,bty='n')

area <- numeric()
for(i in 1:nrow(stan_data$occ)) {
  f <- function(x) 
    sqrt(FI_fnchypergeo(x,as.vector(unlist(stan_data$occ[i,])),N_site))
  area2[i] <- integrate(f,-25,25,subdivisions=1e4)$value
}
area

write.csv(stan_data$occ,'table.csv')

# sensitivity analysis
test.dat <- data.frame(oo=numeric(),   # occupancy pattern: c(2,4,8,16)
                       size = numeric(), # sample size: c(1e1,1e2,1e3)
                       lomega=numeric(), # lomega: -2:2
                       rep = numeric(), # rep number: 1:10
                       grp=numeric(), # group id: 1:600
                       m1=numeric(), # Sp2 present
                       m2=numeric(), # Sp2 absent
                       n=numeric(), # total for Sp1
                       sb=numeric()) # number with both
occ <- 0
grp <- 0
for(i in c(2,4,8,16)) {
  occ <- occ + 1
  for(j in c(1e1,1e2,1e3)) {
    for(k in -2:2) {
      for(l in 1:10) {
        grp <- grp + 1
        rnd <- rFNCHypergeo(j, i, 42 - i, i, exp(k))
        test.dat <- rbind(test.dat,
                          data.frame(oo=occ,size=j,lomega=k,rep=l,
                                     grp=grp,m1=i,m2=42-i,n=i,sb=rnd))
      }
    }
  }
}

stan_code <- 
'
functions {

  // loglikelihood
  real fnchypergeo_lpmf(int[] sb, int[] occ, vector ln_omega, int N_site) {
    vector[num_elements(sb)] ll_num_summand;
    vector[num_elements(sb)] ll_den_summand;
    vector[occ[1] + 1] ll;
    vector[occ[1] + 1] u;
    for(i in 0:occ[1]) {
        u[i+1] = i;
        ll[i+1] = lchoose(occ[1], u[i+1]) + 
                  lchoose(N_site - occ[1], occ[2] - u[i+1]);
    }
    for(i in 1:num_elements(sb)) {
       ll_num_summand[i] = ll[sb[i] + 1] + sb[i]*ln_omega[i];
       ll_den_summand[i] = log_sum_exp(ll + u*ln_omega[i]);
    }    
    return sum(ll_num_summand) - sum(ll_den_summand);
  }

  // jeffreys prior
  real jeffreys_lpdf(vector ln_omega, int[] occ, int N_site, 
                     real jnormalize) {
    vector[num_elements(ln_omega)] jeffreys;
    vector[occ[1] + 1] ll;
    vector[occ[1] + 1] u;
    vector[occ[1] + 1] u2;
    vector[3] l_sum;
    for(i in 0:occ[1]) {
        u[i+1] = i;
        ll[i+1] = lchoose(occ[1], u[i+1]) + 
                  lchoose(N_site - occ[1], occ[2] - u[i+1]);
    }
    u2 = u .* u;
    for(i in 1:num_elements(ln_omega)) {
      l_sum[1] = sum(exp(ll + u*ln_omega[i]));
      l_sum[2] = sum(u .* exp(ll + u*ln_omega[i]));
      l_sum[3] = sum(u2 .* exp(ll + u*ln_omega[i]));
      jeffreys[i] = 0.5*log(l_sum[3]/l_sum[1] - (l_sum[2]/l_sum[1])^2);
    }
    return sum(jeffreys - jnormalize);
  }

}

data {

  int<lower=1> N; // number of pairs
  int<lower=1> N_occ; // number of unique occupancy patterns
  int<lower=1> N_site; // number of sites
  int<lower=1,upper=600> gg[N]; // group assigment for each pair
  int<lower=1,upper=N_occ> oo[N]; // occupancy pattern for each pair
  int<lower=0,upper=N_site> sb[N]; // sites with both species
  int<lower=0,upper=N_site> occ[N_occ,2]; // occupancy table
  int<lower=1,upper=N> oo_l[N_occ]; // pairs index to occ for gg/sb (lower)
  int<lower=1,upper=N> oo_u[N_occ]; // pairs index to occ for gg/sb (upper)
  vector[N_occ] jnormalize; // number of pairs

}

parameters {

  vector[600] ln_omega;

}

transformed parameters {

  vector[N] ln_omega_pr = ln_omega[gg];

}

model {

  // loglikelihood
  for (i in 1:N_occ) {
   // ln_omega ~ jeffreys(occ[i], N_site, jnormalize[i]);
    target += fnchypergeo_lpmf(sb[oo_l[i]:oo_u[i]] | occ[i],
                               ln_omega_pr[oo_l[i]:oo_u[i]],N_site);
  }

}
'
area <- numeric()
for(i in 1:nrow(stan_data$occ)) {
  f <- function(x) 
    sqrt(FI_fnchypergeo(x,as.vector(unlist(stan_data$occ[i,])),N_site))
  area[i] <- integrate(f,-25,25,subdivisions=1e4)$value
}
area

stan_data <- list(N = nrow(test.dat),
                  N_occ = max(test.dat$oo),
                  N_site = 42,
                  gg = test.dat$grp,
                  oo = test.dat$oo,
                  sb = test.dat$sb,
                  occ = unique(test.dat[c('m1','n')]),
                  oo_l = tapply(1:nrow(test.dat),test.dat$oo,min),
                  oo_u = tapply(1:nrow(test.dat),test.dat$oo,max),
                  jnormalize = log(area))

stan_fit3 <- stan(model_code=stan_code,pars=c('ln_omega'),
                 data= stan_data,chains=1)

dog <- summary(stan_fit)$summary

test.dat2 <- unique(test.dat[1:8])
test.dat2$ln_omega <- summary(stan_fit2,'ln_omega')$summary[,1]
test.dat2$bias <- test.dat2$ln_omega - test.dat2$lomega
test.summary <- aggregate(ln_omega~size+lomega+m1,test.dat2,mean)
test.summary$ln_omega_sd <- aggregate(ln_omega~size+lomega+m1,test.dat2,sd)[,4]
test.summary$bias <- aggregate(bias~size+lomega+m1,test.dat2,mean)[,4]
test.summary$bias_se <- aggregate(bias~size+lomega+m1,test.dat2,sd)[,4]/sqrt(10)
test.summary$lomega_plot <- test.summary$lomega + (log10(test.summary$size)-2)/10

pdf('dog.pdf')
  par(mfrow=c(2,2))
  plot.f <-function(temp) {
    plot(temp$lomega_plot,temp$bias,
        xlab=expression(log(omega)),ylab=expression(bias),
        ylim=c(-1,1))
    abline(h=0,col=grey(0.5))
    ltyy <- 0
    for(i in c(10,100,1000)) {
      temp2 <- temp[test.summary$size==i,]
      ltyy <- ltyy + 1
      lines(temp2$lomega_plot,temp2$bias,lty=ltyy)
      for(i in 1:5) {
        segments(temp2$lomega_plot[i],temp2$bias[i]-qt(0.975,9)*temp2$bias_se[i],
                 temp2$lomega_plot[i],temp2$bias[i]+qt(0.975,9)*temp2$bias_se[i])
      }
    }
  }
  plot.f(test.summary[test.summary$m1==2,])
  plot.f(test.summary[test.summary$m1==4,])
  plot.f(test.summary[test.summary$m1==8,])
  plot.f(test.summary[test.summary$m1==16,])
dev.off()

test <- rFNCHypergeo(1000, 1, 41, 1,exp(-2))

ll <- function(x) -sum(fnchypergeo(test,c(1,1),x,41))

optim(0,ll,method='Brent',lower=-10,upper=10)$par


