library(rstan)
library(brms)
library(bayesplot)
library(ggplot2)

ppc_function <- function(stan_fit) {
  xx <- standata(stan_fit)
  theta_vl <- as.matrix(stan_fit, variable = 'theta_vl')
  theta_mu <- as.matrix(stan_fit, variable = 'theta_mu')
  theta_sd <- as.matrix(stan_fit, variable = 'theta_sd')
  rsamp <- sample.int(nrow(theta_mu), size = 20)
  stan_data <- list(N_site = xx$vint4[1],
                    N_dat = length(xx$vint1),
                    nsb = xx$vint3,
                    Y = xx$Y,
                    occ1 = xx$vint1,
                    occ2 = xx$vint2,
                    N_smp = nrow(theta_mu),
                    theta_vl = theta_vl,
                    theta_mu = theta_mu[rsamp,],
                    theta_sd = theta_sd[rsamp,])
  rand_fit <- stan('./stan/nch_ppc.stan', data = stan_data,
                   algorithm = 'Fixed_param', iter = 1, chains = 1)
  return(rand_fit)
}

ppc_plot <- function(stan_ppc, stan_fit) {
  Y <- standata(stan_fit)$Y
  nsb <- standata(stan_fit)$vint3
  y_obs <- numeric()
  for(i in 1:length(Y)) {y_obs <- c(y_obs, rep(Y[i], nsb[i]))}
  y_rep <- as.matrix(extract(stan_ppc,'y_rep')[[1]][1,,])
  return(ppc_rootogram(y_obs, y_rep))
}

files <- c('./Results/bat_winner.rds',
           './Results/bat_winner_NoTrn.rds',
           './Results/bird_winner.rds',
           './Results/bird_winner_NoTrn.rds')
os_quantiles <- list()
ppc_plots <- list()
for(i in 1:4) {
  print(i)
  winner <- readRDS(files[i])
  ppc_fit <- ppc_function(winner)
  chisq_os <- extract(ppc_fit, 'chisq_os')[[1]][1,]
  os_quantiles[[i]] <- quantile(chisq_os, c(0.025, 0.975))
  ppc_plots[[i]] <- ppc_plot(ppc_fit, winner)
}

a <- ppc_plots[[1]] + xlim(c(-0.5,20.5)) + ylim(c(0,125)) + theme(legend.position = "none")
b <- ppc_plots[[2]] + xlim(c(-0.5,20.5)) + ylim(c(0,125)) + theme(legend.position = "none")
c <- ppc_plots[[3]] + ylim(c(0,600)) + theme(legend.position = "none")
d <- ppc_plots[[4]] + ylim(c(0,600)) + theme(legend.position = "none")

pdf('ppc.pdf',width=6,height=6)
plot_grid(a, b, c, d, nrow = 2, labels = LETTERS[1:4])
dev.off()
