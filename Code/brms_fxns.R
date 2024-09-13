
## functions ####

# define custom family for Fisher's Noncentral Hypergeometric distribution
nch <- custom_family(name  = "nch", dpars = "mu", links = "identity", type  = "int",
                     vars  = c("ii[ii_lu[occ_id[n], 1]:ii_lu[occ_id[n], 2]]",
                               "lp[ii_lu[occ_id[n], 1]:ii_lu[occ_id[n], 2]]"))

# produce model data
stan_data_fun <- function(df, nD_cut = 5) {
  df$diet.match[df$diet.match=='Related'] <- 'Intersecting' # related diets excluded
  #df <- df[df$diet.match != "Similar",]
  df$DietGroup <- factor(tolower(df$diet.match))
  df$Habitat <- factor(df$status)
  df$DietPair <- factor(df$diet.pair)
  df$Constant <- factor('Constant') # used for model selection
  df$occ1 <- apply(df[c('presSp1', 'presSp2')], 1, min)
  df$occ2 <- apply(df[c('presSp1', 'presSp2')], 1, max)
  df$N_site <- df$presSp1[1] + df$absSp1[1]
  df$OccSet <- factor(paste(gsub(' ', '0', format(df$occ1)),
                            gsub(' ', '0', format(df$occ2)),
                            gsub(' ', '0', format(df$N_site)), sep = '-'))
  df$sb <- df$presBoth
  df$D <- scale(log(df$dist))[,1]
  df$PhyloD <- cut(df$D, nD_cut, labels=FALSE)
  df$PhyloD <- tapply(df$D, df$PhyloD, mean)[df$PhyloD]
  df <- df[order(df$occ1, df$occ2, df$N_site),]
  
  df2 <- df %>% group_by(occ1, occ2, N_site, OccSet, DietGroup, Habitat, 
                         DietPair, PhyloD, sb, Constant) %>% 
    summarise(N_sb = n())
  return(list(df2, df))
}

# leave-one-out validation
build_loo <- function(stan_fit) {
  # loo calculations take into account the number of pairs represented by each
  # observation by duplicating each log_lik contribution N_sb times; loo
  # calculations are undertaken using the function class because this is far 
  # more memory efficient (though slower)
  dd <- dim(as.array(stan_fit))
  log_lik <- as.data.frame(stan_fit, variable = 'log_lik')
  N_sb <- standata(stan_fit)$weights
  obs_idx <- numeric()
  for(i in seq_along(N_sb)) {
    obs_idx <- c(obs_idx, rep(i, N_sb[i]))
  }
  llfun <- function(data_i, draws = draws) {
    return(draws[, data_i$idx])
  }
  r_eff <- loo::relative_eff(llfun, chain_id = rep(1:dd[2], each = dd[1]), 
                             data = data.frame(idx = obs_idx), 
                             draws = exp(log_lik))
  return(loo::loo(llfun, data = data.frame(idx = obs_idx),
                  draws = log_lik, r_eff = r_eff))
}

# selection of appropriate generated quantities function
nch_stanvars <- function(brm_form) {
  # start: functions to define custom brms family ----
  nch_functions <- stanvar(block = "functions", scode = "
  real nch_lpmf(int sb, real mu, data array[] int ii, data array[] real lp) {
    real log_p = lp[sb - ii[1] + 1] + mu*sb -
                 log_sum_exp(to_vector(lp) + mu*to_vector(ii));
    return (is_inf(log_p) || is_nan(log_p)) ? negative_infinity() : log_p;
  }
  //
  real integrand_re(real re, real xc, array[] real pars, 
                    data array[] real lp, data array[] int x_i) {
      real theta_mu = pars[1];
      array[num_elements(pars) - 1] real re_sd = pars[2:num_elements(pars)];
      int sb = x_i[1]; // observed overlap in occupancy
      array[num_elements(lp)] int ii = x_i[2:num_elements(x_i)];
      real p = exp(nch_lpmf(sb | theta_mu + re, ii, lp) +
                   normal_lpdf(re | 0, sqrt(sum(square(re_sd)))));
      return (is_inf(p) || is_nan(p)) ? 0 : p;
  }
  //
  // integrated log-likelihood
  real log_lik_re(array[] real pars, data array[] real lp, 
                  data array[] int x_i) {
      real p = integrate_1d(integrand_re, negative_infinity(), 
                            positive_infinity(), pars, lp, x_i);
      return (is_inf(p) || is_nan(p)) ? negative_infinity() : log(p);
   }
  //
")
  
  
  nch_tdata <- stanvar(block = "tdata", scode = "
  // the code assumes (and therefore tests) that occ1 and occ2 are in ascending
  // order with occ1 <= occ2; the code yields:
  //    int N_occ: number of occupancy sets (# of occ1/occ2/N_site combinations)
  //    int N_ii: number of elements in ii/lp across all occupancy sets
  //    array[N] int occ_id: occupancy-set assignment for each observation
  //    array[N_occ, 2] int ii_lu: look-up table in ii for each occupancy set
  //    array[N_ii] int ii: allowed co-occurrence values for each occupancy set
  //    array[N_ii] real lp: log(# of permutations) for each element in ii
  //
  int N_occ, N_ii;
  array[N] int occ_id;
  int occ1 = 0, occ2 = 0, N_site = 0;
  int sb_min, sb_max;
  int occ_id_idx = 0, ii_idx = 0;
  for(i in 1:N) {
    if( (vint1[i] < occ1) || ((vint1[i] == occ1) && (vint2[i] < occ2)) ||
        (vint1[i] > vint2[i]) ) {
      reject(\"Error: occ1 & occ2 must be in ascending order w/ occ1<=occ2.\");
    }
    if( (vint1[i] != occ1) || (vint2[i] != occ2) || (vint3[i] != N_site) ) {
      occ_id_idx += 1;
      occ1 = vint1[i];
      occ2 = vint2[i];
      N_site = vint3[i];
      sb_min = occ1 + occ2 > N_site ? occ1 + occ2 - N_site : 0;
      sb_max = occ1;
      ii_idx += sb_max - sb_min + 1;
    }
    occ_id[i] = occ_id_idx;
  }
  N_occ = occ_id_idx;
  N_ii = ii_idx;
  //
  array[N_occ, 2] int ii_lu;
  array[N_ii] int ii;
  array[N_ii] real lp;
  occ1 = 0; occ2 = 0; N_site = 0;
  occ_id_idx = 0; ii_idx = 0;
  for(i in 1:N) {
    if( (vint1[i] != occ1) || (vint2[i] != occ2) || (vint3[i] != N_site) ) {
      occ_id_idx += 1;
      occ1 = vint1[i];
      occ2 = vint2[i];
      N_site = vint3[i];
      sb_min = occ1 + occ2 > N_site ? occ1 + occ2 - N_site : 0;
      sb_max = occ1;
      ii_lu[occ_id_idx, 1] = (occ_id_idx == 1) ? 1 : ii_lu[occ_id_idx - 1, 2] + 1;
      ii_lu[occ_id_idx, 2] = ii_lu[occ_id_idx, 1] + sb_max - sb_min;
      for(n in sb_min:sb_max) {
        ii_idx += 1;
        ii[ii_idx] = n;
        lp[ii_idx] = lchoose(occ1, n) + lchoose(N_site - occ1, occ2 - n);
      }
    }
  }
  //
")
  
  # end: functions to define custom brms family ----
  x <- deparse1(brm_form) # convert formula to text
  N_re <- lengths(regmatches(x, gregexpr("gr", x))) # number of random effects
  intercept_only <- lengths(regmatches(x, gregexpr("~ 1", x)))
  nch_stanvar <- nch_functions + nch_tdata
  
  if(intercept_only == 1) { 
    nch_stanvar <- nch_stanvar + 
      stanvar(block ="genquant", scode =
                "vector[N] theta = rep_vector(Intercept, N);")
  } else {
    nch_stanvar <- nch_stanvar + 
      stanvar(block="genquant", scode =
                "vector[N] theta = Intercept + Xc * b;")
  }
  
  if(N_re > 0) {
    nch_stanvar <- nch_stanvar + 
      stanvar(block="genquant", scode =
                "vector[N] theta_mu = theta;")
  }
  
  nch_stanvar <- nch_stanvar + 
    stanvar(block="genquant", scode =
              "vector[N] log_lik;\n\n  for (n in 1:N) {")
  
  if(N_re >= 1) {
    nch_stanvar <- nch_stanvar +
      stanvar(block = "genquant", scode =
                "theta[n] += r_1_1[J_1[n]] * Z_1_1[n];")
  }
  
  if(N_re == 2) {
    nch_stanvar <- nch_stanvar +
      stanvar(block = "genquant", scode =
                "theta[n] += r_2_1[J_2[n]] * Z_2_1[n];")
  }
  
  if(N_re == 0) {
    nch_stanvar <- nch_stanvar +
      stanvar(block = "genquant", scode =
                "log_lik[n] =
      nch_lpmf(Y[n] | theta[n], ii[ii_lu[occ_id[n], 1]:ii_lu[occ_id[n], 2]],
      lp[ii_lu[occ_id[n], 1]:ii_lu[occ_id[n], 2]]);\n}")
  } else if(N_re == 1) {
    nch_stanvar <- nch_stanvar +
      stanvar(block = "genquant", scode =
                "log_lik[n] = log_lik_re(append_array({theta_mu[n]},
      {sd_1[1, Jby_1[J_1[n]]]}),
      lp[ii_lu[occ_id[n], 1]:ii_lu[occ_id[n], 2]],
      append_array({Y[n]}, ii[ii_lu[occ_id[n], 1]:ii_lu[occ_id[n], 2]]));\n}")
  } else {
    nch_stanvar <- nch_stanvar +
      stanvar(block = "genquant", scode =
                "log_lik[n] = log_lik_re(append_array({theta_mu[n]},
      append_array({sd_1[1, Jby_1[J_1[n]]]}, {sd_2[1, Jby_2[J_2[n]]]})),
      lp[ii_lu[occ_id[n], 1]:ii_lu[occ_id[n], 2]],
      append_array({Y[n]}, ii[ii_lu[occ_id[n], 1]:ii_lu[occ_id[n], 2]]));\n}")
  }
  return(nch_stanvar)
}

# generate model formulas and run random then fixed effects models.
find_best_model <- function(tax, standata){
  
  ## fixed and random effects structures ----
  re_steps <- character()
  for(i in c("(1|gr(Habitat:DietPair, by = Constant))",
             "(1|gr(DietPair, by = Constant))",
             "")) {
    for(j in c("(1|gr(DietPair:PhyloD:Habitat:OccSet, by = Constant))",
               "(1|gr(DietGroup:PhyloD:Habitat:OccSet, by = Constant))",
               "(1|gr(DietPair:Habitat:OccSet, by = Constant))",
               "(1|gr(DietGroup:Habitat:OccSet, by = Constant))",
               "(1|gr(Habitat:OccSet, by = Constant))","")) {
      re_steps <- c(re_steps, paste0(c(i, j),collapse='+')) %>% str_replace("\\+$", "") %>% 
        str_replace("^\\+", "")
    }
  }
  re_steps[length(re_steps)] <- "1"
  f00 <- function(x) {sub(' ','0', format(c(99,x)))[-1]}
  names(re_steps) <- paste0('re', f00(seq_along(re_steps)))
  
  fx_steps <- c("PhyloD + DietGroup + Habitat + PhyloD:DietGroup + PhyloD:Habitat + DietGroup:Habitat + PhyloD:DietGroup:Habitat",
                "PhyloD + DietGroup + Habitat + PhyloD:DietGroup + PhyloD:Habitat + DietGroup:Habitat",
                "PhyloD + DietGroup + Habitat + PhyloD:Habitat + DietGroup:Habitat",
                "PhyloD + DietGroup + Habitat + PhyloD:DietGroup + DietGroup:Habitat",
                "PhyloD + DietGroup + Habitat + DietGroup:Habitat",
                "DietGroup + Habitat + DietGroup:Habitat",
                "PhyloD + DietGroup + Habitat + PhyloD:DietGroup + PhyloD:Habitat",
                "PhyloD + DietGroup + Habitat + PhyloD:Habitat",
                "PhyloD + Habitat + PhyloD:Habitat",
                "PhyloD + DietGroup + Habitat + PhyloD:DietGroup",
                "PhyloD + DietGroup + PhyloD:DietGroup",
                "PhyloD + DietGroup + Habitat",
                "DietGroup + Habitat",
                "PhyloD + Habitat",
                "Habitat",
                "PhyloD + DietGroup",
                "DietGroup",
                "PhyloD",
                "1")
  names(fx_steps) <- paste0('fx', f00(seq_along(fx_steps)))
  ### end define fixed and random effects structures ----
  
  # random effects models
  re_warnings <- list()
  re_loo <- list()
  for(i in seq_along(re_steps)) {
    brm_form <- formula(paste(paste0('sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ ', 
                                     fx_steps[1]), re_steps[i], sep= ' + '))
    print(paste0('model ', i, ' of ', length(re_steps), ':'))
    print(brm_form)
    stan_fit <- brm(brm_form, family = nch, init = 0, data = standata,
                    iter = 3e3, warmup = 1e3, thin = 2, cores = 4, 
                    stanvars = nch_stanvars(brm_form))
    saveRDS(stan_fit, paste0('./stan/', tax ,'_fx01_re', f00(i), '.rds'))
    re_warnings[[paste0(tax, '_fx01_re', f00(i))]] <- warnings()
    re_loo[[paste0(tax, '_fx01_re', f00(i))]] <- build_loo(stan_fit)
  }
  #return(list(standata, re_warnings, re_loo))
  #
  #loo_re <- map2(re_formulas, names(re_formulas), ## saves models and returns elpd evaluations
  #              ~fit_brms_model(formula = as.formula(.x), data = standata, modelname = paste0("fixd01_", .y), tax = tax))
  
  # Select best random effects structure 
  best_randeff <- (loo_compare(re_loo) %>% data.frame() %>% rownames())[1] %>% word(3,3,sep = "_") # get best random effects structure
  message("Best random effects structure is ", best_randeff, ": ", re_steps[[best_randeff]])
  saveRDS(re_loo, paste0("./stan/", tax, "_loo_randeff.rds"))
  
  
  # fixed effects models
  fx_warnings <- list()
  fx_loo <- list()
  for(i in seq_along(fx_steps)[-1]) {
    brm_form <-
      formula(paste(paste0('sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ ',
                           fx_steps[i]), re_steps[best_randeff], sep= ' + '))
    print(paste0('model ', i, ' of ', length(fx_steps), ':'))
    print(brm_form)
    stan_fit <- brm(brm_form, family = nch, init = 0, data = standata, 
                    iter = 3e3, warmup = 1e3, thin = 2, cores = 4,
                    stanvars = nch_stanvars(brm_form))
    saveRDS(stan_fit, paste0('./stan/', tax, '_fx', f00(i), '_', best_randeff,'.rds'))
    fx_warnings[[paste0(tax, '_fx01_re', f00(i))]] <- warnings()
    fx_loo[[paste0(tax, '_fx', f00(i), "_", best_randeff)]] <- build_loo(stan_fit)
  }
  
    # 
  # # fixed effects formulas
  # fx_formulas <- paste0('sb | vint(occ1, occ2, N_sb, N_site) ~ ', re_steps[[best_randeff]]) 
  # fx_formulas <-  map(fx_steps, ~paste(c(fx_formulas, .x), collapse = "+")) %>% 
  #   str_replace("\\+$", "")
  # names(fx_formulas) <- paste0('fixd',gsub(' ','0', format(seq_along(fx_steps)-1)), "_", best_randeff)
  # if(best_randeff %in% c("rand0", "rand00")){fx_formulas <- fx_formulas[-1]}
  # 
  # loo_fx <- map2(fx_formulas, names(fx_formulas), ## saves models and returns elpd evaluations
  #                ~fit_brms_model(formula = as.formula(.x), data = standata, modelname = .y, tax = tax))
  best_model <- (loo_compare(fx_loo) %>% data.frame() %>% rownames())[1] # get best random effects structure
  message("Best fixed effects structure is ", best_model, ": ", 
          paste(paste0('sb | vint(occ1, occ2, N_site) + weights(N_sb) ~ ',
                    fx_steps[word(best_model,2,2,sep = "_")]), re_steps[best_randeff], sep= ' + '))
  
  saveRDS(fx_loo, paste0("./stan/", tax, "_loo_fixdeff.rds"))
  #load best model into environment and examine summary
  
  winner <- readRDS(paste0("./stan/", best_model, ".rds"))
  return(winner)
}





