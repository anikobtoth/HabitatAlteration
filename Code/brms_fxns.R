
## functions ####
# produce model data
stan_data_fun <- function(df, medn, nD_cut = 5) {
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
  df <- df %>% filter(diet.match != "Similar")
  df$diet.match <- factor(tolower(df$diet.match))
  df$status <- factor(df$status)
  df$diet.pair <- factor(df$diet.pair)
  df$fixed <- factor('fixed') # used for model selection
  df$occ1 <- apply(df[c('presSp1','presSp2')],1,min)
  df$occ2 <- apply(df[c('presSp1','presSp2')],1,max)
  df$oc01f <- factor(interaction(gsub(' ','0',format(df$occ1)),
                                 gsub(' ','0',format(df$occ2))))
  df$oc02f <- factor(interaction(gsub(' ','0',format(ceiling(df$occ1/2)*2)),
                                 gsub(' ','0',format(ceiling(df$occ2/2)*2))))
  df$oc03f <- factor(interaction(gsub(' ','0',format(ceiling(df$occ1/3)*3)),
                                 gsub(' ','0',format(ceiling(df$occ2/3)*3))))
  df$ocl1.2f <- 
    factor(interaction(gsub(' ','0',
                            format(floor(1.2^ceiling(log(df$occ1,base=1.2))))),
                       gsub(' ','0',
                            format(floor(1.2^ceiling(log(df$occ2,base=1.2)))))))
  df$ocl1.4f <- 
    factor(interaction(gsub(' ','0',
                            format(floor(1.4^ceiling(log(df$occ1,base=1.4))))),
                       gsub(' ','0',
                            format(floor(1.4^ceiling(log(df$occ2,base=1.4)))))))
  df$ocl1.6f <- 
    factor(interaction(gsub(' ','0',
                            format(floor(1.6^ceiling(log(df$occ1,base=1.6))))),
                       gsub(' ','0',
                            format(floor(1.6^ceiling(log(df$occ2,base=1.6)))))))
  df$sb <- df$presBoth
  df$D <- scale(log(df$dist))[,1]
  df$D_cut <- cut(df$D, nD_cut, labels=FALSE)
  df$D_cut <- tapply(df$D,df$D_cut, mean)[df$D_cut]
  
  df2 <- df %>% group_by(occ1, occ2, oc01f, oc02f, oc03f, ocl1.2f, ocl1.4f, ocl1.6f, diet.match, status, diet.pair, D_cut, sb) %>% 
    summarise(N_sb = n()) %>% mutate(N_site = df$presSp1[1] + df$absSp1[1], 
                                     fixed = factor("fixed"))
  
  return(list(df2,df))
}

# leave-one-out validation
build_loo <- function(x, nsb, ...) {
  # loo calculations take into account the number of pairs represented by each
  # observation by duplicating each log_lik contribution N_sb times; loo
  # calculations are undertaken using the function class because this is far 
  # more memory efficient (though slower)
  dd <- dim(as.array(x))
  log_lik <- as.data.frame(x)
  log_lik <- log_lik[grep('log_lik_gen',names(log_lik))]
  obs_idx <- numeric()
  for(i in seq_along(nsb)) {
    obs_idx <- c(obs_idx, rep(i, nsb[i]))
  }
  llfun <- function(data_i, draws = draws) {
    return(draws[,data_i$idx])
  }
  r_eff <- loo::relative_eff(llfun, chain_id = rep(1:dd[2], each = dd[1]), 
                             data = data.frame(idx = obs_idx), 
                             draws = exp(log_lik))
  return(loo::loo(llfun, data = data.frame(idx = obs_idx),
                  draws = log_lik, r_eff = r_eff))
}

# selection of appropriate generated quantities function
select_genquant <- function(formula) {
  n_rand <- gregexpr("gr", formula) 
  n_rand <- ifelse(n_rand[1] == -1, 0, length(n_rand))
  
  n_fixd <- length(gregexpr("+", formula)) + 1 - n_rand
  
  if(n_fixd > 0) {
    # used to calculate log-likelihood when there are no random effects
    if(n_rand == 0) {genquant <- stanvar(block = "genquant", scode =
                                           "//
        vector[N] theta_vl = Intercept + Xc * b;
        vector[N] log_lik_gen;
        //
        for (i in 1:N) {
          log_lik_gen[i] = 
            nch_ll_lpmf(Y[i] | theta_vl[i], ll[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]],
                        ii[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]]);
         }
        //
        ")}
    # used to calculate integrated log-likelihood when there is 1 random effect
    # and fixed effects
    if(n_rand == 1) {genquant <- stanvar(block = "genquant", scode =
                                           "//
          vector[N] theta_mu = Intercept + Xc * b;
          vector[N] theta_vl = theta_mu;
          vector[N] theta_sd;
          vector[N] log_lik_gen;
          //
          for (i in 1:N) {
            theta_vl[i] += r_1_1[J_1[i]] * Z_1_1[i];
            theta_sd[i] = sd_1[1, Jby_1[J_1[i]]];
            log_lik_gen[i] = 
                log_lik_int1D(append_array({theta_mu[i]}, {theta_sd[i]}),
                              ll[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]], 
                              append_array({Y[i]},
                                ii[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]]));
           }
          //
          ")}
    # used to calculate integrated log-likelihood when there are 2 random effects
    if(n_rand == 2) {genquant <- stanvar(block = "genquant", scode =
                                           "//
              vector[N] theta_mu = Intercept + Xc * b;
              vector[N] theta_vl = theta_mu;
              vector[N] theta_sd1;
              vector[N] theta_sd2;
              vector[N] log_lik_gen;
              //
              for (i in 1:N) {
                theta_vl[i] += r_1_1[J_1[i]] * Z_1_1[i] + r_2_1[J_2[i]] * Z_2_1[i];
                theta_sd1[i] = sd_1[1, Jby_1[J_1[i]]];
                theta_sd2[i] = sd_2[1, Jby_2[J_2[i]]];
                log_lik_gen[i] = 
                  log_lik_int2D(append_array({theta_mu[i]}, 
                                     append_array({theta_sd1[i]}, {theta_sd2[i]})),
                                ll[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]],
                                append_array({Y[i]},
                                      ii[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]]));
               }
              //
              ")}
  }else{
    # used to calculate integrated log-likelihood when there is 1 random effect
    # and no fixed effects
    if(n_rand == 1) {genquant <- stanvar(block = "genquant", scode =
                                           "//
              vector[N] theta_mu = rep_vector(Intercept, N);
              vector[N] theta_vl = theta_mu;
              vector[N] theta_sd;
              vector[N] log_lik_gen;
              //
              for (i in 1:N) {
                theta_vl[i] += r_1_1[J_1[i]] * Z_1_1[i];
                theta_sd[i] = sd_1[1, Jby_1[J_1[i]]];
                log_lik_gen[i] = 
                    log_lik_int1D(append_array({theta_mu[i]}, {theta_sd[i]}),
                                  ll[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]], 
                                  append_array({Y[i]},
                                    ii[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]]));
               }
              //
              ")}
    # used to calculate log likelihood when there are 2 random effects
    # just for intercept only model
    if(n_rand == 2) {genquant <- stanvar(block = "genquant", scode =
                                           "//
            vector[N] theta_mu = rep_vector(Intercept, N);
            vector[N] theta_vl = theta_mu;
            vector[N] theta_sd1;
            vector[N] theta_sd2;
            vector[N] log_lik_gen;
            //
            for (i in 1:N) {
              theta_vl[i] += r_1_1[J_1[i]] * Z_1_1[i] + r_2_1[J_2[i]] * Z_2_1[i];
              theta_sd1[i] = sd_1[1, Jby_1[J_1[i]]];
              theta_sd2[i] = sd_2[1, Jby_2[J_2[i]]];
              log_lik_gen[i] = 
                log_lik_int2D(append_array({theta_mu[i]}, 
                                   append_array({theta_sd1[i]}, {theta_sd2[i]})),
                              ll[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]],
                              append_array({Y[i]},
                                    ii[ll_idx[occ_id[i], 1]:ll_idx[occ_id[i], 2]]));
             }
            //
            ")}
  }
  return(genquant)
}

# fit a single model, save it to file  and return loo validation object
fit_brms_model <- function(formula, data, modelname, tax, c = 4){
  # start: functions to define custom brms family ----
  # vint1 and vint2: numbers of sites occupied by pairs (occ1, occ2)
  # vint3: number of pairs represented by observation (N_sb)
  # vint4: number of sites (N_site); MUST BE IDENTICAL FOR ALL PAIRS IN ANALYSIS
  nch <- custom_family('nch', dpars = c('mu'), links = c('identity'), 
                       type = 'int',
                       vars = c('vint1[n]','vint2[n]','vint3[n]','vint4[n]'))
  
  nch_functions <- stanvar(block = 'functions', scode = ######
                             "// log-likehood function used in model block
real nch_lpmf(int sb, real mu, int occ1, int occ2, int N_sb, int N_site) {
    int sb_min = (occ1 + occ2) > N_site ? (occ1 + occ2) - N_site : 0;
    int sb_max = occ2 > occ1 ? occ1 : occ2;
    vector[sb_max - sb_min + 1] ii;
    vector[sb_max - sb_min + 1] ll;
    for(i in sb_min:sb_max) {
      ii[i - sb_min + 1] = i;
      ll[i - sb_min + 1] = lchoose(occ1, i) + 
                           lchoose(N_site - occ1, occ2 - i);
    }
    return (ll[sb - sb_min + 1] + mu*sb - log_sum_exp(ll + mu*ii))*N_sb;
 }
//
// log-likelihood function used in generated quantities block;
// it works faster than nch_lpmf because it uses log(# permutations) values (ll)
// that have been precalculated in the transformed data block
real nch_ll_lpmf(int sb, real theta, data array[] real ll, 
                 data array[] int ii) {
    real log_p = ll[sb - ii[1] + 1] + theta*sb -
           log_sum_exp(to_vector(ll) + theta*to_vector(ii));
    return (is_inf(log_p) || is_nan(log_p)) ? negative_infinity() : log_p;
 }
//
// integrand used to calculate likelihood for models with 1 random effect
real integrand_1D(real re, real xc, array[] real pars, 
                  data array[] real ll, data array[] int x_i) {
  real theta_mu = pars[1];
  real theta_sd = pars[2];
  int sb = x_i[1];
  array[num_elements(ll)] int ii = x_i[2:num_elements(x_i)];
  real p = exp(normal_lpdf(re | 0, theta_sd) +
                 nch_ll_lpmf(sb | theta_mu + re, ll, ii));
  return (is_inf(p) || is_nan(p)) ? 0 : p;
 }
//
// integrated log-likelihood for models with 1 random effect
real log_lik_int1D(array[] real pars, data array[] real ll, 
                   data array[] int x_i) {
  real p = integrate_1d(integrand_1D, negative_infinity(), 
                        positive_infinity(), pars, ll, x_i);
  return (is_inf(p) || is_nan(p)) ? negative_infinity() : log(p);
 }
//
// integrand used to calculate likelihood for models with 2 random effects;
// double integral is reduced to a single integral by noting that the
// convolution of two independent, normally distributed random variables is
// normally distributed with a variance equal to the sum of the two variances
real integrand_2D(real re1_plus_re2, real xc, array[] real pars, 
                         data array[] real ll, data array[] int x_i) {
    real theta_mu = pars[1];
    real theta_sd1 = pars[2];
    real theta_sd2 = pars[3];
    int sb = x_i[1]; // observed overlap in occupancy
    array[num_elements(ll)] int ii = x_i[2:num_elements(x_i)];
    real p = exp(nch_ll_lpmf(sb | theta_mu + re1_plus_re2, ll, ii) +
             normal_lpdf(re1_plus_re2 | 0, sqrt(theta_sd1^2 + theta_sd2^2)));
    return (is_inf(p) || is_nan(p)) ? 0 : p;
 }
//
// integrated likelihood for models with 2 random effects
real log_lik_int2D(array[] real pars1, data array[] real ll, 
                   data array[] int x_i) {
    real p = integrate_1d(integrand_2D, negative_infinity(), 
                          positive_infinity(), pars1, ll, x_i);
    return (is_inf(p) || is_nan(p)) ? negative_infinity() : log(p);
 }
//
")
  
 
  nch_tdata <- stanvar(block = "tdata", scode =  #####
                         "// Data is transformed in the transformed data block are used to speed up the
// calculations of the integrated log-likelihood in generated quantities block;
// Key transformed items include:
// N_occ: number of unique occupancy patterns {min(occ1, occ2), max(occ1, occ2)}
// occ_id: occupancy pattern ID assigned to each observation
// ii: vector representing all possible overlapping values for each occupancy
//     pattern given occ1, occ2, and N_site
// ll: vector representing log(# of permutations) for each element of ii
// ll_idx: index to ii and ll for each occupancy pattern ID
//
// EVERY observation in the analysis MUST have the same number of sites
int N_site = vint4[1];
//
// calculate total # of occupancy patterns (N_occ) and assign each observation
// to an occupancy pattern ID (occ_id)
array[max(append_array(vint1, vint2)), 
      max(append_array(vint1, vint2))] int flag = 
          rep_array(0, max(append_array(vint1, vint2)),
                       max(append_array(vint1, vint2)));
int occ_min;
int occ_max;
int idx = 0;
array[N] int occ_id; // occupancy pattern ID for each observation
for(i in 1:N) {
  occ_min = vint1[i] < vint2[i] ? vint1[i] : vint2[i];
  occ_max = vint1[i] > vint2[i] ? vint1[i] : vint2[i];
  if(flag[occ_min, occ_max] == 0) {
    idx += 1;
    flag[occ_min, occ_max] = idx;
  }
  occ_id[i] = flag[occ_min, occ_max];
 }
int N_occ = idx;
//
// calculate quantities for occupancy sets
array[N_occ] int occ1; // minimum occupancy for pair
array[N_occ] int occ2; // maximum occupancy for pair
array[N_occ] int sb_min; // minimum possible overlap given N_site, occ1, occ2
array[N_occ] int sb_max; // maximum possible overlap given N_site, occ1, occ2
array[N_occ, 2] int ll_idx; // precalculated values for log(# permutations)
flag = rep_array(0, max(append_array(vint1, vint2)),
                    max(append_array(vint1, vint2)));
idx = 0;
for(i in 1:N) {
  occ_min = vint1[i] < vint2[i] ? vint1[i] : vint2[i];
  occ_max = vint1[i] > vint2[i] ? vint1[i] : vint2[i];
  if(flag[occ_min, occ_max] == 0) {
    idx += 1;
    flag[occ_min, occ_max] = idx;
    occ1[idx] = occ_min;
    occ2[idx] = occ_max;
    sb_min[idx] = occ_min + occ_max > N_site ? occ_min + occ_max - N_site : 0;
    sb_max[idx] = occ_min;
    ll_idx[idx, 1] = (idx == 1) ? 1 : ll_idx[idx - 1, 2] + 1;
    ll_idx[idx, 2] = ll_idx[idx, 1] + sb_max[idx] - sb_min[idx];
  }
 }
//
// calculate numbers of permutations for each occupancy set
array[ll_idx[N_occ, 2]] real ll; // log(# permutations) values
array[ll_idx[N_occ, 2]] int  ii; // corresponding sb values
idx = 0;
for(i in 1:N_occ) {
  for(j in sb_min[i]:sb_max[i]) {
    idx += 1;
    ll[idx] = lchoose(occ1[i], j) + 
              lchoose(N_site - occ1[i], occ2[i] - j);
    ii[idx] = j;
  }
 }
//
")
  
  
  genquant_function <- select_genquant(formula)
  
  # end: functions to define custom brms family ----
  
  print(as.formula(formula))
  
  stan_fit <- brm(as.formula(formula), family = nch, data = data, init = 0, cores = c,
        stanvars = nch_functions + nch_tdata + genquant_function)
  
  saveRDS(stan_fit, paste0('./stan/', tax, "_", modelname,'.rds'))
  
  return(build_loo(stan_fit, data$N_sb))
}

# generate model formulas and run random then fixed effects models.
find_best_model <- function(tax, standata){
  
  ## fixed and random effects structures ----
  re_steps <-c("",
               "(1|gr(status:diet.pair, by = status)) + (1|gr(status:diet.pair:oc01f, by = status:diet.match))",
               "(1|gr(status:diet.pair, by = fixed)) + (1|gr(status:diet.pair:oc01f, by = status:diet.match))",
               "(1|gr(diet.pair, by = fixed)) + (1|gr(status:diet.pair:oc01f, by = status:diet.match))",
               "(1|gr(status:diet.pair:oc01f, by = status:diet.match))",
               "(1|gr(status:diet.pair, by = status)) + (1|gr(status:diet.pair:oc01f, by = status))",
               "(1|gr(status:diet.pair, by = fixed)) + (1|gr(status:diet.pair:oc01f, by = status))",
               "(1|gr(diet.pair, by = fixed)) + (1|gr(status:diet.pair:oc01f, by = status))",
               "(1|gr(status:diet.pair:oc01f, by = status))",
               "(1|gr(status:diet.pair, by = status)) + (1|gr(status:diet.pair:oc01f, by = fixed))",
               "(1|gr(status:diet.pair, by = fixed)) +  (1|gr(status:diet.pair:oc01f, by = fixed))",
               "(1|gr(diet.pair, by = fixed)) +  (1|gr(status:diet.pair:oc01f, by = fixed))",
               "(1|gr(status:diet.pair:oc01f, by = fixed))",
               "(1|gr(status:diet.pair, by = status)) + (1|gr(status:oc01f, by = status))",
               "(1|gr(status:diet.pair, by = fixed)) + (1|gr(status:oc01f, by = status))",
               "(1|gr(diet.pair, by = fixed)) + (1|gr(status:oc01f, by = status))",
               "(1|gr(status:oc01f, by = status))",
               "(1|gr(status:diet.pair, by = status)) + (1|gr(status:oc01f, by = fixed))",                 
               "(1|gr(status:diet.pair, by = fixed)) + (1|gr(status:oc01f, by = fixed))",
               "(1|gr(diet.pair, by = fixed)) + (1|gr(status:oc01f, by = fixed))",
               "(1|gr(status:oc01f, by = fixed))",
               "(1|gr(status:diet.pair, by = status))",                                                   
               "(1|gr(status:diet.pair, by = fixed))",                                                    
               "(1|gr(diet.pair, by = fixed))")
  names(re_steps) <- paste0('rand',gsub(' ','0', format(seq_along(re_steps)-1))) 
  # random effects formulas
  re_formulas <- 'sb | vint(occ1, occ2, N_sb, N_site) ~ status*diet.match*D_cut'
  re_formulas <-  map(re_steps, ~paste(c(re_formulas, .x), collapse = "+")) %>% 
    str_replace("\\+$", "")
  names(re_formulas) <- paste0('rand',gsub(' ','0', format(seq_along(re_steps)-1))) 
  
  
  fx_steps <- c("",
                "D_cut + diet.match + status + D_cut:diet.match + D_cut:status + diet.match:status + D_cut:diet.match:status",
                "D_cut + diet.match + status + D_cut:diet.match + D_cut:status + diet.match:status",
                "D_cut + diet.match + status + D_cut:status + diet.match:status",
                "D_cut + diet.match + status + D_cut:diet.match + diet.match:status",
                "D_cut + diet.match + status + diet.match:status",
                "diet.match + status + diet.match:status",
                "D_cut + diet.match + status + D_cut:diet.match + D_cut:status",
                "D_cut + diet.match + status + D_cut:status",
                "D_cut + status + D_cut:status",
                "D_cut + diet.match + status + D_cut:diet.match",
                "D_cut + diet.match + D_cut:diet.match",
                "D_cut + diet.match + status",
                "diet.match + status",
                "D_cut + status",
                "status",
                "D_cut + diet.match",
                "diet.match",
                "D_cut")
  
  ## Code above generates re_steps, re_formulas, and fx_steps----
  
  loo_re <- map2(re_formulas, names(re_formulas), ## saves models and returns elpd evaluations
                 ~fit_brms_model(formula = as.formula(.x), data = standata, modelname = paste0("fixd01_", .y), tax = tax))
  best_randeff <- (loo_compare(loo_re) %>% data.frame() %>% rownames())[1]  # get best random effects structure
  message("Best random effects structure is ", best_randeff, ": ", re_steps[[best_randeff]])
  
  saveRDS(loo_re, paste0("./", tax, "_loo_randeff.rds"))
  # fixed effects formulas
  fx_formulas <- paste0('sb | vint(occ1, occ2, N_sb, N_site) ~ ', re_steps[[best_randeff]]) 
  fx_formulas <-  map(fx_steps, ~paste(c(fx_formulas, .x), collapse = "+")) %>% 
    str_replace("\\+$", "")
  names(fx_formulas) <- paste0('fixd',gsub(' ','0', format(seq_along(fx_steps)-1)), "_", best_randeff)
  if(best_randeff %in% c("rand0", "rand00")){fx_formulas <- fx_formulas[-1]}
  
  loo_fx <- map2(fx_formulas, names(fx_formulas), ## saves models and returns elpd evaluations
                 ~fit_brms_model(formula = as.formula(.x), data = standata, modelname = .y, tax = tax))
  best_model <- (loo_compare(loo_fx) %>% data.frame() %>% rownames())[1] # get best random effects structure
  message("Best fixed effects structure is ", best_model, ": ", fx_formulas[[best_model]])
  
  saveRDS(loo_fx, paste0("./", tax, "_loo_fixdeff.rds"))
  #load best model into environment and examine summary
  
  winner <- readRDS(paste0("./stan/", paste(tax, best_model, sep = "_"), ".rds"))
  return(winner)
}





