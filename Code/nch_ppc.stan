functions {
  //
  // PMF for noncentral hypergeometric distribution (NCH)
  real nch_lpmf(int sb, array[] real ll, array[] int ii, real theta) {
    return ll[sb - ii[1] + 1] + theta*sb -
           log_sum_exp(to_vector(ll) + theta*to_vector(ii));
  }
  //
  // mean of PMF
  real nch_m(array[] real ll, array[] int ii, real theta) {
    vector[num_elements(ii)] summand;
    for(i in 1:num_elements(ii)) {
      summand[i] = log(ii[i]) + nch_lpmf(ii[i] | ll, ii, theta);
    }
    return exp(log_sum_exp(summand));
  }
  //
  // variance of PMF
  real nch_v(array[] real ll, array[] int ii, real theta) {
    vector[num_elements(ii)] summand;
    real mu = nch_m(ll, ii, theta);
    for(i in 1:num_elements(ii)) {
      summand[i] = log((ii[i] - mu)^2) + nch_lpmf(ii[i] | ll, ii, theta);
    }
    return exp(log_sum_exp(summand));
  }
  //
  // generate random draws from NCH and compute sum of squared pearson residuals
  real pearson_rng(array[] real ll, array[] int ii, real theta, int nvals) {
    real theta_mean = nch_m(ll, ii, theta);
    real theta_var = nch_v(ll, ii, theta);
    vector[nvals] chisq_sm;
    int n_ii = num_elements(ii);
    vector[n_ii] log_lik;
    array[nvals] real uvals = uniform_rng(rep_vector(0, nvals), 
                                          rep_vector(1, nvals));
    for(i in 1:n_ii) {
        log_lik[i] = nch_lpmf(ii[i] | ll, ii, theta);
    }
    vector[n_ii] cdf = cumulative_sum(exp(log_lik));
    array[nvals] int rvals;
    for(i in 1:nvals) {
      for(j in 1:n_ii) {
        if(cdf[j] > uvals[i]) {
          rvals[i] = ii[j];
          break;
        }
      }
      chisq_sm[i] = (rvals[i] - theta_mean)^2/theta_var;
    }
    return sum(chisq_sm);
  }
  //
  // random number generator for NCH with random effects
  array[] int nch_re_rng(array[] real ll, array[] int ii, real theta_mu, 
                   real theta_sd, int nvals) {
    int n_ii = num_elements(ii);
    vector[n_ii] log_lik;
    vector[n_ii] cdf;
    array[nvals] real uvals = uniform_rng(rep_vector(0, nvals),
                                          rep_vector(1, nvals));
    array[nvals] real theta = normal_rng(rep_vector(theta_mu, nvals), 
                                         rep_vector(theta_sd, nvals));
    array[nvals] int rvals;
    for(i in 1:nvals) {
      for(j in 1:n_ii) {
        log_lik[j] = nch_lpmf(ii[j] | ll, ii, theta[i]);
      }
      cdf = cumulative_sum(exp(log_lik));
      for(j in 1:n_ii) {
        if(cdf[j] > uvals[i]) {
          rvals[i] = ii[j];
          break;
        }
      }
    }
    return rvals;
  }
}
//
//
data {
  //
  int N_site; // # of sites
  int N_dat; // # of binned observations (dat)
  int nsb[N_dat]; // # of pairwise comparisons per binned observation
  int Y[N_dat]; // observed number of co-occurrences
  int occ1[N_dat]; // occupancy of species 1
  int occ2[N_dat]; // occupancy of species 2
  int N_smp; // # of posterior samples
  matrix[N_smp, N_dat] theta_vl; // theta value
  matrix[20, N_dat] theta_mu; // theta mean
  matrix[20, N_dat] theta_sd; // theta sd
}
//
//
transformed data {
  // calculate total # of occupancy pairings (N_occ) and assign each
  // observation to an occupancy pairing ID (occ_id)
  array[max(append_array(occ1, occ2)), 
        max(append_array(occ1, occ2))] int flag = 
            rep_array(0, max(append_array(occ1, occ2)),
                         max(append_array(occ1, occ2)));
  int occ_min;
  int occ_max;
  int idx = 0;
  array[N_dat] int occ_id; // occupancy pattern ID for each observation
  for(i in 1:N_dat) {
    occ_min = occ1[i] < occ2[i] ? occ1[i] : occ2[i];
    occ_max = occ1[i] > occ2[i] ? occ1[i] : occ2[i];
    if(flag[occ_min, occ_max] == 0) {
      idx += 1;
      flag[occ_min, occ_max] = idx;
    }
    occ_id[i] = flag[occ_min, occ_max];
  }
  int N_occ = idx;
  //
  // calculate quantities for occupancy sets
  array[N_occ] int sb_min; // minimum possible # of co-occurrences
  array[N_occ] int sb_max; // maximum possible # of co-occurrences
  array[N_occ] int occ_s1; // occurrence for occupancy set species 1
  array[N_occ] int occ_s2; // occurrence for occupancy set species 2
  array[N_occ, 2] int ll_idx; // index to occupancy set
  flag = rep_array(0, max(append_array(occ1, occ2)),
                      max(append_array(occ1, occ2)));
  idx = 0;
  for(i in 1:N_dat) {
    occ_min = occ1[i] < occ2[i] ? occ1[i] : occ2[i];
    occ_max = occ1[i] > occ2[i] ? occ1[i] : occ2[i];
    if(flag[occ_min, occ_max] == 0) {
      idx += 1;
      flag[occ_min, occ_max] = idx;
      occ_s1[idx] = occ_min;
      occ_s2[idx] = occ_max;
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
      ll[idx] = lchoose(occ_s1[i], j) + 
                lchoose(N_site - occ_s1[i], occ_s2[i] - j);
      ii[idx] = j;
    }
  }
  // 
  // calculate random index
  array[N_dat, 2] int rng_idx;
  for(i in 1:N_dat) {
    rng_idx[i,1] = (i == 1) ? 1 : rng_idx[i - 1, 2] + 1;
    rng_idx[i,2] = rng_idx[i,1] + nsb[i] - 1;
  }
  //
}
//
//
generated quantities {
  //
  real theta_m;
  real theta_v;
  vector[N_smp] chisq_ob = rep_vector(0, N_smp); // observed sum-of-squares
  vector[N_smp] chisq_sm = rep_vector(0, N_smp); // simulated sum-of-squares
  vector[N_smp] chisq_os; // sum-of-squares ratio
  //
  for(i in 1:N_smp) {
    for(j in 1:N_dat) {
      theta_m = nch_m(ll[ll_idx[occ_id[j], 1]:ll_idx[occ_id[j], 2]],
                      ii[ll_idx[occ_id[j], 1]:ll_idx[occ_id[j], 2]], 
                      theta_vl[i, j]);
      theta_v = nch_v(ll[ll_idx[occ_id[j], 1]:ll_idx[occ_id[j], 2]],
                      ii[ll_idx[occ_id[j], 1]:ll_idx[occ_id[j], 2]], 
                      theta_vl[i, j]);
      chisq_ob[i] += (Y[j] - theta_m)^2/theta_v*nsb[i];
      chisq_sm[i] += pearson_rng(ll[ll_idx[occ_id[j], 1]:ll_idx[occ_id[j], 2]], 
                                 ii[ll_idx[occ_id[j], 1]:ll_idx[occ_id[j], 2]],
                                 theta_vl[i, j], nsb[i]);
    }
    chisq_os[i] = chisq_ob[i]/chisq_sm[i];
  }
  //
  array[20, rng_idx[N_dat, 2]] int y_rep;
  for(i in 1:20) {
    for(j in 1:N_dat) {
      y_rep[i, rng_idx[j, 1]:rng_idx[j, 2]] =
        nch_re_rng(ll[ll_idx[occ_id[j], 1]:ll_idx[occ_id[j], 2]],
                   ii[ll_idx[occ_id[j], 1]:ll_idx[occ_id[j], 2]],
                   theta_mu[i, j], theta_sd[i, j], nsb[j]);
    }
  }
}
