data {
  int nobs;
  vector[nobs] z;
  vector[nobs] mu;
  vector[nobs] dmu;
}

parameters {
  real<lower=0> H0;
  real<lower=0,upper=1> omega_m;
  real<lower=-1, upper=1> omega;
}

transformed parameters {
  vector[nobs] d_L = 1.0/H0  *(z+z .*z * ((3-10*omega+3*square(omega)+10*omega*omega_m+6*square(omega)*omega_m-9*square(omega)*square(omega_m))/(4*(1-3*omega+3*omega*omega_m)))) ./(1+z *((1-2*omega -3*square(omega)+2*omega*omega_m+12*square(omega)*omega_m-9*square(omega)*square(omega_m))/(2*(1-3*omega+3*omega*omega_m))));
  vector[nobs] mu_t = 5*log10(d_L/10);
}

model {
  /* Flat in m and b. */

  mu ~ normal(mu_t, dmu);
}
