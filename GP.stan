data {
  int<lower=1> N1;
  int<lower=1> D;
  vector[D] x1[N1];
  vector[N1] y1;
  int<lower=1> N2;
  vector[D] x2[N2];
}
transformed data {
  vector[N1] mu = rep_vector(0, N1);
  real delta = 1e-6;
}
parameters {
  vector<lower=0>[D] rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
model {
  real sq_sigma = square(sigma);
  matrix[N1, N1] L_K;
{
  matrix[N1, N1] K;
  for (i in 1:(N1 - 1)) {
    K[i, i] = alpha + sq_sigma;
    for (j in (i + 1):N1) {
      vector[D] v;
      for (d in 1:D) {
        v[d] = -rho[d] * square(x1[i,d] - x1[j,d]);
  }
  K[i,j] = alpha*exp(sum(v));
  K[j, i] = K[i, j];
}
}
K[N1, N1] = alpha + sq_sigma;
L_K = cholesky_decompose(K);
}
rho ~ inv_gamma(5, 5);
alpha ~ normal(0, 1);
sigma ~ normal(0, 1);
y1 ~ multi_normal_cholesky(mu, L_K);
}
