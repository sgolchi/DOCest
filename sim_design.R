

## generating 100 designs for the training step
sim_des = NULL
for (i in 1:100) {
  OR = runif(100, 0.6, 1)
  p0 = runif(100, 0.25, 0.7)
  cases = data.frame(cbind(OR, p0))
  names(cases) = c('OR', 'p0')
  hc = hclust(dist(cases)^2, method = "ward")
  memb = cutree(hc, k=n)
  cent = NULL
  for(k in 1:n){
    cent = rbind(cent, colMeans(cases[memb == k, , drop = FALSE]))
  }
  sim_des = rbind(sim_des, cbind(iter = rep(i, n), cent))
}

sim_des = as.data.frame(sim_des)

df = NULL
for (l in 1:nrow(sim_des)) {
  ps = psup(p0 = sim_des$p0[l], OR = sim_des$OR[l], t = 1, N = 1000)
  mu = mean(ps)
  sig2 = var(ps)
  alpha = uniroot(asolve, c(.01, 1000), mu = mu, sig2 = sig2)$root
  beta = alpha*(1-mu)/mu
  df = rbind(df, data.frame(alpha = alpha, beta = beta))
}
sim_des$alpha = df$alpha
sim_des$beta = df$beta

#save(sim_des, file = 'simulation_training_sets.Rdata')


