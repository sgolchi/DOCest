load('training.Rdata')
load('predict.Rdata')

cases = X1
df = NULL
for (l in 1:nrow(cases)) {
  ps = psup_dist(cases$OR[l], pi_c = c(cases$p0.1[l], cases$p0.2[l], cases$p0.3[l], 1 - (cases$p0.1[l] + cases$p0.2[l] + cases$p0.3[l])), 1000)
  mu = mean(ps)
  sig2 = var(ps)
  alpha = uniroot(asolve, c(.01, 1000), mu = mu, sig2 = sig2)$root
  beta = alpha*(1-mu)/mu
  df = rbind(df, data.frame(alpha = alpha, beta = beta, pi_c = cbind(cases$p0.1[l], cases$p0.2[l], cases$p0.3[l], 1 - (cases$p0.1[l] + cases$p0.2[l] + cases$p0.3[l])), OR = cases$OR[l],
                            tpow = mean(ps>.95)))
}

#save(df, file = 'sims_po_Dec16.Rdata')

#####################################################################
#load('sims_po_Dec16.Rdata')

#' Prediction function
#' @param rho GP correlation parameters vector 
#' @param alpha GP variance parameter
#' @param sigma observation standard error
#' @param X1 training set inputs
#' @param X2 prediction set
#' @param y1 training set output
#' @return vector of predicted outputs at X2
#' 
pred = function(rho, alpha, sigma, X1, X2, y1) {
  N1 = nrow(X1)
  N2 = nrow(X2)
  D = ncol(X1)
  K = matrix(NA, N1, N1)
  for (i in 1:(N1 - 1)) {
    K[i, i] = alpha + sigma^2;
    for (j in (i + 1):N1) {
      v = c()
      for (d in 1:D) v[d] = -rho[d] * (X1[i,d] - X1[j,d])^2
      K[i,j] = alpha*exp(sum(v));
      K[j, i] = K[i, j];
    }
  }
  K[N1, N1] = alpha + sigma^2
  L_K = Rfast::cholesky(K)
  K2 = matrix(NA, N2, N2)
  for (i in 1:(N2 - 1)) {
    K2[i, i] = alpha# + sigma^2
    for (j in (i + 1):N2) {
      v = c()
      for (d in 1:D) v[d] = -rho[d] * (X2[i,d] - X2[j,d])^2
      K2[i,j] = alpha*exp(sum(v))
      K2[j, i] = K2[i, j]
    }
  }
  K2[N2, N2] = alpha# + sigma^2
  L_K2 = Rfast::cholesky(K2)
  k_x1_x2 = matrix(NA, N1, N2)
  for (i in 1:N1) {
    for (j in 1:N2) {
      v = c()
      for (d in 1:D) v[d] = -rho[d] * (X1[i,d] - X2[j,d])^2
      k_x1_x2[i,j] = alpha*exp(sum(v))
    }
  }
  #Kinv = chol2inv(K)
  Kinv = solve(K)
  K22 = K2 - t(k_x1_x2)%*%Kinv%*%k_x1_x2
  preds = t(k_x1_x2)%*%Kinv%*%y1
  out = list(mean = preds, cov = K22)
  return(out)
}

# GP fit for alpha's
y1 = df$alpha#((log(df$alpha) - mean(log(df$alpha)))/sd(log(df$alpha)))
data = list(N1 = nrow(X1), D = ncol(X1), x1 = X1, y1 = y1, N2 = nrow(X2), x2 = X2)
fit = stan(file='GP.stan', data=data, chains=0)
fit = stan('fit'=fit, 'data'=data, warmup=1000, iter=2000, chains=2)
fitss = rstan::extract(fit)

# GP fit for beta's
df$beta[df$beta==0] = 1e-4
y11 = df$beta#((log(df$beta) - mean(log(df$beta)))/sd(log(df$beta)))
data2 = list(N1 = nrow(X1), D = ncol(X1), x1 = X1, y1 = y11, N2 = nrow(X2), x2 = X2)
#fit = stan(file='GP.stan', data=data, chains=0)
fit2 = stan('fit'=fit, 'data'=data2, warmup=1000, iter=2000, chains=2)
fitss2 = rstan::extract(fit2)

# A sample from the posterior
samp = sample(1:2000, 100)
apreds = NULL
bpreds = NULL
for (i in 1:100) {
  ap = pred(fitss$rho[samp[i],], fitss$alpha[samp[i]], fitss$sigma[samp[i]], X1, X2, y1)
  apreds = cbind(apreds, t(MASS::mvrnorm(n = 100, ap$mean, ap$cov)))
  bp = pred(fitss2$rho[samp[i],], fitss2$alpha[samp[i]], fitss2$sigma[samp[i]], X1, X2, y11)
  bpreds = cbind(bpreds, t(MASS::mvrnorm(n = 100, bp$mean, bp$cov)))
}

alphahats = apreds
# discarding non-positive parameter values
for (i in 1:400) alphahats[i,alphahats[i,]<0] = sample(alphahats[i,alphahats[i,]>0], length(alphahats[i,alphahats[i,]<0]))
betahats = bpreds
# discarding non-positive parameter values
for (i in 1:400) betahats[i,betahats[i,]<0] = sample(betahats[i,betahats[i,]>0], length(betahats[i,betahats[i,]<0]))

# power and probability of concluding futility estimates 
pow = matrix(NA, 400, 10000)
fut = matrix(NA, 400, 10000)
for (i in 1:400) {
  for (j in 1:10000) {
    pow[i, j] = 1 - pbeta(.95, alphahats[i,j], betahats[i,j])
    fut[i,j] = pbeta(.05, alphahats[i,j], betahats[i,j])
   }
}


df1 = data.frame(pow_hat = apply(pow, 1, mean), pow_low = apply(pow, 1, quantile, 0.025), 
                 pow_up = apply(pow, 1, quantile, 0.975),
                 fut_hat = apply(fut, 1, mean), fut_low = apply(fut, 1, quantile, 0.025), 
                 fut_up = apply(fut, 1, quantile, 0.975),
                 alpha_hat = apply(alphahats, 1, mean), 
                 alpha_low = apply(alphahats, 1, quantile, 0.025), 
                 alpha_up = apply(alphahats, 1, quantile, 0.975),
                 beta_hat = apply(betahats, 1, mean),
                 beta_low = apply(betahats, 1, quantile, 0.025), 
                 beta_up = apply(betahats, 1, quantile, 0.975),
                 OR = X2[,1], p0 = X2[,2:4])

## plotting
 # ggplot(df1, aes(x = p0.p0.1, y = pow_hat)) + 
 #   geom_errorbar(data = df1, aes(x = p0.p0.1, ymin=pow_low, ymax=pow_up), color = 'grey60', width=0.008) +
 #   geom_point(color = 'grey40', alpha = 1, size = 1.25) +
 #   geom_point(data = df, aes(x = pi_c.1, y = tpow), color ='grey40', shape = 15) +
 #   facet_grid(.~ OR, labeller = label_both) + ylab('probability of concluding superiority') + xlab(expression(p[1])) +
 #   geom_hline(yintercept = 0.8, col = 'grey40', alpha = 0.7) +
 #   theme_bw() +
 #   theme(axis.title = element_text(size = 14, face = 'bold'),
 #         axis.text = element_text(size = 14),
 #         strip.text = element_text(size = 14))
 #   
 # ggsave('power_pred_obs03.png')
 # 
 # ggplot(df1, aes(x = p0.p0.1, y = pow_hat)) + 
 #   geom_errorbar(data = df1, aes(x = p0.p0.1, ymin=pow_low, ymax=pow_up), color = 'grey60', width=0.008) +
 #   geom_point(color = 'grey40', alpha = 1, size = 1.25) +
 #   #geom_line(data = df1, aes(x = p0.p0.1, y = pow_low, color = as.factor(OR))) +
 #   # geom_line(data = df1, aes(x = p0.p0.1, y = pow_up, color = as.factor(OR))) +
 #   #geom_point(data = df, aes(x = pi_c.1, y = tpow), color ='grey40', shape = 15) +
 #   facet_grid(.~ OR, labeller = label_both) + ylab('probability of concluding superiority') + xlab(expression(p[1])) +
 #   geom_hline(yintercept = 0.8, col = 'darkred', alpha = 0.7) +
 #   theme_bw() +
 #   theme(axis.title = element_text(size = 14, face = 'bold'),
 #         axis.text = element_text(size = 14))
 # 
 # ggsave('power_90.png')
 # 
 # df01 = df1 %>% 
 #   filter(OR %in% c(0.85, 0.95))
 # ggplot(df01, aes(x = p0.p0.1, y = fut_hat)) + 
 #   geom_errorbar(data = df01, aes(x = p0.p0.1, ymin=fut_low, ymax=fut_up), color = 'grey60', width=0.008) +
 #   geom_point(color = 'grey40', alpha = 1, size = 1.25) +
 #   #geom_line(data = df1, aes(x = p0.p0.1, y = pow_low, color = as.factor(OR))) +
 #   # geom_line(data = df1, aes(x = p0.p0.1, y = pow_up, color = as.factor(OR))) +
 #   #geom_point(data = df, aes(x = pi_c.1, y = tpow), color ='grey40', shape = 15) +
 #   facet_grid(.~ OR, labeller = label_both) + ylab('probability of concluding futility') + xlab(expression(p[1])) +
 #   #geom_hline(yintercept = 0.8, col = 'darkred', alpha = 0.7) +
 #   theme_bw() + ylim(c(0, 0.219)) +
 #   theme(axis.title = element_text(size = 14, face = 'bold'),
 #         axis.text = element_text(size = 14))
 # ggsave('fut_01.png')
 # 
 # ggplot(df1, aes(x = p0.p0.1, y = alpha_hat)) + 
 #   geom_errorbar(data = df1, aes(x = p0.p0.1, ymin=alpha_low, ymax=alpha_up), color = 'grey60') +
 #   geom_point(color = 'grey40') +
 #   # geom_line(data = df1, aes(x = p0.p0.1, y = pow_low, color = as.factor(OR))) +
 #   # geom_line(data = df1, aes(x = p0.p0.1, y = pow_up, color = as.factor(OR))) +
 #   geom_point(data = df, aes(x = pi_c.1, y = alpha), color ='grey40', shape = 15) +
 #   facet_grid(.~ OR, labeller = label_both) + ylab(expression(alpha)) + 
 #   xlab(expression(p[1])) +
 #   theme_bw() +
 #   theme(axis.title = element_text(size = 14, face = 'bold'),
 #         axis.text = element_text(size = 14),
 #         strip.text = element_text(size = 14))
 # ggsave('alpha.png')
 # 
 # ggplot(df1, aes(x = p0.p0.1, y = beta_hat)) + 
 #   geom_errorbar(data = df1, aes(x = p0.p0.1, ymin=beta_low, ymax=beta_up), color = 'grey60') +
 #   geom_point(color = 'grey40') +
 #   # geom_line(data = df1, aes(x = p0.p0.1, y = pow_low, color = as.factor(OR))) +
 #   # geom_line(data = df1, aes(x = p0.p0.1, y = pow_up, color = as.factor(OR))) +
 #   geom_point(data = df, aes(x = pi_c.1, y = beta), color ='grey40', shape = 15) +
 #   facet_grid(.~ OR, labeller = label_both) + ylab(expression(beta)) + 
 #   xlab(expression(p[1])) +
 #   theme_bw() +
 #   theme(axis.title = element_text(size = 14, face = 'bold'),
 #         axis.text = element_text(size = 14),
 #         strip.text = element_text(size = 14))
 # ggsave('beta.png')


