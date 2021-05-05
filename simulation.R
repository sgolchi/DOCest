load('simulation_training_sets.Rdata')

## Construct test set and obtain trial simulation based measures (true responses) for 
# the test set
OR = round(seq(.6, 1.05, .05), 2)
p0 = seq(0.25, 0.7, .05)
test_set = expand.grid(OR, p0)
names(test_set) = (c('OR', 'p0'))

#' Relationship between beta shape parameter, the mean and variance
#' @param a beta shape parameter
#' @param mu mean of the beta distribution
#' @param sig2 variance of the beta distribution
#' @return zero for correct parameters

asolve = function(a, mu, sig2) (a^2*(1-mu)/mu)/((a + a*(1-mu)/mu)^2 * (a + a*(1-mu)/mu + 1)) - sig2


#' probability of superiority 
#'
#' @param p0 true base risk
#' @param OR true odds ratio
#' @param t OR threshold, default is 1
#' @param N trial sample size
#' @return probability that OR<t
#'
psup = function(p0, OR, t=1, N) {
  p1 = (OR*(p0/(1-p0)))/(1+OR*(p0/(1-p0)))
  ps = NULL
  for (m in 1:10000) {
    c = rbinom(N, 1, .5)
    y0 = rbinom(sum(c), 1, p0)
    y1 = rbinom(N-sum(c), 1, p1)
    a0 = sum(y0) + 1
    b0 = N - sum(y0) + 1
    a1 = sum(y1) + 1
    b1 = N - sum(y1) + 1
    samp0 = rbeta(1000, a0, b0)
    samp1 = rbeta(1000, a1, b1)
    OR_post = (samp1/(1-samp1))/(samp0/(1-samp0))
    ps = c(ps, mean(OR_post<t))
  }
  return(ps)
}

# estimate power by trial simulation
cases = expand.grid(OR, p0)
names(cases) = c('OR', 'p0')
df = NULL
for (l in 1:nrow(cases)) {
  ps = psup(p0 = cases$p0[l], OR = cases$OR[l], t = 1, N = 1000)
  mu = mean(ps)
  sig2 = var(ps)
  alpha = uniroot(asolve, c(.01, 1000), mu = mu, sig2 = sig2)$root
  beta = alpha*(1-mu)/mu
  df = rbind(df, data.frame(alpha = alpha, beta = beta, p0 = cases$p0[l], OR = cases$OR[l],
                            tpow = mean(ps>.95)))
}
save(df, file = 'truth.Rdata')
truth = df


fit = stan(file='GP0.stan', data=data, chains=0)
##############################################################
### simulation function for assessing surrogates trained over the 100 training sets generated
sim = function(i, train_set, test_set) {
  x1 = cbind(train_set$p0[train_set$iter==i], train_set$OR[train_set$iter==i])
  x2 = cbind(test_set$p0, test_set$OR)
  y1 = train_set$alpha[train_set$iter==i]#(log(alpha1) - mean(log(alpha1)))/sd(log(alpha1))
  data = list(N1 = nrow(x1), D = 2, x1 = x1, y1 = y1, N2 = nrow(x2), x2 = x2)
  fit = stan('fit'=fit, 'data'=data, warmup=1000, iter=2000, chains=1)
  fitss = rstan::extract(fit)
  y11 = train_set$beta[train_set$iter==i]#(log(beta1) - mean(log(beta1)))/sd(log(beta1))
  data2 = list(N1 = nrow(x1), D = 2, x1 = x1, y1 = y11, N2 = nrow(x2), x2 = x2)
  fit2 = stan('fit'=fit, 'data'=data2, warmup=1000, iter=2000, chains=1)
  fitss2 = rstan::extract(fit2)
  samp = sample(1:1000, 100)
  apreds = NULL
  bpreds = NULL
  for (i in 1:100) {
    ap = pred(fitss$rho[samp[i],], fitss$alpha[samp[i]], fitss$sigma[samp[i]], x1, x2, y1)
    apreds = cbind(apreds, t(tryCatch(MASS::mvrnorm(n = 100, ap$mean, ap$cov),error = function(e){ next })))
    bp = pred(fitss2$rho[samp[i],], fitss2$alpha[samp[i]], fitss2$sigma[samp[i]], x1, x2, y11)
    bpreds = cbind(bpreds, t(tryCatch(MASS::mvrnorm(n = 100, bp$mean, bp$cov),error = function(e){ next })))
  }
  alphahats = apreds#exp(apreds*sd(log(df$alpha)) + mean(log(df$alpha)))
  for (i in 1:100) alphahats[i,alphahats[i,]<0] = sample(alphahats[i,alphahats[i,]>0], length(alphahats[i,alphahats[i,]<0]))
  betahats = bpreds#exp(bpreds*sd(log(df$beta)) + mean(log(df$beta)))
  for (i in 1:100) betahats[i,betahats[i,]<0] = sample(betahats[i,betahats[i,]>0], length(betahats[i,betahats[i,]<0]))
  pow = matrix(NA, 100, 10000)
  for (i in 1:100) for (j in 1:10000) pow[i, j] = 1 - pbeta(.95, alphahats[i,j], betahats[i,j])
  RMSE = c()
  bias = c()
  pse = c()
  cov = c()
  for(k in 1:nrow(x2)) {
    RMSE[k] = sqrt(mean((pow[k,]-test_set$tpow[k])^2))
    bias[k] = mean((pow[k,])-test_set$tpow[k])
    pse[k] = sd(pow[k,])
    cov[k] = test_set$tpow[k]>quantile(pow[k,], 0.025) & test_set$tpow[k]<quantile(pow[k,], 0.975)
  }
  out = data.frame(RMSE = RMSE, bias = bias, pse = pse, cov = cov)
  return(out)
}

load('truth.Rdata')

#sim01 = NULL
for (l in 1:100){
  s1 = sim(l, train_set = train_set, test_set = test_set)
  sim01 = rbind(sim01, cbind(iter = rep(l,100), s1))
}

#save(sim01, file = 'sim01.Rdata')

##########################################################
### plots
# load('sim01.Rdata')
# dfp = cbind(OR = rep(test_set$OR, 100), p0 = rep(test_set$p0, 100), sim01)
# dfp = dfp %>% filter(OR != 0.6 & OR!=1.05)
# ggplot(dfp, aes(x = as.factor(p0), y = RMSE)) + geom_boxplot(color = 'grey40') + 
#   facet_grid(.~OR, labeller = label_both) + theme_bw() + xlab(expression(p[0])) +
#   theme(axis.title = element_text(size = 14, face = 'bold'),
#         axis.text.y = element_text(size = 14),
#         axis.text.x = element_text(angle = 90),
#         legend.title = element_text(size = 14, face = 'bold'))
# ggsave(file = 'RMSE_box.png', width = 12)
# ggplot(dfp, aes(x = as.factor(p0), y = bias)) + geom_boxplot(color = 'grey40') + 
#   facet_grid(.~OR, labeller = label_both) + theme_bw() + xlab(expression(p[0])) +
#   geom_hline(yintercept = 0, color = 'darkred', alpha = 0.5) +
#   theme(axis.title = element_text(size = 14, face = 'bold'),
#         axis.text.y = element_text(size = 14),
#         axis.text.x = element_text(angle = 90),
#         legend.title = element_text(size = 14, face = 'bold'))
# ggsave(file = 'bias_box.png', width = 12)
# ggplot(dfp, aes(x = as.factor(p0), y = pse)) + geom_boxplot(color = 'grey40') + 
#   facet_grid(.~OR, labeller = label_both) + theme_bw() + xlab(expression(p[0])) +
#   ylab('PSE') +
#   theme(axis.title = element_text(size = 14, face = 'bold'),
#         axis.text.y = element_text(size = 14),
#         axis.text.x = element_text(angle = 90),
#         legend.title = element_text(size = 14, face = 'bold'))
# ggsave(file = 'pse_box.png', width = 12)
# 
# dfs = dfp %>%
#   group_by(p0, OR) %>%
#   summarize(RMSE = mean((RMSE)), BIAS = mean(bias), pse = mean(pse), coverage = mean(cov))
# dim(dfs)
# 
# ggplot(dfs, aes(x = p0, y = RMSE, color = as.factor(OR))) + geom_point() + 
#   geom_line(aes(linetype = as.factor(OR))) 
# 
# ggplot(dfs, aes(x = p0, y = BIAS, color = as.factor(OR))) + geom_point() + 
#   geom_line(aes(linetype = as.factor(OR))) 
# 
# ggplot(dfs, aes(x = p0, y = OR, z = BIAS)) + geom_contour_filled(alpha= 0.85) + 
#   geom_contour(color = "grey", size = 0.1) +  xlab(expression(p[0])) +
#   labs(fill = 'bias') +
#   theme_minimal() +
#   theme(axis.title = element_text(size = 14, face = 'bold'),
#         axis.text = element_text(size = 14),
#         legend.title = element_text(size = 14, face = 'bold'),
#         legend.text = element_text(size = 14))
# ggsave(file = 'bias.png')
# ggplot(dfs, aes(x = p0, y = OR, z = RMSE)) + geom_contour_filled(alpha= 0.85) + 
#   geom_contour(color = "grey", size = 0.1) + xlab(expression(p[0])) +
#   labs(fill = 'RMSE') +
#   theme_minimal() +
#   theme(axis.title = element_text(size = 14, face = 'bold'),
#         axis.text = element_text(size = 14),
#         legend.title = element_text(size = 14, face = 'bold'),
#         legend.text = element_text(size = 14))
# ggsave(file = 'RMSE.png')
# ggplot(dfs, aes(x = p0, y = OR, z = pse)) + geom_contour_filled(alpha= 0.85) + 
#   geom_contour(color = "grey", size = 0.1) + xlab(expression(p[0])) +
#   #scale_fill_brewer(palette = 'Set1') + 
#   labs(fill = 'PSE') +
#   theme_minimal()+
#   theme(axis.title = element_text(size = 14, face = 'bold'),
#         axis.text = element_text(size = 14),
#         legend.title = element_text(size = 14, face = 'bold'),
#         legend.text = element_text(size = 14))
# ggsave(file = 'pse.png')
# 
# 
# summary = sim01 %>%
#   group_by(iter) %>%
#   summarize(cover = mean(cov), RMSE = mean(RMSE), bias = mean(bias), pse = mean(pse))
# 
