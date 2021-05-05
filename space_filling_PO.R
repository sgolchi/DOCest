library(pscl)
library(MCMCpack)
library(mnormt)
library(abind)


#' Constraint function 
#' @param sample a candidate point, a vector with size 
#' @param l lower limit for each dimension
#' @param u upper limit for each dimension
#' @return deviation vector of the given point from the target space defined by constraint
#' 
mixture = function(sample, l = c(50, 5, 1, 0.5)/100, u = c(90, 30, 5, 2.5)/100) 
  return(c(abs(sum(sample)-1), (l - sample), (sample - u)))

#' adaptive sequence of constraint parameters in SCMC
adapt_seq = function(nu, nut, N, sample, Wt = Wt, constraint = mixture) {
  wt<-c()
  for (i in 1:N) {
    term = constraint(sample[i,])
    cons1 = sum(pnorm(-term/nu, log=T))
    cons2 = sum(pnorm(-term/nut, log=T))
    wt[i] = cons1 - cons2
  }
  Wt = Wt*exp(wt)
  Wt = Wt/sum(Wt)
  ESS = ifelse(sum(is.na(Wt))==N, 0, 1/sum(Wt^2))
  return(ESS - (N/2))
}

#' Initial unrestricted sample for SCMC
unrestricted = function(N, rge) {
  x = cbind(rep(N, nrow(rge)), rge)
  unif = function(x) return(runif(x[1], x[2], x[3]))
  samp = apply(x, 1, unif)
  return(samp)
}

#' constrained log-posterior
logpost = function(sample, t, nuseq, constraint = mixture){
  term = constraint(sample)
  return(sum(pnorm(-term/nuseq[t], log=T)))
}

#' Gibbs step in SCMC
Gibbs = function(x, q, d, a, t=t, logpost.=logpost(), lpdent = lpdent, nuseq) {
  delta = rnorm(1,0,q[d])
  newx = x
  newx[d] = newx[d] + delta
  lpnum = logpost(newx,t, nuseq)
  ratio = lpnum - lpdent
  prob = min(1, exp(ratio))
  u = runif(1)
  if (u <= prob) {
    x = newx
    lpdent = lpnum
    a = a + 1
  }
  return(list(x = x, a = a, lpdent = lpdent))
}

#' SCMC sampler
SCMC = function(D = 4, l = c(50, 5, 1, 0.5)/100, u = c(90, 30, 5, 2.5)/100, 
                epsilon = 1e-3, constraint = mixture) {
  t = 1
  N = 100000
  ESS = c()
  nuseq = c(Inf)
  L = 25
  b = seq(1.5, .1, length = L)
  nuseq_0 = c(Inf, b^7)
  nuseq_T = 1e-4
  rge = cbind(l - epsilon, u + epsilon)
  a = array(0, dim = c(2,D))
  samplet = array(dim = c(N, 1, D))
  lpdent = array(dim = c(N, 1))
  Wt = array(dim = c(N, 1))
  samplet[,1,] = unrestricted(N, rge)
  for (i in 1:N) lpdent[i,1] = logpost(samplet[i,t,], t, nuseq)
  Wt[,1] = rep(1/N, N)
  repeat {
    t = t + 1
    newsample = samplet[,t-1,]
    newlpdent = lpdent[,t-1]
    newWt = Wt[,t-1]
    as = adapt_seq(nuseq_T, nut=nuseq[t-1], N = N, sample = newsample, Wt = newWt)
    if (as>0) {
      nuseq[t] = nuseq_T
    } else {
      #if (t<3) low = min(.1, nuseq_0[t]) else 
      low = nuseq_T
      if (t != 2) up = nuseq[t-1] else up = 1e5
      nextt = uniroot(adapt_seq, interval = c(low,up), nut = nuseq[t-1], N = N, 
                      sample = newsample, Wt = newWt)$root
      nuseq[t] = nextt
    }
    wt = c()
    for (i in 1:N) {
      term = constraint(newsample[i,])
      constraint1 = sum(pnorm(-term/nuseq[t], log=T)) 
      constraint2 = sum(pnorm(-term/nuseq[t-1], log=T))
      wt[i] = constraint1 - constraint2
    }
    newWt = newWt*exp(wt)
    newWt = newWt/sum(newWt)
    ESS[t] = 1/sum(newWt^2)
    index = sample(1:N, N, prob = newWt, replace = T)
    newsample = newsample[index,] 
    newlpdent = newlpdent[index]
    newWt = rep(1/N,N)
    q = apply(newsample, 2, sd)/(t^2)
    a = abind(a, rep(0, D), along = 1)
    for(j in 1:10) {
      for (i in 1:N) {
        for (d in 1:D) {
          out = Gibbs(newsample[i,], q, d, a[t,d], t=t, logpost=logpost, 
                      lpdent=newlpdent[i], nuseq = nuseq)
          newsample[i,] = out$x
          a[t,d] = out$a
          newlpdent[i] = out$lpdent
        }
      }
    }
    samplet = abind(samplet, newsample, along=2)
    lpdent = abind(lpdent, newlpdent, along=2)
    Wt = abind(Wt, newWt, along=2)
    if (nuseq[t] <= nuseq_T) break
  }
  t_final = t
  out = samplet[,t_final,]
  return(out)
}

#out = SCMC()
#save(out,file='S_mixture_experiment.Rdata')



############## matrix plots #############
# S0 = S[,9,]
# panel.hist <- function(x, ...)
# {
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(usr[1:2], 0, 1.5) )
#   h <- hist(x, plot = FALSE)
#   breaks <- h$breaks; nB <- length(breaks)
#   y <- h$counts; y <- y/max(y)
#   rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
# }
# df = data.frame(S0)
# pairs(df,col = 'darkgrey', upper.panel = NULL)

############################################

#load('S_mixture_experiment.Rdata')
S = out

# Clustering design of Lekivetz and Jones (2014)
FFFdesign = function(S, n) {
  hc = hclust(dist(S)^2, method = "ward")
  memb = cutree(hc, k=n)
  cent = NULL
  for(k in 1:n){
    cent = rbind(cent, colMeans(S[memb == k, , drop = FALSE]))
  }
  return(cent)
}


p01 = FFFdesign(S, 20)
p02 = FFFdesign(S, 100)

OR1 = c(.7, .8, 0.9, 1)
OR2 = c(0.65, 0.75, 0.85, 0.95)

# training set
X1 = data.frame(OR = rep(OR1, each = nrow(p01)), p0 = do.call(rbind, rep(list(p01), length(OR1))))

# prediction set
X2 = data.frame(OR = rep(OR2, each = nrow(p02)), p0 = do.call(rbind, rep(list(p02), length(OR2))))


#save(X1, file = 'training.Rdata')
#save(X2, file = 'predict.Rdata')
