len  <- c(530, 840, 930) #lengths of the transcripts
M  <- matrix(c(1, 1, 1, 0, 1, 1, 0, 0, 1), byrow=TRUE, ncol=3) 
k  <- c(1666, 896, 81) #counts for each set
lf  <- 50

NITER=1024
#NITER = 4096

std  <- function(X) {
  stdev  <- sd(X)
  return(stdev/ sqrt(length(X)))
}

GetInitialParams  <- function(k, M, len) {
  # Compute the initial parameters m which
  # correspond to the maximum likelihood estimates
  # given the model
  mt  <- ((k / rowSums(M)) %*% M) / len
  return(mt)
}

GetExpectedObs  <- function(k, M, mt) {
  ms  <- t((M %*% t(mt)))
  e  <- (k / ms) %*% M
  return(e)
}

RunEM <- function(k, M, len) {
  trace_em  <- matrix(nrow=ncol(M), ncol=NITER)
  mt  <- GetInitialParams(k, M, len)
  trace_em[, 1]  <- mt
  for (i in 2:NITER) {
    mt  <- (mt/len) * GetExpectedObs(k, M, mt)
    trace_em[, i]  <- mt
  }
  
  return(trace_em)
}

SampleObs  <-  function(mt, M, k) {
  X  <- matrix(byrow=TRUE, nrow=nrow(M), ncol=ncol(M))
  p  <- matrix(rep(mt, 3), byrow=TRUE, ncol=ncol(M)) * M
  
  for (i in 1:nrow(X)) {
    prob  <- p[i, ]/sum(p[i, ])
    size  <- k[i]
    X[i, ]  <- t(rmultinom(1, size, prob))
  }
  
  return(X)
}

SampleParams  <- function(X, len) {
  a  <- 1.2
  b  <- 0.001
  mt  <- rep(0, ncol(X))
  r  <- colSums(X)
  
  for (i in 1:ncol(X)) {
    mt[i]  <- rgamma(1, shape=a+r[i], rate=b+len[i]) 
  }
  
  return(mt)
}

RunGibbs  <- function(mt, M, k, len) {
  trace_gibbs  <- matrix(nrow=ncol(M), ncol=NITER)
  X  <- SampleObs(mt, M, k)
  m  <- SampleParams(X, len)
  
  trace_gibbs[, 1]  <- m
  for (i in 2:NITER) {
    X  <- SampleObs(m, M, k)
    m  <- SampleParams(X, len)
    trace_gibbs[, i]  <- m
  }
  
  return(trace_gibbs)
}

PlotTraces <- function(trace_em, trace_gibbs) {
  nt  <- nrow(trace_em)
  #nt  <- 1
  tr  <- as.factor(unlist(Map(rep, paste("t", 1:nt, sep=""), NITER)))
  #tr  <- rep("t1", NITER)
  unname(tr)
  em.traces  <- as.vector(t(trace_em))
  gibbs.traces  <- as.vector(t(trace_gibbs))
  traces  <- c(gibbs.traces, em.traces)
  alg  <- c(rep("gibbs", nt*NITER), rep("em", nt*NITER))
  g  <- data.frame(vals=traces, tr=tr, alg=alg, n=1:NITER, each=NITER)
  ggplot(data=g, aes(x=n, y=vals)) + geom_line(aes(colour=tr, linetype=alg)) + 
    scale_linetype_manual(values = c("dashed", "solid"))
}

GetNormLength  <- function(lt) {
  mf  <- 180
  sdf  <- 30
  lr  <- 30
  
  poss.fl  <- seq(lr, lt)
  nl  <- sum(dnorm(poss.fl, mean=mf, sd=sdf) * (lt-poss.fl+1))
  
  return(floor(nl))
}

# firstly run EM algorithm to ge initial estimates for mu
# use the EM values as initial values for the Gibbs sampling
# then plot the traces

#first normalise the lengths
#len  <- sapply(len, GetNormLength)
print(len)
trace_em  <- RunEM(k, M, len)
mt  <- trace_em[, NITER]
trace_gibbs  <- RunGibbs(mt, M, k, len)

# TODO some more experimenting with 4 transcripts with t3 and t4 being identical







