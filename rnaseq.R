len  <- c(530, 840, 930) #lengths of the transcripts
M  <- matrix(c(1, 1, 1, 0, 1, 1, 0, 0, 1), byrow=TRUE, ncol=3) 
k  <- c(1666, 896, 81) #counts for each set

NITER=1024

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
  t  <- as.vector(Map(rep, paste("t", 1:nt, sep=""), 1024))
  em.traces  <- as.vector(t(trace_em))
  gibbs.traces  <- as.vector(t(trace_gibbs))
  traces  <- c(gibbs.traces, em.traces)
  alg  <- c(rep("gibbs", nt*NITER), rep("em", nt*NITER))
  g  <- data.frame(vals=traces, t=t, alg=alg, n=1:1024, each=1024)
  ggplot(data=g, aes(x=n, y=vals)) + geom_line(aes(colour=t, linetype=alg)) + scale_linetype_manual(values = c("dashed", "solid"))
  
}







