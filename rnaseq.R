std  <- function(X) {
  # Computes the standard error for either a matrix X in
  # which case it computer it along the rows or a vector X
  if (is.null(dim(X))) {
    stdev  <- sd(X)
    se  <- stdev/ sqrt(length(X))
  } else {
    stdev  <- apply(X, 1, sd)
    se  <- stdev / sqrt(length(X[1, ]))
  }
  return(se)
}

GetInitialParams  <- function(k, M, len) {
  # Computes the initial parameters m which
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

PlotTraces <- function(trace_em, trace_gibbs, tr) {
  nt  <- nrow(trace_gibbs)
  em.traces  <- as.vector(t(trace_em))
  gibbs.traces  <- as.vector(t(trace_gibbs))
  traces  <- c(gibbs.traces, em.traces)
  alg  <- c(rep("Gibbs", nt*NITER), rep("EM", nt*NITER))
  g  <- data.frame(mu=traces, tr=tr, alg=alg, n=1:NITER, each=NITER)
  ggplot(data=g, aes(x=n, y=mu)) + geom_line(aes(colour=tr, linetype=alg)) + 
    scale_linetype_manual(values = c("dashed", "solid")) + xlab("iteration")
}

GetEffLength  <- function(lt) {
  # Returns the expected number of positions that 
  # a fragment can start from along the transcript
  mf  <- 180
  sdf  <- 25
  lr  <- 30
  
  poss.fl  <- seq(lr, lt)
  nl  <- sum(dnorm(poss.fl, mean=mf, sd=sdf) * (lt-poss.fl+1))
  
  return(floor(nl))
}

RunExp1  <- function() {
  # Runs the first experiment with the 3 initial transcripts
  # and plots the traces obtained from EM and Gibbs algorithms
  NITER=1024
  len  <- c(530, 840, 930) #lengths of the transcripts
  M  <- matrix(c(1, 1, 1, 0, 1, 1, 0, 0, 1), byrow=TRUE, ncol=3) 
  k  <- c(1666, 896, 81) #counts for each set
  trace_em  <- RunEM(k, M, len)
  mt  <- trace_em[, NITER]
  trace_gibbs  <- RunGibbs(mt, M, k, len)
  
  nt  <- nrow(trace_em)
  tr  <- as.factor(unlist(Map(rep, paste("transcript", 1:nt, sep=""), NITER)))
  PlotTraces(trace_em, trace_gibbs, tr)
  
}

RunExp2  <- function() {
  # Runs the second experiment with the 2 identical transcripts
  # and plots the traces obtained from EM and Gibbs algorithms
  NITER=1024
  len  <- c(530, 840, 930, 930) #lengths of the transcripts
  M  <- matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1), byrow=TRUE, ncol=4) 
  k  <- c(1666, 896, 81) #counts for each set
  trace_em  <- RunEM(k, M, len)
  mt  <- trace_em[, NITER]
  trace_gibbs  <- RunGibbs(mt, M, k, len)
  ident_em  <- rbind(trace_em[3:4, ], trace_em[3, ]+trace_em[4,])
  ident_gibbs  <- rbind(trace_gibbs[3:4, ], trace_gibbs[3, ]+trace_gibbs[4, ])
  tr  <- as.factor(c(rep("transcript 3", NITER), rep("transcript 4", NITER), rep("transcript 3 + transcript 4", NITER)))
  PlotTraces(ident_em, ident_gibbs, tr) 
}

