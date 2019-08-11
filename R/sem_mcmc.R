mcmcv5 <- function(ct, u1, u2, X, ...)
{
  
  n  <- nrow(X)  # sample size
  p  <- ncol(X)  # dimension
  q1 <- ncol(u1)
  q2 <- ncol(u2)
  
  
  # Survival response
  time         <- ct[, 1]
  status       <- ct[, 2]
  censored.id  <- which(status == 0)
  n.censored   <- length(censored.id)  # number of censored observations
  X.censored   <- X[censored.id, ]
  X.observed   <- X[-censored.id, ]
  logt <- logtime <- log(time)   # for coding convenience, since the whole code is written with y
  logt.censored <- logt[censored.id]
  logt.observed <- logt[-censored.id]
  
  
  alpha.t  <- runif(n = 1, min = -1, max = 1)  # intercept term
  alpha.u1 <- runif(n = q1, min = -1, max = 1)  # intercept term
  alpha.u2 <- runif(n = q2, min = -1, max = 1)  # intercept term
  beta.t   <- runif(n = p, min = -1, max = 1)  # regression parameter 
  phi.t    <- 1
  phi.u1   <- rep(1, q1)
  phi.u2   <- rep(1, q2)
  
  sigma.t.square    <- 1
  sigma.u1.square   <- 1
  sigma.u2.square   <- 1
  sigma.eta1.square <- 1
  sigma.eta2.square <- 1
  
  eta2 <- rnorm(n = 1, mean = 0, sd = 1)
  eta1 <- rnorm(n = 1, mean = eta2, sd = 1)
  
  nburnin <- nburnin
  nmc     <- nmc
  niter   <- nburnin + nmc
  effsamp <- (niter - nburnin)/nthin
  
  # output
  beta.tout   <- matrix(0, nrow = p, ncol = effsamp)
  alpha.tout  <- rep(0, effsamp)
  phi.tout    <- rep(0, effsamp)
  alpha.u1out <- matrix(0, nrow = q1, ncol = effsamp)
  alpha.u2out <- matrix(0, nrow = q2, ncol = effsamp)
  phi.u1out   <- matrix(0, nrow = q1, ncol = effsamp)
  phi.u2out   <- matrix(0, nrow = q2, ncol = effsamp)
  eta1.out    <- rep(0, effsamp)
  eta2.out    <- rep(0, effsamp)
  
  sigma.t.square.out  <- rep(0, effsamp)
  sigma.u1.square.out <- rep(1, effsamp)
  sigma.u2.square.out <- rep(1, effsamp)
  
  logt.out          <- matrix(0, nrow = n, ncol = effsamp)
  logt.hat.out      <- matrix(0, nrow = n, ncol = effsamp)
  loglikelihood.out <- rep(0, effsamp)
  likelihood.out    <- matrix(0, nrow = n, ncol = effsamp)
  
  for(iter in 1:niter)  # MCMC
  {
    
    ## Update survival latent variable ##
    mean.impute       <- alpha.t + as.matrix(X.censored) %*% as.matrix(beta.t) + eta1 * phi.t
    sd.impute         <- sqrt(sigma.t.square)
    time.censored     <- rtnorm(n.censored, mean = mean.impute, sd = sd.impute, lower = logt.censored)
    logt[censored.id] <- time.censored  # truncated at log(time) for censored data
    
    
    # Sample $ \beta.t $
    A            <- crossprod(x = X) + chol2inv(chol(100 * diag(p)))
    Ainv         <- chol2inv(chol(A))
    Sigma.beta.t <- sigma.t.square * Ainv
    mean.beta.t  <- as.vector(Ainv %*% t(X) %*% (logt - alpha.t - eta1 * phi.t))
    beta.t       <- as.vector(mvrnorm(n = 1, mu = mean.beta.t, Sigma = Sigma.beta.t))
    
    # Sample $ \alpha_t $
    A             <- crossprod(x = rep(1, n)) + chol2inv(chol(diag(1)))
    Ainv          <- chol2inv(chol(A))
    Sigma.alpha.t <- sigma.t.square * Ainv
    mean.alpha.t  <- as.vector(Ainv %*% t(rep(1, n)) %*% (logt - X %*% beta.t - eta1 * phi.t))
    alpha.t       <- as.vector(mvrnorm(n = 1, mu = mean.alpha.t, Sigma = Sigma.alpha.t))
    
    # Sample $ \phi_t $
    A           <- crossprod(x = rep(eta1, n)) + chol2inv(chol(diag(1)))
    Ainv        <- chol2inv(chol(A))
    Sigma.phi.t <- Ainv
    mean.phi.t  <- as.vector(Ainv %*% t(rep(eta1, n)) %*% (logt - alpha.t - X %*% beta.t))
    phi.t       <- as.vector(mvrnorm(n = 1, mu = mean.phi.t, Sigma = Sigma.phi.t))
    
    # Sample $ \alpha_u1 $
    A             <- crossprod(x = rep(1, n)) + chol2inv(chol(diag(1)))
    Ainv          <- chol2inv(chol(A))
    Sigma.alpha.u <- sigma.u1.square * Ainv
    mean.alpha.u  <- as.vector(Ainv %*% t(rep(1, n)) %*% (u1 - eta1 * phi.u1))
    alpha.u1      <- as.vector(rnorm(n = q1, mean = mean.alpha.u, sd = sqrt(Sigma.alpha.u)))
    
    # Sample $ \alpha_u2 $
    A             <- crossprod(x = rep(1, n)) + chol2inv(chol(diag(1)))
    Ainv          <- chol2inv(chol(A))
    Sigma.alpha.u <- sigma.u2.square * Ainv
    mean.alpha.u  <- as.vector(Ainv %*% t(rep(1, n)) %*% (u2 - eta2 * phi.u2))
    alpha.u2      <- as.vector(rnorm(n = q2, mean = mean.alpha.u, sd = sqrt(Sigma.alpha.u)))
    
    # Sample $ \phi_u1 $
    A           <- crossprod(x = rep(eta1, n)) + chol2inv(chol(diag(1)))
    Ainv        <- chol2inv(chol(A))
    Sigma.phi.u <- Ainv
    mean.phi.u  <- as.vector(Ainv %*% t(rep(eta1, n)) %*% (u1 - alpha.u1))
    phi.u1      <- as.vector(rnorm(n = q1, mean = mean.phi.u, sd = sqrt(Sigma.phi.u)))
    
    # Sample $ \phi_u2 $
    A           <- crossprod(x = rep(eta2, n)) + chol2inv(chol(diag(1)))
    Ainv        <- chol2inv(chol(A))
    Sigma.phi.u <- Ainv
    mean.phi.u  <- as.vector(Ainv %*% t(rep(eta2, n)) %*% (u2 - alpha.u2))
    phi.u2      <- as.vector(rnorm(n = q2, mean = mean.phi.u, sd = sqrt(Sigma.phi.u)))
    
    
    # Sample $ \eta1 $
    sum.phi.eta1 <- 0
    for(k in 1:q1)
    {
      sum.phi.eta1 <- sum.phi.eta1 + sum(phi.u1[k] * u1[, k])
    }
    Sigma.eta <- 1/(1/sigma.eta1.square + sum(phi.u1^2)/sigma.u1.square)
    mean.eta  <- Sigma.eta * (eta2/sigma.eta1.square + sum.phi.eta1/sigma.u1.square +
                                phi.t * sum(logt)/sigma.t.square - sum(phi.u1 * alpha.u1)/sigma.u1.square -
                                phi.t * alpha.t/sigma.t.square -
                                (t(beta.t) %*% t(X) %*% rep(phi.t, n))/sigma.t.square)
    eta1      <- as.vector(mvrnorm(n = 1, mu = mean.eta, Sigma = Sigma.eta))
    
    # Sample $ \eta2 $
    sum.phi.eta2 <- 0
    for(l in 1:q2)
    {
      sum.phi.eta2 <- sum.phi.eta2 + sum(phi.u2[l] * u2[, l])
    }
    Sigma.eta <- 1/(1/sigma.eta1.square + 1/sigma.eta2.square + sum(phi.u2^2)/sigma.u2.square)
    mean.eta  <- Sigma.eta * (eta1/sigma.eta1.square + sum.phi.eta2/sigma.u2.square +
                                sum(phi.u2 * alpha.u2)/sigma.u2.square)
    eta2      <- as.vector(mvrnorm(n = 1, mu = mean.eta, Sigma = Sigma.eta))
    
    
    
    
    # Sample $ \sigma_t^2
    E1 <- crossprod(x = logt - alpha.t - X %*% beta.t  - eta1 * phi.t)
    E2 <- crossprod(x = beta.t)
    E3 <- crossprod(x = alpha.t)
    # E4 <- crossprod(x = phi.t)
    sigma.t.square <- 1/stats::rgamma(1, shape = (n + p + 1)/2, scale = 2/(E1 + E2 + E3))
    # if(sigma.t.square.propesed < 10)
    # {
    #   sigma.t.square <- sigma.t.square.propesed
    # }
    
    # # Sample $ \sigma_u1^2 
    # E1 <- crossprod(x = u1 - alpha.u1 - eta1 * phi.u1)
    # E3 <- crossprod(x = alpha.u1)
    # E4 <- crossprod(x = phi.u1)
    # sigma.u1.square <- 1/stats::rgamma(1, (n + 1)/2, scale = 2/(E1 + E3 + E4))
    # 
    # # Sample $ \sigma_u2^2 
    # E1 <- crossprod(x = u2 - alpha.u2 - eta2 * phi.u2)
    # E3 <- crossprod(x = alpha.u2)
    # E4 <- crossprod(x = phi.u2)
    # sigma.u2.square <- 1/stats::rgamma(1, (n + 1)/2, scale = 2/(E1 + E3 + E4))
    
    
    logt.hat <- alpha.t + X %*% beta.t + eta1 * phi.t
    
    loglikelihood <- sum(status * dnorm(logt, mean = alpha.t + X %*% beta.t + eta1 * phi.t, 
                                        sd = sqrt(sigma.t.square), log = TRUE) + 
                           (1 - status) * pnorm(logt, mean = alpha.t + X %*% beta.t + eta1 * phi.t, 
                                                sd = sqrt(sigma.t.square), 
                                                lower.tail = FALSE, log.p = TRUE)) 
    likelihood    <- exp(loglikelihood)
    
    
    # stores the MCMC samples after discarding burnin samples
    if (iter > nburnin)
    {
      beta.tout[, (iter - nburnin)/nthin]   <- beta.t
      alpha.tout[(iter - nburnin)/nthin]    <- alpha.t
      phi.tout[(iter - nburnin)/nthin]      <- phi.t
      alpha.u1out[, (iter - nburnin)/nthin] <- alpha.u1
      alpha.u2out[, (iter - nburnin)/nthin] <- alpha.u2
      phi.u1out[, (iter - nburnin)/nthin]   <- phi.u1
      phi.u2out[, (iter - nburnin)/nthin]   <- phi.u2
      eta1.out[(iter - nburnin)/nthin]      <- eta1
      eta2.out[(iter - nburnin)/nthin]      <- eta2
      
      sigma.t.square.out[(iter - nburnin)/nthin]  <- sigma.t.square
      sigma.u1.square.out[(iter - nburnin)/nthin] <- sigma.u1.square
      sigma.u2.square.out[(iter - nburnin)/nthin] <- sigma.u2.square
      
      logt.out[, (iter - nburnin)/nthin]        <- logt
      logt.hat.out[, (iter - nburnin)/nthin]    <- logt.hat
      loglikelihood.out[(iter - nburnin)/nthin] <- loglikelihood
      likelihood.out[, (iter - nburnin)/nthin]  <- likelihood
    }
  }
  
  
  
  # Posterior Mean
  pMean.beta.t   <- apply(beta.tout, 1, mean)
  pMean.alpha.t  <- mean(alpha.tout)
  pMean.phi.t    <- mean(phi.tout)
  pMean.alpha.u1 <- apply(alpha.u1out, 1, mean)
  pMean.alpha.u2 <- apply(alpha.u2out, 1, mean)
  pMean.phi.u1   <- apply(phi.u1out, 1, mean)
  pMean.phi.u2   <- apply(phi.u2out, 1, mean)
  pMean.eta1     <- mean(eta1.out)
  pMean.eta2     <- mean(eta2.out)
  
  pMean.sigma.t.square  <- median(sigma.t.square.out)
  pMean.sigma.u1.square <- mean(sigma.u1.square.out)
  pMean.sigma.u2.square <- mean(sigma.u2.square.out)
  
  pMean.logt     <- apply(logt.out, 1, mean)
  pMean.logt.hat <- apply(logt.hat.out, 1, mean)
  pLoglikelihood <- mean(loglikelihood.out)
  plikelihood    <- apply(likelihood.out, 1, mean)
  
  
  loglikelihood.posterior <- sum(status * dnorm(pMean.logt, mean = pMean.alpha.t + X %*% pMean.beta.t + 
                                                  pMean.eta1 * pMean.phi.t, 
                                                sd = sqrt(pMean.sigma.t.square), log = TRUE) + 
                                   (1 - status) * pnorm(pMean.logt, mean = pMean.alpha.t + X %*% pMean.beta.t + 
                                                          pMean.eta1 * pMean.phi.t, 
                                                        sd = sqrt(pMean.sigma.t.square), 
                                                        lower.tail = FALSE, log.p = TRUE))
  
  DIC  <- -4 * pLoglikelihood + 2 * loglikelihood.posterior
  lppd <- sum(log(plikelihood))
  WAIC <- -2 * (lppd - 2 * (loglikelihood.posterior - pLoglikelihood))
  # WAIC <- -2 * (lppd - 2 * (lppd - pLoglikelihood))
  
  
  result = list(pMean.beta.t = pMean.beta.t, pMean.alpha.t = pMean.alpha.t, pMean.phi.t = pMean.phi.t,
                pMean.alpha.u1 = pMean.alpha.u1, pMean.alpha.u2 = pMean.alpha.u2, 
                pMean.phi.u1 = pMean.phi.u1, pMean.phi.u2 = pMean.phi.u2, 
                pMean.eta1 = pMean.eta1, pMean.eta2 = pMean.eta2, 
                pMean.sigma.t.square = pMean.sigma.t.square, 
                pMean.sigma.u1.square = pMean.sigma.u1.square, pMean.sigma.u2.square = pMean.sigma.u2.square,
                alpha.t.samples = alpha.tout, phi.t.samples = phi.tout, 
                beta1.t.samples = beta.tout[1, ], beta2.t.samples = beta.tout[2, ],
                beta.t.samples = beta.tout,
                alpha.u1.samples = alpha.u1out, alpha.u2.samples = alpha.u2out,
                phi.u1.samples = phi.u1out, phi.u2.samples = phi.u2out, 
                eta1.samples = eta1.out, eta2.samples = eta2.out,
                sigma.t.square.samples = sigma.t.square.out, 
                sigma.u1.square.samples = sigma.u1.square.out, sigma.u2.square.samples = sigma.u2.square.out,
                pMean.logt.hat = pMean.logt.hat, DIC = DIC, WAIC = WAIC)
  
  return(result)
}