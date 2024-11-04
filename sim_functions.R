
#####

## auxiliary function to estimate regression model and population estimates on response and sample

est.regr.function <- function(stp.pi, resp.y, stp.y, X, pop.gem, N, num.cfs, sigma.est = c('unweighted', 'Hajek', 'Hajekplus')) {
  
  y <- stp.y
  w <- 1 / stp.pi
  r <- (!is.na(resp.y))
  
  stp.L <- lm.wfit(x = X, y = y, w = w)
  resp.L <- lm.wfit(x = X[r, , drop = FALSE], y = y[r], w = w[r])
  
  stp.cf <- stp.L$coefficients
  resp.cf <- resp.L$coefficients
  
  resp.estregr <- resp.cf %*% pop.gem
  stp.estregr <- stp.cf %*% pop.gem
  
  est.pop.gem <- apply(X, 2, function(m) sum(m * w) / N)
  resp.estHj <- resp.cf %*% est.pop.gem * (N / sum(w))
  stp.estHj <- stp.cf %*% est.pop.gem * (N / sum(w))
  
  if (sigma.est == 'unweighted') {
    s2p <- sum(resp.L$residuals^2)/(sum(r) - length(resp.cf))
  } else if (sigma.est == 'Hajek') {
    s2p <- sum(w[r] * resp.L$residuals^2)/sum(w[r])
  } else if (sigma.est == 'Hajekplus') {
    s2p <- (sum(w[r] * resp.L$residuals^2)/sum(w[r])) * sum(r)/(sum(r) - length(resp.cf))
  } else {
    stop('Unknown option sigma.est!')
  }
  
  res <- c(as.numeric(resp.cf),
           as.numeric(resp.estregr),
           as.numeric(resp.estHj),
           as.numeric(stp.cf),
           as.numeric(stp.estregr),
           as.numeric(stp.estHj),
           as.numeric(est.pop.gem),
           s2p)
  names(res) <- c(sprintf('beta%d.est.resp', num.cfs),
                  'ygem.estregr.resp',
                  'ygem.estHj.resp',
                  sprintf('beta%d.est.stp', num.cfs),
                  'ygem.estregr.stp',
                  'ygem.estHj.stp',
                  sprintf('xgem%d.estHT.stp', num.cfs),
                  's2e.resp')
  return(res)
  
}


### auxiliary function to run data augmentation algorithm of Tanner & Wong
### as used by Kim & Yang (2017) (note: Kim & Yang use fpc = TRUE)
impTW <- function(y, X, w, tau.inv, n, it = 25L, fpc = TRUE) {
  
  r <- which(!is.na(y))
  nr <- which(is.na(y))
  
  L <- lm.wfit(x = X[r, ,drop=FALSE], y = y[r], w = w[r])
  cf <- L$coefficients
  s2n <- sum(w[r] * L$residuals^2)/sum(w[r])
  
  for (i in 1:it) {
    
    # I step
    y[nr] <- X[nr, ,drop=FALSE] %*% cf + rnorm(n = length(nr), mean = 0, sd = sqrt(s2n))
    
    # P step
    L <- lm.wfit(x = X, y = y, w = w)
    cf <- L$coefficients
    e <- L$residuals
    s2n <- sum(w*e^2)/sum(w)
    
    if (fpc) {
      V0 <- tau.inv %*% ((t(X) %*% diag(w*(w-1)*e^2) %*% X)) %*% tau.inv
    } else {
      V0 <- tau.inv %*% ((t(X) %*% diag((w^2)*e^2) %*% X)) %*% tau.inv
    }
    
    cf <- cf + mvrnorm(n = 1, mu = rep(0, length(cf)), Sigma = V0)
    
  }
  
  return(y)
  
}


## function to run simulation study for a given scenario

runSimulation <- function(  S                  # number of simulations
                          , S_at_a_time        # number of simulations to process at once
                          , N                  # population size
                          , cfs                # number of regression coefficients
                          , use.intercept = TRUE # include intercept in regression model for y?
                          , m = 100L           # number of multiple imputations
                          , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                          , model.u = NULL
                          , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                          , model.y = '(pop.x - 1)/2 + pop.eps'
                          , model.logit.pi = NULL
                          , model.pi = NULL
                          , model.logit.phi = NULL
                          , model.phi = NULL
                          , sigma.est = 'Hajekplus' # choice of estimator for sigma^2 (see est.regr.function)
                          , draw.u = 'normal'   # how to draw u for imputations ('normal' or 'hotdeck')?
                          , informative = FALSE
                          , runKY = FALSE      # also apply method of Kim & Yang (2017)?
                          , fpcKY = TRUE       # use finite population correction in method of Kim & Yang (2017)?
                          , approxKY = FALSE   # use approximate version of method of Kim & Yang (2017)?
                          , runTradMI = FALSE  # also apply 'traditional multiple imputation'?
                          , writeExcel = FALSE # write output to Excel file?
                          , nameExcel          # name of Excel file with output
                          , seed = NULL        # random seed
) {
  
  start <- Sys.time()
  
  if (!is.null(seed)) set.seed(seed)
  
  # number of simulations to hold in memory at once
  # (assumption: S is a multiple of S_at_a_time)
  S_at_a_time <- min(S, S_at_a_time)
  K <- S / S_at_a_time
  
  # numbering of regression coefficients
  num.cfs <- seq(0, cfs - 1, 1)
  
  
  est.regr <- NULL
  ygem.pop <- NULL
  n <- NULL
  p <- NULL
  estMI <- NULL
  
  # populations: x values (fixed!)
  pop.x <- eval(parse(text = model.x))
  xgem.pop <- mean(pop.x)
  
  if (use.intercept) {
    pop.gem <- c(1, xgem.pop)
  } else {
    pop.gem <- xgem.pop
  }
  
  # populations: u values (fixed!)
  if (!is.null(model.u)) {
    pop.u <- eval(parse(text = model.u))
  }
  
  # populations: response mechanism P(response) = phi
  if (!is.null(model.logit.phi)) {
    pop.phi <- exp(eval(parse(text = model.logit.phi)))
    pop.phi <- pop.phi / (1 + pop.phi)
  } else if (!is.null(model.phi)) {
    pop.phi <- eval(parse(text = model.phi))
  }
  if (length(pop.phi) == 1) pop.phi <- rep(pop.phi, N)
  
  # response indicator
  resp.pop <- (runif(n = N) < pop.phi)
  
  if (!informative) {
    # populations: inclusion probabilities
    if (!is.null(model.logit.pi)) {
      pop.pi <- exp(eval(parse(text = model.logit.pi)))
      pop.pi <- pop.pi / (1 + pop.pi)
    } else if (!is.null(model.pi)) {
      pop.pi <- eval(parse(text = model.pi))
    }
    if (length(pop.pi) == 1) pop.pi <- rep(pop.pi, N)
  }
  
  for (k in 1:K) {
    
    print(sprintf('Running simulation rounds %d to %d of %d...',
                  (k-1)*S_at_a_time + 1, k*S_at_a_time, S))
    
    # populations: disturbance terms
    pop.eps <- sapply(1:S_at_a_time, function(s) {
      eval(parse(text = model.eps))
    })
    
    # populations: y values
    pop.y <- eval(parse(text = model.y))
    
    if (informative) {
      # populations: inclusion probabilities
      if (!is.null(model.logit.pi)) {
        pop.pi <- sapply(1:S_at_a_time, function(s) {
          pis <- exp(eval(parse(text = model.logit.pi)))
          pis <- pis / (1 + pis)
          if (length(pis) == 1) pis <- rep(pis, N)
          return(pis)
        })
      } else if (!is.null(model.pi)) {
        pop.pi <- sapply(1:S_at_a_time, function(s) {
          pis <- eval(parse(text = model.pi))
          if (length(pis) == 1) pis <- rep(pis, N)
          return(pis)
        })
      }
    }
    
    ## determine response values
    resp.pop.y <- sapply(1:S_at_a_time, function(s) {
      ifelse(resp.pop, pop.y[ ,s], NA)
    })
    
    # population means
    ygem.pop0 <- colMeans(pop.y)
    
    ## draw samples
    if (informative) {
      stp <- sapply(1:S_at_a_time, function(s) which(sampling::UPpoisson(pik = pop.pi[ ,s]) > 0))
      stp.y <- sapply(1:S_at_a_time, function(s) pop.y[stp[[s]], s])
      stp.x <- sapply(1:S_at_a_time, function(s) pop.x[stp[[s]]])
      stp.pi <- sapply(1:S_at_a_time, function(s) pop.pi[stp[[s]], s])
    } else {
      stp <- sapply(1:S_at_a_time, function(s) which(sampling::UPpoisson(pik = pop.pi) > 0))
      stp.y <- sapply(1:S_at_a_time, function(s) pop.y[stp[[s]], s])
      stp.x <- sapply(1:S_at_a_time, function(s) pop.x[stp[[s]]])
      stp.pi <- sapply(1:S_at_a_time, function(s) pop.pi[stp[[s]]])
    }
    
    # response values per sample
    resp <- sapply(1:S_at_a_time, function(s) resp.pop[stp[[s]]])
    resp.y <- sapply(1:S_at_a_time, function(s) ifelse(resp[[s]], stp.y[[s]], NA))
    
    
    ## estimate regression model and population estimates on response and sample
    
    est.regr0 <- sapply(1:S_at_a_time, function(s) {
      n <- length(stp.x[[s]])
      
      if (use.intercept) {
        Xn <- cbind(rep(1, n), stp.x[[s]])
      } else {
        Xn <- matrix(stp.x[[s]], nrow = n, ncol = 1)
      }
      
      est.regr.function(stp.pi = stp.pi[[s]], resp.y = resp.y[[s]],
               stp.y = stp.y[[s]], X = Xn, pop.gem = pop.gem,
               N = N, num.cfs = num.cfs, sigma.est = sigma.est)
    })
    
    
    ## perform multiple imputation
    
    cores <- detectCores()
    cl <- makeCluster(cores[1] - 1)
    registerDoParallel(cl)
    
    estMI0 <- foreach(i = 1:S_at_a_time,
                      .combine = rbind,
                      .export = c('impTW'),
                      .packages = c("MASS")) %dopar% {
                        
                        print(i)
                        if (!is.null(seed)) set.seed(seed + (k-1)*S_at_a_time + i)
                        
                        ygem.est.MI.m <- rep(NA_real_, m)
                        Vbinnen <- rep(NA_real_, m)
                        
                        n <- length(stp.y[[i]])
                        r <- resp[[i]]
                        w <- 1/stp.pi[[i]]
                        
                        if (use.intercept) {
                          Xn <- cbind(rep(1, n), stp.x[[i]])
                        } else {
                          Xn <- matrix(stp.x[[i]], nrow = n, ncol = 1)
                        }
                        Xq <- Xn[!r, , drop = FALSE]
                        Xr <- Xn[r, , drop = FALSE]
                        
                        Vn <- solve(t(Xn) %*% (Xn / stp.pi[[i]]))
                        Vq_inv <- t(Xq) %*% (Xq / stp.pi[[i]][!r])
                        Vq <- solve(Vq_inv)
                        Vq_unw_inv <- t(Xq) %*% Xq
                        Vr <- solve(t(Xr) %*% (Xr / stp.pi[[i]][r]))
                        
                        Kq <- Vn %*% t(Xq / stp.pi[[i]][!r])
                        fq <- N/(w[!r]*n)
                        s2e.p <- est.regr0['s2e.resp', i]
                        Aq <- t(Xq) %*% (Xq / (stp.pi[[i]][!r])) * (N / n) * s2e.p
                        
                        resid <- lm.wfit(x = Xr, y = stp.y[[i]][r], w = w[r])$residuals
                        if (informative) {
                          Ar <- t(Xr) %*% (Xr * resid^2 / (stp.pi[[i]][r]^2))
                        } else {
                          Ar <- t(Xr) %*% (Xr / (stp.pi[[i]][r]^2)) * s2e.p
                        }
                        
                        bp <- est.regr0[sprintf('beta%d.est.resp', num.cfs), i, drop = FALSE]
                        xgemn <- est.regr0[sprintf('xgem%d.estHT.stp', num.cfs), i, drop = FALSE]
                        
                        BRG <- Vr %*% Ar %*% Vr - Vn %*% (Ar + Aq) %*% Vn +
                          s2e.p * Vn %*% ((N/n)*Vq_unw_inv - Vq_inv) %*% Vn
                        if (informative) {
                          lambda <- as.numeric(t(pop.gem) %*% ((2/(sum(w))) * (Vr - Vn) %*% (t(Xr) %*% ((Xr - rep(1, nrow(Xr)) %*% t(pop.gem)) * resid * w[r]^2))) %*% bp) /
                            as.numeric(t(pop.gem) %*% BRG %*% pop.gem)
                        } else {
                          lambda <- 0
                        }
                        BRGlambda <- BRG * (1 + lambda)
                        
                        Wp <- Vq %*% solve(Vn) %*% (BRGlambda / s2e.p - Kq %*% (t(Kq) * fq)) %*% solve(Vn) %*% Vq
                        
                        # check whether Wp is a positive semidefinite matrix
                        suitable <- all(eigen(Wp)$values >= 0)
                        if (suitable) {
                          bn <- matrix(bp, nrow = cfs, ncol = m) +
                            t(mvrnorm(n = m, mu = rep(0, cfs), Sigma = Wp * s2e.p))
                          gamma <- 1
                        } else { # w = 0
                          bn <- matrix(bp, nrow = cfs, ncol = m)
                          
                          if (informative) {
                            gamma <- (as.numeric(t(pop.gem) %*% (Vr %*% Ar %*% Vr - Vn %*% (Ar + s2e.p * Vq_inv) %*% Vn) %*% pop.gem) + as.numeric(t(pop.gem) %*% ((2/(sum(w))) * (Vr - Vn) %*% (t(Xr) %*% ((Xr - rep(1, nrow(Xr)) %*% t(pop.gem)) * resid * w[r]^2))) %*% bp)) /
                              as.numeric(t(pop.gem) %*% (Vn %*% (2*Aq - s2e.p * (N/n) * Vq_unw_inv) %*% Vn) %*% pop.gem)
                          } else {
                            gamma <- as.numeric(t(pop.gem) %*% (Vr %*% Ar %*% Vr - Vn %*% (Ar + s2e.p * Vq_inv) %*% Vn) %*% pop.gem) /
                              as.numeric(t(pop.gem) %*% (Vn %*% (2*Aq - s2e.p * (N/n) * Vq_unw_inv) %*% Vn) %*% pop.gem)
                          }
                        }
                        
                        if (gamma >= 0) {
                          y <- resp.y[[i]]
                          
                          for (j in 1:m) {
                            # generate imputations
                            if (draw.u == 'normal') {
                              y[!r] <- Xq %*% bn[ ,j] + rnorm(n = sum(!r), mean = 0, sd = sqrt(s2e.p)) * sqrt(fq * gamma)
                            } else if (draw.u == 'hotdeck') {
                              y[!r] <- Xq %*% bn[ ,j] + sample(x = resid, size = sum(!r), replace = TRUE, prob = w[r]/sum(w[r])) * sqrt(fq * gamma)
                            }
                            
                            # estimate population mean from imputed sample
                            # L <- lm.wfit(x = Xn, y = y, w = w)
                            # ygem.est.MI.m[j] <- L$coefficients %*% xgemn * (N / sum(w))
                            ygem.est.MI.m[j] <- weighted.mean(y, w = w)
                            
                            # compute contribution of within variance
                            Vbinnen[j] <- sum((w^2 - w) * (y - ygem.est.MI.m[j])^2) / (sum(w)^2)
                          }
                          
                          ygem.est.MI <- mean(ygem.est.MI.m)
                          U.MI <- mean(Vbinnen)
                          B.MI <- (1 + 1/m) * var(ygem.est.MI.m)
                        } else { # no MI possible
                          ygem.est.MI <- NA
                          U.MI <- NA
                          B.MI <- NA
                          lambda <- NA
                          gamma <- NA
                        }
                        
                        returnlist <- c(ygem.est.MI = ygem.est.MI,
                                        V.MI = U.MI + B.MI,
                                        U.MI = U.MI,
                                        B.MI = B.MI,
                                        s2e.p = ifelse(is.na(gamma), NA, s2e.p[[1]]),
                                        suitable = suitable,
                                        lambda = lambda,
                                        gamma = gamma)
                        
                        if (runKY) { # also apply method of Kim & Yang (2017)
                          
                          ygem.est.MI.KY.m <- rep(NA_real_, m)
                          Vbinnen.KY <- rep(NA_real_, m)
                          
                          if (approxKY) { # use fast approximation to Kim & Yang's method
                            
                            if (fpcKY) {
                              Wp.KY <- (Vr %*% Ar %*% Vr)/s2e.p - Vn
                            } else {
                              Wp.KY <- (Vr %*% Ar %*% Vr)/s2e.p
                            }
                            
                            # check whether Wp is a positive semidefinite matrix
                            suitable.KY <- all(eigen(Wp.KY)$values >= 0)
                            if (suitable.KY) {
                              bn.KY <- matrix(bp, nrow = cfs, ncol = m) +
                                t(mvrnorm(n = m, mu = rep(0, cfs), Sigma = Wp.KY * s2e.p))
                            } else { # w = 0
                              bn.KY <- matrix(bp, nrow = cfs, ncol = m)
                            }
                            
                            yimp.KY <- resp.y[[i]]
                            
                            for (j in 1:m) {
                              # generate imputations
                              if (draw.u == 'normal') {
                                yimp.KY[!r] <- Xq %*% bn.KY[ ,j] + rnorm(n = sum(!r), mean = 0, sd = sqrt(s2e.p))
                              } else if (draw.u == 'hotdeck') {
                                yimp.KY[!r] <- Xq %*% bn.KY[ ,j] + sample(x = resid, size = sum(!r), replace = TRUE, prob = w[r]/sum(w[r]))
                              }
                              
                              # estimate population mean from imputed sample
                              ygem.est.MI.KY.m[j] <- weighted.mean(yimp.KY, w = w)
                              
                              # compute contribution of within variance
                              Vbinnen.KY[j] <- sum((w^2 - w) * (yimp.KY - ygem.est.MI.KY.m[j])^2) / (sum(w)^2)
                            }
                            
                          } else { # use Kim & Yang's method in its original form
                            
                            tau.inv <- solve(t(Xn) %*% diag(w) %*% Xn)
                            suitable.KY <- TRUE
                            
                            for (j in 1:m) {
                              # generate imputations
                              yimp.KY <- impTW(y = resp.y[[i]],
                                               X = Xn,
                                               w = w,
                                               tau.inv = tau.inv,
                                               n = n, it = 25L, fpc = fpcKY)
                              
                              # estimate population mean from imputed sample
                              ygem.est.MI.KY.m[j] <- weighted.mean(yimp.KY, w = w)
                              
                              # compute contribution of within variance
                              Vbinnen.KY[j] <- sum((w^2 - w) * (yimp.KY - ygem.est.MI.KY.m[j])^2) / (sum(w)^2)
                            }
                            
                          }
                          
                          ygem.est.MI.KY = mean(ygem.est.MI.KY.m)
                          U.MI.KY <- mean(Vbinnen.KY)
                          B.MI.KY <- (1 + 1/m) * var(ygem.est.MI.KY.m)
                          
                          returnlist <- c(returnlist,
                                          ygem.est.MI.KY = ygem.est.MI.KY,
                                          V.MI.KY = U.MI.KY + B.MI.KY,
                                          U.MI.KY = U.MI.KY,
                                          B.MI.KY = B.MI.KY,
                                          suitable.KY = suitable.KY)
                          
                        }
                        
                        if (runTradMI) { # also apply 'traditional multiple imputation'
                          
                          ygem.est.MI.trad.m <- rep(NA_real_, m)
                          Vbinnen.trad <- rep(NA_real_, m)
                          
                          Xn.trad <- cbind(Xn, w)
                          Xq.trad <- Xn.trad[!r, , drop = FALSE]
                          Xr.trad <- Xn.trad[r, , drop = FALSE]
                          
                          resp.L.trad <- lm.fit(x = Xr.trad, y = stp.y[[i]][r])
                          bp.trad <- resp.L.trad$coefficients
                          resid.trad <- resp.L.trad$residuals
                          s2p.trad <- var(resid.trad)*(sum(r)-1)/(sum(r) - length(bp.trad))
                          Wp.trad <- solve(t(Xr.trad) %*% Xr.trad)
                          
                          bn.trad <- matrix(bp.trad, nrow = length(bp.trad), ncol = m) +
                            t(mvrnorm(n = m, mu = rep(0, length(bp.trad)), Sigma = Wp.trad * s2p.trad))
                          
                          yimp.trad <- resp.y[[i]]
                          
                          for (j in 1:m) {
                            # generate imputations
                            if (draw.u == 'normal') {
                              yimp.trad[!r] <- Xq.trad %*% bn.trad[ ,j] + rnorm(n = sum(!r), mean = 0, sd = sqrt(s2p.trad))
                            } else if (draw.u == 'hotdeck') {
                              yimp.trad[!r] <- Xq.trad %*% bn.trad[ ,j] + sample(x = resid.trad, size = sum(!r), replace = TRUE, prob = w[r]/sum(w[r]))
                            }
                            
                            # estimate population mean from imputed sample
                            ygem.est.MI.trad.m[j] <- weighted.mean(yimp.trad, w = w)
                            
                            # compute contribution of within variance
                            Vbinnen.trad[j] <- sum((w^2 - w) * (yimp.trad - ygem.est.MI.trad.m[j])^2) / (sum(w)^2)
                          }
                          
                          ygem.est.MI.trad = mean(ygem.est.MI.trad.m)
                          U.MI.trad <- mean(Vbinnen.trad)
                          B.MI.trad <- (1 + 1/m) * var(ygem.est.MI.trad.m)
                          
                          returnlist <- c(returnlist,
                                          ygem.est.MI.trad = ygem.est.MI.trad,
                                          V.MI.trad = U.MI.trad + B.MI.trad,
                                          U.MI.trad = U.MI.trad,
                                          B.MI.trad = B.MI.trad,
                                          s2p.trad = s2p.trad)
                          
                        }
                        
                        return(returnlist)
                      }
    
    stopCluster(cl)
    
    est.regr <- cbind(est.regr, est.regr0)
    
    estMI <- rbind(estMI, estMI0)
    
    ygem.pop <- c(ygem.pop, ygem.pop0)
    
    n <- c(n, sapply(stp, length))
    p <- c(p, sapply(resp.y, function(s) sum(!is.na(s))))
    
    rm(pop.eps, pop.y, resp.pop.y,
       resp, resp.y,
       stp, stp.pi, stp.x, stp.y,
       est.regr0, estMI0, ygem.pop0)
    
    if (informative) {
      rm(pop.pi)
    }
    
  }
  
  
  ##### compute simulated properties of estimators
  
  # estimated population mean based on Hajek estimator (MI)
  bias.ygem.est.MI.Hj <- mean(estMI[ ,'ygem.est.MI'] - ygem.pop, na.rm = TRUE)
  var.ygem.est.MI.Hj <- var(estMI[ ,'ygem.est.MI'] - ygem.pop, na.rm = TRUE)
  
  resMI0 <- data.frame(method = 'MI + Hajek estimator',
                       bias_ygem = bias.ygem.est.MI.Hj,
                       var_ygem = var.ygem.est.MI.Hj,
                       V_MI = mean(estMI[ , 'V.MI'], na.rm = TRUE),
                       U_MI = mean(estMI[ , 'U.MI'], na.rm = TRUE),
                       B_MI = mean(estMI[ , 'B.MI'], na.rm = TRUE),
                       s2e.p = mean(estMI[ , 's2e.p'], na.rm = TRUE),
                       lambda_gem = mean(estMI[ , 'lambda'], na.rm = TRUE),
                       lambda_min = min(estMI[ , 'lambda'], na.rm = TRUE),
                       lambda_max = max(estMI[ , 'lambda'], na.rm = TRUE),
                       gamma_gem = mean(estMI[ , 'gamma'], na.rm = TRUE),
                       gamma_min = min(estMI[ , 'gamma'], na.rm = TRUE),
                       gamma_max = max(estMI[ , 'gamma'], na.rm = TRUE),
                       gamma_used = sum(estMI[ , 'suitable'] == 0, na.rm = FALSE),
                       gamma_neg = sum(is.na(estMI[ , 'gamma'])))
  
  if (runKY) {
    
    bias.ygem.est.MI.Hj.KY <- mean(estMI[ ,'ygem.est.MI.KY'] - ygem.pop, na.rm = TRUE)
    var.ygem.est.MI.Hj.KY <- var(estMI[ ,'ygem.est.MI.KY'] - ygem.pop, na.rm = TRUE)
    
    resMI0 <- rbind(resMI0,
                    data.frame(method = ifelse(approxKY,
                                               'MI + Hajek estimator, Kim & Yang (approx)',
                                               'MI + Hajek estimator, Kim & Yang'),
                               bias_ygem = bias.ygem.est.MI.Hj.KY,
                               var_ygem = var.ygem.est.MI.Hj.KY,
                               V_MI = mean(estMI[ , 'V.MI.KY'], na.rm = TRUE),
                               U_MI = mean(estMI[ , 'U.MI.KY'], na.rm = TRUE),
                               B_MI = mean(estMI[ , 'B.MI.KY'], na.rm = TRUE),
                               s2e.p = NA,
                               lambda_gem = NA,
                               lambda_min = NA,
                               lambda_max = NA,
                               gamma_gem = NA,
                               gamma_min = NA,
                               gamma_max = NA,
                               gamma_used = sum(estMI[ , 'suitable.KY'] == 0, na.rm = FALSE),
                               gamma_neg = NA))
    
  }
  
  if (runTradMI) {
    
    bias.ygem.est.MI.Hj.trad <- mean(estMI[ ,'ygem.est.MI.trad'] - ygem.pop, na.rm = TRUE)
    var.ygem.est.MI.Hj.trad <- var(estMI[ ,'ygem.est.MI.trad'] - ygem.pop, na.rm = TRUE)
    
    resMI0 <- rbind(resMI0,
                    data.frame(method = 'MI + Hajek estimator, traditional MI',
                               bias_ygem = bias.ygem.est.MI.Hj.trad,
                               var_ygem = var.ygem.est.MI.Hj.trad,
                               V_MI = mean(estMI[ , 'V.MI.trad'], na.rm = TRUE),
                               U_MI = mean(estMI[ , 'U.MI.trad'], na.rm = TRUE),
                               B_MI = mean(estMI[ , 'B.MI.trad'], na.rm = TRUE),
                               s2e.p = mean(estMI[ , 's2p.trad'], na.rm = TRUE),
                               lambda_gem = NA,
                               lambda_min = NA,
                               lambda_max = NA,
                               gamma_gem = NA,
                               gamma_min = NA,
                               gamma_max = NA,
                               gamma_used = NA,
                               gamma_neg = NA))
    
  }
  
  resMI <- t(resMI0[,-1])
  colnames(resMI) <- resMI0[,1]
  
  end <- Sys.time()
  duration <- end - start
  
  
  if (writeExcel) { # write results to Excel
    
    wb <- createWorkbook()
    
    addWorksheet(wb, 'results')
    
    writeData(wb, 'results',
              'N', xy = c(1, 1))
    writeData(wb, 'results',
              N, xy = c(2, 1))
    
    writeData(wb, 'results',
              'x', xy = c(1, 2))
    writeData(wb, 'results',
              model.x, xy = c(2, 2))
    
    writeData(wb, 'results',
              'eps', xy = c(1, 3))
    writeData(wb, 'results',
              model.eps, xy = c(2, 3))
    
    if (!is.null(model.u)) {
      writeData(wb, 'results',
                'u', xy = c(1, 4))
      writeData(wb, 'results',
                model.u, xy = c(2, 4))
    }
    
    writeData(wb, 'results',
              'y', xy = c(1, 5))
    writeData(wb, 'results',
              model.y, xy = c(2, 5))
    
    if (!is.null(model.logit.pi)) {
      writeData(wb, 'results',
                'logit(pi)', xy = c(1, 6))
      writeData(wb, 'results',
                model.logit.pi, xy = c(2, 6))
    } else if (!is.null(model.pi)) {
      writeData(wb, 'results',
                'pi', xy = c(1, 6))
      writeData(wb, 'results',
                model.pi, xy = c(2, 6))
    }
    
    if (!is.null(model.logit.phi)) {
      writeData(wb, 'results',
                'logit(phi)', xy = c(1, 7))
      writeData(wb, 'results',
                model.logit.phi, xy = c(2, 7))
    } else if (!is.null(model.phi)) {
      writeData(wb, 'results',
                'phi', xy = c(1, 7))
      writeData(wb, 'results',
                model.phi, xy = c(2, 7))
    }
    
    writeData(wb, 'results',
              'simulations', xy = c(1, 9))
    writeData(wb, 'results',
              S, xy = c(2, 9))
    writeData(wb, 'results',
              'imputations', xy = c(1, 10))
    writeData(wb, 'results',
              m, xy = c(2, 10))
    writeData(wb, 'results',
              'running time', xy = c(1, 11))
    writeData(wb, 'results',
              format(round(duration, digits = 1)), xy = c(2, 11))
    
    writeData(wb, 'results',
              'mean n', xy = c(1, 13))
    writeData(wb, 'results',
              mean(n), xy = c(2, 13))
    writeData(wb, 'results',
              'st.dev. n', xy = c(1, 14))
    writeData(wb, 'results',
              sd(n), xy = c(2, 14))
    writeData(wb, 'results',
              'mean p', xy = c(1, 15))
    writeData(wb, 'results',
              mean(p), xy = c(2, 15))
    writeData(wb, 'results',
              'st.dev. p', xy = c(1, 16))
    writeData(wb, 'results',
              sd(p), xy = c(2, 16))
    
    
    writeData(wb, 'results',
              resMI, xy = c(4, 4), rowNames = TRUE)
    
    addStyle(wb, 'results',
             createStyle(numFmt = '0.000000'),
             rows = 5:(3+nrow(resMI)-1), cols = 5:(5+as.integer(runKY)+as.integer(runTradMI)),
             gridExpand = TRUE)
    addStyle(wb, 'results',
             createStyle(numFmt = '0'),
             rows = (3+nrow(resMI)) + (0:1), cols = 5:(5+as.integer(runKY)+as.integer(runTradMI)),
             gridExpand = TRUE)
    
    addStyle(wb, 'results',
             createStyle(halign = 'right'),
             rows = 1:16, cols = 2,
             gridExpand = TRUE)
    addStyle(wb, 'results',
             createStyle(numFmt = '0.0', halign = 'right'),
             rows = c(11, 13:16), cols = 2,
             gridExpand = TRUE)
    
    
    setColWidths(wb, 'results',
                 1:6, widths = "auto")
    
    
    saveWorkbook(wb, file = paste0(nameExcel, '.xlsx'),
                 overwrite = TRUE)
  }
  
  return(resMI)
  
}
