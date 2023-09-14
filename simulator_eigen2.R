# load class methods
source("C:/Users/enriq/Desktop/UNIVERSIDAD/Máster en Big Data/TFM/Code/datos_multivariantes_simulados/simulation_model_workshop_Belen/funData-master/R/funDataClass.R")
source("C:/Users/enriq/Desktop/UNIVERSIDAD/Máster en Big Data/TFM/Code/datos_multivariantes_simulados/simulation_model_workshop_Belen/funData-master/R/funDataMethods.R")

# load engines
efFourier <- function(argvals, M, linear = FALSE)
{
  Phi <- matrix(NA, nrow = M, ncol = length(argvals))
  
  Phi[1,] <- sqrt(1 / diff(range(argvals)))
  
  if(M == 1)
    return(funData(argvals, Phi))
  
  for(m in 2:M)
  {
    if(m %% 2 == 0) # m even
      Phi[m, ] <- sqrt(2 / diff(range(argvals))) * cos((m %/% 2) * (2*pi * (argvals - min(argvals)) / diff(range(argvals)) - pi))
    else # m odd
      Phi[m, ] <- sqrt(2 / diff(range(argvals))) * sin((m %/% 2) * (2*pi * (argvals - min(argvals)) / diff(range(argvals)) - pi))
  }
  
  if(linear) # overwrite Phi[M, ], add linear function and orthonormalize (Gram-Schmidt)
  {
    if(any(range(argvals) != c(0,1)))
      stop("efFourier, option linear: not yet implemented for argvals != [0,1]!")
    
    # orthogonalize (exact function)
    Phi[M, ] <- argvals - 1/2  + rowSums(apply(matrix(seq_len(((M-1) %/% 2))), 1, function(k) (-1)^k / (pi*k) * sin(k * (2*pi*argvals - pi))))
    # normalize
    Phi[M, ] <-  Phi[M, ] / sqrt(1/3 - 1/4 - 1 / (2*pi^2)* sum( 1 / (seq_len(((M-1) %/% 2 )))^2 ))
  }
  
  return(funData(argvals, Phi))
}
eFun <- function(argvals, M, ignoreDeg = NULL, type)
{
  if(! all(is.numeric(argvals), length(argvals) > 0))
    stop("Parameter 'argvals' must be numeric.")
  
  if(! all(is.numeric(M), length(M) == 1, M > 0))
    stop("Parameter 'M' must be passed as a positive number.") 
  
  if(!(is.null(ignoreDeg ) | all(is.numeric(ignoreDeg), ignoreDeg > 0)))
    stop("Parameter 'ignoreDeg' must be either NULL or a vector of positive numbers.") 
  
  if(! all(is.character(type), length(type) == 1))
    stop("Parameter 'type' must be passed as a string.")
  
  
  
  ret <- switch(type,
                Poly = efPoly(argvals, M),
                PolyHigh = {
                  if(is.null(ignoreDeg ))
                    stop("eFun, type = PolyHigh: specify ignoreDeg !")
                  
                  efPoly(argvals, M + length(ignoreDeg))[-ignoreDeg]
                },
                Fourier = efFourier(argvals, M, linear = FALSE),
                FourierLin = efFourier(argvals, M, linear = TRUE),
                Wiener = efWiener(argvals, M),
                stop("Choose either Poly, PolyHigh, Fourier, FourierLin or Wiener"))
  return(ret)
}
simMultiSplit <- function(argvals, M, eFunType, ignoreDeg = NULL, eValType, N)
{
  # consistency check
  if(any( c(length(M), length(eFunType), length(eValType)) != 1) )
    stop("argvals, M, eFunType, eValType must all be of length 1!")
  
  # number of elements
  p <- length(argvals)
  
  # "rearrange" argvalss
  x <- vector("list", length = length(argvals))
  splitVals <- rep(NA, length(argvals) + 1)
  
  x[[1]] <- unlist(argvals[[1]]) # convert to vector, if argvals[[1]] is a list
  splitVals[1:2] <- c(0, length(x[[1]]))
  
  for(i in 2:p)
  {
    x[[i]] <- unlist(argvals[[i]]) # convert to vector, if argvals[[i]] is a list
    x[[i]] <- argvals[[i]] - min(argvals[[i]]) + max(x[[i-1]])
    splitVals[i+1] <- splitVals[i]+length(x[[i]])
  }
  
  # generate "big" orthonormal system
  f <-  eFun(unlist(x), M, ignoreDeg = ignoreDeg, type = eFunType)
  
  # sample sign randomly
  s <- sample(c(-1,1), p, 0.5)
  
  # result object
  trueFuns  <- vector("list", p)
  
  for(j in seq_len(p))
    trueFuns[[j]] <- funData(argvals[[j]],  s[j] * f@X[,(1 + splitVals[j]):splitVals[j+1]])
  
  return(multiFunData(trueFuns))
}
eVal <- function(M, type)
{
  if(! all(is.numeric(M), length(M) == 1, M > 0))
    stop("Parameter 'M' must be passed as a positive number.") 
  
  if(! all(is.character(type), length(type) == 1))
    stop("Parameter 'type' must be passed as a string.")
  
  ret <- switch(type,
                linear = ((M+1) - (seq_len(M))) / M,
                exponential = exp(-((seq_len(M)-1) / 2)),
                wiener = 1/(pi/2 * (2 * (seq_len(M)) - 1))^2,
                stop("Choose either linear, exponential or wiener"))
  return(ret)
}
# load interface

simMultiFunData <- function(type, argvals, M,
                            eFunType, ignoreDeg = NULL,
                            eValType, N, 
                            mufs = NULL)
{
  if(! all(is.character(type), length(type) == 1))
    stop("Parameter 'type' must be passed as a string.")
  
  if(! (is.list(argvals) & all(is.numeric(unlist(argvals)))) )
    stop("Parameter 'argvals' must be passed as a list of numerics.")
  
  if(! all(is.numeric(unlist(M))))
    stop("Parameter 'M' must contain only numerics.") 
  
  if(! all(is.character(unlist(eFunType))))
    stop("Parameter 'eFunType' must contain only strings.")
  
  if(!(is.null(ignoreDeg ) | all(is.numeric(ignoreDeg), ignoreDeg > 0)))
    stop("Parameter 'ignoreDeg' must be either NULL or a vector of positive numbers.") 
  
  if(! all(is.character(eValType), length(eValType) == 1))
    stop("Parameter 'eValType' must be passed as a string.")
  
  if(! all(is.numeric(N), length(N) == 1, N > 0))
    stop("Parameter 'N' must be passed as a positive number.") 
  
  # generate eigenfunctions
  trueFuns <- switch(type,
                     split = simMultiSplit(argvals, M, eFunType, ignoreDeg, eValType, N),
                     weighted = simMultiWeight(argvals, M, eFunType, ignoreDeg, eValType, N),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data.")
  )
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(Mtotal, eValType)
  scores <- t(replicate(N, stats::rnorm(Mtotal, sd = sqrt(eVal(Mtotal, eValType)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    # add the mean here:
    if(!is.null(mufs)){
      X = X + rep(mufs[[j]](argvals[[j]]), rep(nrow(X), ncol(X)))
    }
    
    if(N == 1)
      dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  
  return(list(simData = multiFunData(simData),
              trueFuns = trueFuns,
              trueVals = trueVals))
}

addError <- function(funDataObject, sd){
  if(! all(is.numeric(sd), length(sd) == 1, sd > 0))
    stop("Parameter 'sd' must be passed as a positive number.") 
  
  ME <- array(stats::rnorm(prod(dim(funDataObject@X)), mean = 0, sd = sd),
              dim(funDataObject@X))
  
  return(funDataObject + funData(funDataObject@argvals, ME))
}

addErrorMult <- function(funDataObject, sd){
  return(multiFunData(mapply(function(dat, sd){addError(dat, sd)},
                             funDataObject, sd)))
}


# model0, no outliers, 2 variants (different mean functions)

multimodel0_0 <- function(no_functions = 100,
                          argvals = list(seq(0,1, length.out = 50),
                                         seq(0,1, length.out = 50),
                                         seq(0, 1, length.out = 50)),
                          basisfunction = "Fourier",
                          eigenfntype = "split",
                          nbasis = 9,
                          eigenvalue_decay = "linear", 
                          mean_fns = list(mu1 = function(t){return(4*t)},
                                          mu2 = function(t){return(30*t* ((1-t)^(3/2)) )},
                                          mu3 = function(t){return(5 *((t-1)^2))}), 
                          error_params = c(.1, .3)){
  
  ssm1 <- simMultiFunData(N = no_functions,
                          argvals = argvals,
                          eFunType = basisfunction,
                          eValType = eigenvalue_decay,
                          type = eigenfntype, 
                          M = nbasis,
                          mufs = mean_fns)
  
  # now add error
  ssm1$simData <- addErrorMult(ssm1$simData,
                               sd = runif(length(argvals),
                                          error_params[1],
                                          error_params[2]))
  
  return(ssm1$simData)
}


multimodel0_1 <- function(no_functions = 100,
                          argvals = list(seq(0,1, length.out = 50),
                                         seq(0,1, length.out = 50),
                                         seq(0, 1, length.out = 50)),
                          basisfunction = "Fourier",
                          eigenfntype = "split",
                          nbasis = 9,
                          eigenvalue_decay = "linear", 
                          mean_fns = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                          mu2 = function(t){return(5 * cos(2*pi*t))},
                                          mu3 = function(t){return(5 * ((t-1)^2) )}), 
                          error_params = c(.1, .3)){
  
  ssm1 <- simMultiFunData(N = no_functions,
                          argvals = argvals,
                          eFunType = basisfunction,
                          eValType = eigenvalue_decay,
                          type = eigenfntype, 
                          M = nbasis,
                          mufs = mean_fns)
  
  # now add error
  ssm1$simData <- addErrorMult(ssm1$simData,
                               sd = runif(length(argvals),
                                          error_params[1],
                                          error_params[2]))
  
  return(ssm1$simData)
}


# Model 1: todo
multimodel1 <- function(no_functions = 100,
                          argvals = list(seq(0,1, length.out = 50),
                                         seq(0,1, length.out = 50),
                                         seq(0, 1, length.out = 50)),
                          basisfunction = "Fourier",
                          eigenfntype = "split",
                          nbasis = 9,
                          eigenvalue_decay = "linear",
                          mean_fns = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                          mu2 = function(t){return(5 * cos(2*pi*t))},
                                          mu3 = function(t){return(5 * ((t-1)^2) )}),
                          error_params = c(.1, .3),
                          kprob = 0.5,
                          n_outliers = 20,
                        specialp = c())


{
  
  # generate eigenfunctions
  trueFuns <- switch(eigenfntype,
                     split = simMultiSplit(argvals = argvals,
                                           M = nbasis,
                                           eFunType = basisfunction,
                                           ignoreDeg = NULL,
                                           eValType = eigenvalue_decay,
                                           N = no_functions),
                     weighted = simMultiWeight(argvals = argvals,
                                               M = nbasis,
                                               eFunType = basisfunction,
                                               ignoreDeg = NULL,
                                               eValType = eigenvalue_decay,
                                               N = no_functions),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
  scores <- t(replicate(no_functions,
                        stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
                                                            type = eigenvalue_decay)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    X <- scores %*% trueFuns[[j]]@X
    # add the mean here:
    if(!is.null(mean_fns)){
      #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
      X = X + rep(mean_fns[[j]](argvals[[j]]), each = no_functions)
    }
    if(length(specialp) == 0){
      if(n_outliers > 0){
        qcoeff <- rbinom(n_outliers, 1, kprob) 
        qcoeff[qcoeff == 0] <- -1 
        X[1:n_outliers, ] <- X[1:n_outliers, ] + 8*qcoeff 
      }
    }else{
      if(!(j %in% specialp)){
        if(n_outliers > 0){
          qcoeff <- rbinom(n_outliers, 1, kprob) 
          qcoeff[qcoeff == 0] <- -1 
          X[1:n_outliers, ] <- X[1:n_outliers, ] + 8*qcoeff 
        }
      }
    }
    
    if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  simData = multiFunData(simData)
  simData <- addErrorMult(simData,
                          sd = runif(length(argvals),
                                     error_params[1],
                                     error_params[2]))
  
  return(simData)
}


# multimodel1_1 <- function(no_functions = 100,
#                           argvals = list(seq(0,1, length.out = 50),
#                                          seq(0,1, length.out = 50),
#                                          seq(0, 1, length.out = 50)),
#                           basisfunction = "Fourier",
#                           eigenfntype = "split",
#                           nbasis = 9,
#                           eigenvalue_decay = "linear",
#                           mean_fns = list(mu1 = function(t){return(5 * sin(2*pi*t))},
#                                           mu2 = function(t){return(5 * cos(2*pi*t))},
#                                           mu3 = function(t){return(5 * ((t-1)^2) )}),
#                           error_params = c(.1, .3),
#                           kprob = 0.5,
#                           n_outliers = 20,
#                           specialp = c(1))
# 
# 
# {
#   
#   # generate eigenfunctions
#   trueFuns <- switch(eigenfntype,
#                      split = simMultiSplit(argvals = argvals,
#                                            M = nbasis,
#                                            eFunType = basisfunction,
#                                            ignoreDeg = NULL,
#                                            eValType = eigenvalue_decay,
#                                            N = no_functions),
#                      weighted = simMultiWeight(argvals = argvals,
#                                                M = nbasis,
#                                                eFunType = basisfunction,
#                                                ignoreDeg = NULL,
#                                                eValType = eigenvalue_decay,
#                                                N = no_functions),
#                      stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
#   
#   # number of eigenfunctions generated
#   Mtotal <- nObs(trueFuns)
#   
#   # number of elements in multivariate functional basis
#   p <- length(trueFuns)
#   
#   # generate eigenvalues and scores
#   trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
#   scores <- t(replicate(no_functions,
#                         stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
#                                                             type = eigenvalue_decay)))))
#   
#   # generate individual observations
#   simData  <- vector("list", p)
#   
#   for(j in seq_len(p))
#   {
#     #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
#     X <- scores %*% trueFuns[[j]]@X
#     # add the mean here:
#     if(!is.null(mean_fns)){
#       #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
#       X = X + rep(mean_fns[[j]](argvals[[j]]), each = no_functions)
#     }
#     
#     if(!(j %in% specialp)){
#       if(n_outliers > 0){
#         qcoeff <- rbinom(n_outliers, 1, kprob) 
#         qcoeff[qcoeff == 0] <- -1 
#         X[1:n_outliers, ] <- X[1:n_outliers, ] + 8*qcoeff 
#       }
#     }
#     if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
#     
#     simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
#   } 
#   simData = multiFunData(simData)
#   simData <- addErrorMult(simData,
#                           sd = runif(length(argvals),
#                                      error_params[1],
#                                      error_params[2]))
#   
#   return(simData)
# }
# Model 2 Also magnitude model:
multimodel2 <- function(no_functions = 100,
                          argvals = list(seq(0,1, length.out = 50),
                                         seq(0,1, length.out = 50),
                                         seq(0, 1, length.out = 50)),
                          basisfunction = "Fourier",
                          eigenfntype = "split",
                          nbasis = 9,
                          eigenvalue_decay = "linear",
                          mean_fns = list(mu1 = function(t){return(4*t)},
                                          mu2 = function(t){return(30*t* ((1-t)^(3/2)) )},
                                          mu3 = function(t){return(5 *((t-1)^2))}),
                          error_params = c(.1, .3),
                          kprob = 0.5,
                          n_outliers = 20, 
                          specialp = c())
  
 
  {
  
  # generate eigenfunctions
  trueFuns <- switch(eigenfntype,
                     split = simMultiSplit(argvals = argvals,
                                           M = nbasis,
                                           eFunType = basisfunction,
                                           ignoreDeg = NULL,
                                           eValType = eigenvalue_decay,
                                           N = no_functions),
                     weighted = simMultiWeight(argvals = argvals,
                                               M = nbasis,
                                               eFunType = basisfunction,
                                               ignoreDeg = NULL,
                                               eValType = eigenvalue_decay,
                                               N = no_functions),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
  scores <- t(replicate(no_functions,
                        stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
                                                            type = eigenvalue_decay)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    X <- scores %*% trueFuns[[j]]@X
    # add the mean here:
    if(!is.null(mean_fns)){
      #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
      X = X + rep(mean_fns[[j]](argvals[[j]]), each = no_functions)
    }
    if(length(specialp) == 0){
      if(n_outliers > 0){
        qcoeff <- rbinom(n_outliers, 1, kprob) 
        qcoeff[qcoeff == 0] <- -1 
        X[1:n_outliers, ] <- X[1:n_outliers, ] + 8*qcoeff 
      }
    }
    else{
      if(!(j %in% specialp)){
        if(n_outliers > 0){
          qcoeff <- rbinom(n_outliers, 1, kprob) 
          qcoeff[qcoeff == 0] <- -1 
          X[1:n_outliers, ] <- X[1:n_outliers, ] + 8*qcoeff 
        }
      }
    }
    
    
    if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  simData = multiFunData(simData)
  simData <- addErrorMult(simData,
                               sd = runif(length(argvals),
                                          error_params[1],
                                          error_params[2]))
  
  return(simData)
}



# multimodel2_1 <- function(no_functions = 100,
#                           argvals = list(seq(0,1, length.out = 50),
#                                          seq(0,1, length.out = 50),
#                                          seq(0, 1, length.out = 50)),
#                           basisfunction = "Fourier",
#                           eigenfntype = "split",
#                           nbasis = 9,
#                           eigenvalue_decay = "linear",
#                           mean_fns = list(mu1 = function(t){return(4*t)},
#                                           mu2 = function(t){return(30*t* ((1-t)^(3/2)) )},
#                                           mu3 = function(t){return(5 *((t-1)^2))}),
#                           error_params = c(.1, .3),
#                           kprob = 0.5,
#                           n_outliers = 20,
#                           specialp = c(1))
# 
# 
# {
#   
#   # generate eigenfunctions
#   trueFuns <- switch(eigenfntype,
#                      split = simMultiSplit(argvals = argvals,
#                                            M = nbasis,
#                                            eFunType = basisfunction,
#                                            ignoreDeg = NULL,
#                                            eValType = eigenvalue_decay,
#                                            N = no_functions),
#                      weighted = simMultiWeight(argvals = argvals,
#                                                M = nbasis,
#                                                eFunType = basisfunction,
#                                                ignoreDeg = NULL,
#                                                eValType = eigenvalue_decay,
#                                                N = no_functions),
#                      stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
#   
#   # number of eigenfunctions generated
#   Mtotal <- nObs(trueFuns)
#   
#   # number of elements in multivariate functional basis
#   p <- length(trueFuns)
#   
#   # generate eigenvalues and scores
#   trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
#   scores <- t(replicate(no_functions,
#                         stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
#                                                             type = eigenvalue_decay)))))
#   
#   # generate individual observations
#   simData  <- vector("list", p)
#   
#   for(j in seq_len(p))
#   {
#     #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
#     X <- scores %*% trueFuns[[j]]@X
#     # add the mean here:
#     if(!is.null(mean_fns)){
#       #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
#       X = X + rep(mean_fns[[j]](argvals[[j]]), each = no_functions)
#     }
#     if(!(j %in% specialp)){
#       if(n_outliers > 0){
#         qcoeff <- rbinom(n_outliers, 1, kprob) 
#         qcoeff[qcoeff == 0] <- -1 
#         X[1:n_outliers, ] <- X[1:n_outliers, ] + 8*qcoeff 
#       }
#     }
#     if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
#     simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
#   } 
#   simData = multiFunData(simData)
#   simData <- addErrorMult(simData,
#                           sd = runif(length(argvals),
#                                      error_params[1],
#                                      error_params[2]))
#   
#   return(simData)
# }

multimodel3 <- function(no_functions = 100,
                          argvals = list(seq(0,1, length.out = 50),
                                         seq(0,1, length.out = 50),
                                         seq(0, 1, length.out = 50)),
                          basisfunction = "Fourier",
                          eigenfntype = "split",
                          nbasis = 9,
                          eigenvalue_decay = "linear",
                          mean_fns = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                          mu2 = function(t){return(5 * cos(2*pi*t))},
                                          mu3 = function(t){return(5 * ((t-1)^2) )}),
                          error_params = c(.1, .3),
                          kprob = 0.5,
                          n_outliers = 20, 
                          a = 0, b = 0.9, l = 0.1, 
                          specialp = c())


{
  
  # generate eigenfunctions
  trueFuns <- switch(eigenfntype,
                     split = simMultiSplit(argvals = argvals,
                                           M = nbasis,
                                           eFunType = basisfunction,
                                           ignoreDeg = NULL,
                                           eValType = eigenvalue_decay,
                                           N = no_functions),
                     weighted = simMultiWeight(argvals = argvals,
                                               M = nbasis,
                                               eFunType = basisfunction,
                                               ignoreDeg = NULL,
                                               eValType = eigenvalue_decay,
                                               N = no_functions),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
  scores <- t(replicate(no_functions,
                        stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
                                                            type = eigenvalue_decay)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    X <- scores %*% trueFuns[[j]]@X
    # add the mean here:
    if(!is.null(mean_fns)){
      #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
      X = X + rep(mean_fns[[j]](argvals[[j]]), each = no_functions)
    }
    if(length(specialp) == 0){
      if(n_outliers > 0){
        qcoeff <- rbinom(n_outliers, 1, kprob) 
        qcoeff[qcoeff == 0] <- -1 
        indicator <- sapply(runif(n_outliers, a, b),
                            function(x) (argvals[[j]] >= x)*(argvals[[j]] <= x + l) )
        X[1:n_outliers, ] <- X[1:n_outliers, ] + t(indicator*rep(qcoeff*8,
                                                                 each = length(argvals[[j]])))
      }
    }else{
      if(!(j %in% specialp)){
        if(n_outliers > 0){
          qcoeff <- rbinom(n_outliers, 1, kprob) 
          qcoeff[qcoeff == 0] <- -1 
          indicator <- sapply(runif(n_outliers, a, b),
                              function(x) (argvals[[j]] >= x)*(argvals[[j]] <= x + l) )
          X[1:n_outliers, ] <- X[1:n_outliers, ] + t(indicator*rep(qcoeff*8,
                                                                   each = length(argvals[[j]])))
        }
      }
    }
    
    if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  simData = multiFunData(simData)
  simData <- addErrorMult(simData,
                          sd = runif(length(argvals),
                                     error_params[1],
                                     error_params[2]))
  
  return(simData)
}
multimodel4 <- function(no_functions = 100,
                        argvals = list(seq(0,1, length.out = 50),
                                       seq(0,1, length.out = 50),
                                       seq(0, 1, length.out = 50)),
                        basisfunction = "Fourier",
                        eigenfntype = "split",
                        nbasis = 9,
                        eigenvalue_decay = "linear",
                        mean_fns = list(mu1 = function(t){return(4*t)},
                                        mu2 = function(t){return(30*t* ((1-t)^(3/2)) )},
                                        mu3 = function(t){return(5 *((t-1)^2))}),
                        error_params = c(.1, .3),
                        kprob = 0.5,
                        n_outliers = 20, 
                        a = 0, b = 0.9, l = 0.1, 
                        specialp = c())


{
  
  # generate eigenfunctions
  trueFuns <- switch(eigenfntype,
                     split = simMultiSplit(argvals = argvals,
                                           M = nbasis,
                                           eFunType = basisfunction,
                                           ignoreDeg = NULL,
                                           eValType = eigenvalue_decay,
                                           N = no_functions),
                     weighted = simMultiWeight(argvals = argvals,
                                               M = nbasis,
                                               eFunType = basisfunction,
                                               ignoreDeg = NULL,
                                               eValType = eigenvalue_decay,
                                               N = no_functions),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
  scores <- t(replicate(no_functions,
                        stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
                                                            type = eigenvalue_decay)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    X <- scores %*% trueFuns[[j]]@X
    # add the mean here:
    if(!is.null(mean_fns)){
      #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
      X = X + rep(mean_fns[[j]](argvals[[j]]), each = no_functions)
    }
    if(length(specialp) == 0){
      if(n_outliers > 0){
        qcoeff <- rbinom(n_outliers, 1, kprob) 
        qcoeff[qcoeff == 0] <- -1 
        indicator <- sapply(runif(n_outliers, a, b),
                            function(x) (argvals[[j]] >= x)*(argvals[[j]] <= x + l) )
        X[1:n_outliers, ] <- X[1:n_outliers, ] + t(indicator*rep(qcoeff*8,
                                                                 each = length(argvals[[j]])))
      }
    }else{
      if(!(j %in% specialp)){
        if(n_outliers > 0){
          qcoeff <- rbinom(n_outliers, 1, kprob) 
          qcoeff[qcoeff == 0] <- -1 
          indicator <- sapply(runif(n_outliers, a, b),
                              function(x) (argvals[[j]] >= x)*(argvals[[j]] <= x + l) )
          X[1:n_outliers, ] <- X[1:n_outliers, ] + t(indicator*rep(qcoeff*8,
                                                                   each = length(argvals[[j]])))
        }
      }
    }
    
    if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  simData = multiFunData(simData)
  simData <- addErrorMult(simData,
                          sd = runif(length(argvals),
                                     error_params[1],
                                     error_params[2]))
  
  return(simData)
}


multimodel5 <- function(no_functions = 100,
                        argvals = list(seq(0,1, length.out = 50),
                                       seq(0,1, length.out = 50),
                                       seq(0, 1, length.out = 50)),
                        basisfunction = "Fourier",
                        eigenfntype = "split",
                        nbasis = 9,
                        eigenvalue_decay = "linear",
                        mean_fns = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                        mu2 = function(t){return(5 * cos(2*pi*t))},
                                        mu3 = function(t){return(5 * ((t-1)^2) )}),
                        mean_fns2 = list(mu1 = function(t){return(5 * sin(2*pi*(t-0.3)))},
                                         mu2 = function(t){return(5 * cos(2*pi*(t-0.2)))},
                                         mu3 = function(t){return(5 * ((.1 - t)^2))}),
                        error_params = c(.1, .3),
                        n_outliers = 10,
                        specialp = c())


{
  
  # generate eigenfunctions
  trueFuns <- switch(eigenfntype,
                     split = simMultiSplit(argvals = argvals,
                                           M = nbasis,
                                           eFunType = basisfunction,
                                           ignoreDeg = NULL,
                                           eValType = eigenvalue_decay,
                                           N = no_functions),
                     weighted = simMultiWeight(argvals = argvals,
                                               M = nbasis,
                                               eFunType = basisfunction,
                                               ignoreDeg = NULL,
                                               eValType = eigenvalue_decay,
                                               N = no_functions),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
  scores <- t(replicate(no_functions,
                        stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
                                                            type = eigenvalue_decay)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    X <- scores %*% trueFuns[[j]]@X
    # add the mean here:
    if(!is.null(mean_fns)){
      #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
      # X[(n_outliers+1):N, ] <- X[(n_outliers+1):N, ] + rep(mufs[[j]](argvals[[j]]),
      #                                                      rep(N-n_outliers, ncol(X)))
      X[(n_outliers+1):no_functions, ] = X[(n_outliers+1):no_functions, ] +  
        rep(mean_fns[[j]](argvals[[j]]),
            each = no_functions - n_outliers)
    }
    if(length(specialp) == 0){
      if(n_outliers > 0){
        X[1:n_outliers, ] <- X[1:n_outliers, ] + rep(mean_fns2[[j]](argvals[[j]]),
                                                     each = n_outliers)
      }
    }else{
      if(!(j %in% specialp)){
        if(n_outliers > 0){
          X[1:n_outliers, ] <- X[1:n_outliers, ] + rep(mean_fns2[[j]](argvals[[j]]),
                                                       each = n_outliers)
        }
      }else{
        X[1:n_outliers, ] <- X[1:n_outliers, ] + rep(mean_fns[[j]](argvals[[j]]),
                                                     each = n_outliers)
      }
    }
    
    if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  simData = multiFunData(simData)
  simData <- addErrorMult(simData,
                          sd = runif(length(argvals),
                                     error_params[1],
                                     error_params[2]))
  
  return(simData)
}


multimodel6 <- function(no_functions = 100,
                        argvals = list(seq(0,1, length.out = 50),
                                       seq(0,1, length.out = 50),
                                       seq(0, 1, length.out = 50)),
                        basisfunction = "Fourier",
                        eigenfntype = "split",
                        nbasis = 9,
                        eigenvalue_decay = "linear",
                        mean_fns = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                        mu2 = function(t){return(5 * cos(2*pi*t))},
                                        mu3 = function(t){return(5 * ((t-1)^2) )}),
                        mean_fns2 = list(u1 = function(t){return(2 * sin(4*pi*t))},
                                        u2 = function(t){return(2 * cos(4*pi*t))},
                                        u3 = function(t){return(2 * cos(8*pi*t))}),
                        error_params = c(.1, .3),
                        n_outliers = 10,
                        specialp = c())


{
  
  # generate eigenfunctions
  trueFuns <- switch(eigenfntype,
                     split = simMultiSplit(argvals = argvals,
                                           M = nbasis,
                                           eFunType = basisfunction,
                                           ignoreDeg = NULL,
                                           eValType = eigenvalue_decay,
                                           N = no_functions),
                     weighted = simMultiWeight(argvals = argvals,
                                               M = nbasis,
                                               eFunType = basisfunction,
                                               ignoreDeg = NULL,
                                               eValType = eigenvalue_decay,
                                               N = no_functions),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
  scores <- t(replicate(no_functions,
                        stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
                                                            type = eigenvalue_decay)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    X <- scores %*% trueFuns[[j]]@X
    # add the mean here:
    if(!is.null(mean_fns)){
      #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
      # X[(n_outliers+1):N, ] <- X[(n_outliers+1):N, ] + rep(mufs[[j]](argvals[[j]]),
      #                                                      rep(N-n_outliers, ncol(X)))
      runc <- runif(1,-2.1,2.1)
      X[(n_outliers+1):no_functions, ] = X[(n_outliers+1):no_functions, ] +  
        rep(mean_fns[[j]](argvals[[j]]), each = no_functions - n_outliers) + 
        runc
    }
    if(length(specialp) == 0){
      if(n_outliers > 0){
        X[1:n_outliers, ] <- X[1:n_outliers, ] + 
          rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) + 
          rep(mean_fns2[[j]](argvals[[j]]), each = n_outliers)
      }
    }else{
      if(!(j %in% specialp)){
        if(n_outliers > 0){
          X[1:n_outliers, ] <- X[1:n_outliers, ] + 
            rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) + 
            rep(mean_fns2[[j]](argvals[[j]]), each = n_outliers)
        }
      }else{
        X[1:n_outliers, ] <- X[1:n_outliers, ] + 
          rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) + 
          runc
      }
    }
    
    if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  simData = multiFunData(simData)
  simData <- addErrorMult(simData,
                          sd = runif(length(argvals),
                                     error_params[1],
                                     error_params[2]))
  
  return(simData)
}

multimodel7 <- function(no_functions = 100,
                        argvals = list(seq(0,1, length.out = 50),
                                       seq(0,1, length.out = 50),
                                       seq(0, 1, length.out = 50)),
                        basisfunction = "Fourier",
                        eigenfntype = "split",
                        nbasis = 9,
                        eigenvalue_decay = "linear",
                        mean_fns = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                        mu2 = function(t){return(4*t )},
                                        mu3 = function(t){return(5 * ((t-1)^2) )}),
                        mean_fns2 = list(u1 = function(t){return(2 * sin(8*pi*t))},
                                         u2 = function(tt, nnn){
                                           runc <- runif(nnn, .25, .75)
                                           sapply(tt, function(x){
                                             2 * sin(4*(x + runc)*pi)
                                           })},
                                         u3 = function(t){return(2 * cos(8*pi*t))}),
                        error_params = c(.1, .3),
                        n_outliers = 10,
                        specialp = c())


{
  
  # generate eigenfunctions
  trueFuns <- switch(eigenfntype,
                     split = simMultiSplit(argvals = argvals,
                                           M = nbasis,
                                           eFunType = basisfunction,
                                           ignoreDeg = NULL,
                                           eValType = eigenvalue_decay,
                                           N = no_functions),
                     weighted = simMultiWeight(argvals = argvals,
                                               M = nbasis,
                                               eFunType = basisfunction,
                                               ignoreDeg = NULL,
                                               eValType = eigenvalue_decay,
                                               N = no_functions),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
  scores <- t(replicate(no_functions,
                        stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
                                                            type = eigenvalue_decay)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    X <- scores %*% trueFuns[[j]]@X
    # add the mean here:
    if(!is.null(mean_fns)){
      #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
      # X[(n_outliers+1):N, ] <- X[(n_outliers+1):N, ] + rep(mufs[[j]](argvals[[j]]),
      #                                                      rep(N-n_outliers, ncol(X)))
      runc <- runif(1,-2.1,2.1)
      X[(n_outliers+1):no_functions, ] = X[(n_outliers+1):no_functions, ] +  
        rep(mean_fns[[j]](argvals[[j]]), each = no_functions - n_outliers) + 
        runc
    }
    if(length(specialp) == 0){
      if(n_outliers > 0){
        if(j == 2){
          X[1:n_outliers, ] <- X[1:n_outliers, ] + 
            rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) + 
            mean_fns2[[j]](argvals[[j]], n_outliers)
        }else{
          X[1:n_outliers, ] <- X[1:n_outliers, ] + 
            rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) + 
            rep(mean_fns2[[j]](argvals[[j]]), each = n_outliers)
        }
        
      }
    }else{
      if(!(j %in% specialp)){
        if(n_outliers > 0){
          if(j == 2){
            X[1:n_outliers, ] <- X[1:n_outliers, ] + 
              rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) + 
              mean_fns2[[j]](argvals[[j]], n_outliers)
          }else{
            X[1:n_outliers, ] <- X[1:n_outliers, ] + 
              rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) + 
              rep(mean_fns2[[j]](argvals[[j]]), each = n_outliers)
          }
        }
      }else{
        X[1:n_outliers, ] <- X[1:n_outliers, ] + 
          rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) + 
          runc
      }
    }
    
    if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  simData = multiFunData(simData)
  simData <- addErrorMult(simData,
                          sd = runif(length(argvals),
                                     error_params[1],
                                     error_params[2]))
  
  return(simData)
}



multimodel8 <- function(no_functions = 100,
                        argvals = list(seq(0,1, length.out = 50),
                                       seq(0,1, length.out = 50),
                                       seq(0, 1, length.out = 50)),
                        basisfunction = "Fourier",
                        eigenfntype = "split",
                        nbasis = 9,
                        eigenvalue_decay = "linear",
                        mean_fns = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                        mu2 = function(t){return(5 * cos(2*pi*t))},
                                        mu3 = function(t){return(5 * ((t-1)^2) )}),
                        mean_fns2 = list(u1 = function(t, n){return(rep(5 * sin(2*pi*t), each = n) *
                                                                      (2+rexp(n, 2)))},
                                         u2 = function(t, n){return(rep(5 * cos(2*pi*t), each = n) * 
                                                                      (2+rexp(n, 2)))},
                                         u3 = function(t, n ){return(rep((5 * (t-1)^2), each = n) * 
                                                                       (2+rexp(n, 2)) -6)}),
                        error_params = c(.1, .3),
                        n_outliers = 10,
                        specialp = c())


{
  
  # generate eigenfunctions
  trueFuns <- switch(eigenfntype,
                     split = simMultiSplit(argvals = argvals,
                                           M = nbasis,
                                           eFunType = basisfunction,
                                           ignoreDeg = NULL,
                                           eValType = eigenvalue_decay,
                                           N = no_functions),
                     weighted = simMultiWeight(argvals = argvals,
                                               M = nbasis,
                                               eFunType = basisfunction,
                                               ignoreDeg = NULL,
                                               eValType = eigenvalue_decay,
                                               N = no_functions),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
  scores <- t(replicate(no_functions,
                        stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
                                                            type = eigenvalue_decay)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    X <- scores %*% trueFuns[[j]]@X
    # add the mean here:
    if(!is.null(mean_fns)){
      #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
      # X[(n_outliers+1):N, ] <- X[(n_outliers+1):N, ] + rep(mufs[[j]](argvals[[j]]),
      #                                                      rep(N-n_outliers, ncol(X)))
      X[(n_outliers+1):no_functions, ] = X[(n_outliers+1):no_functions, ] +  
        rep(mean_fns[[j]](argvals[[j]]),
            each = no_functions - n_outliers)
    }
    if(length(specialp) == 0){
      if(n_outliers > 0){
        # X[1:n_outliers, ] <- X[1:n_outliers, ] + rep(mean_fns2[[j]](argvals[[j]]),
        #                                              each = n_outliers)
        X[1:n_outliers, ] <- X[1:n_outliers, ] + 
          rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) +
          mean_fns2[[j]](argvals[[j]], n_outliers) 
      }
    }else{
      if(!(j %in% specialp)){
        if(n_outliers > 0){
          X[1:n_outliers, ] <- X[1:n_outliers, ] + 
            rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) +
            mean_fns2[[j]](argvals[[j]], n_outliers) 
        }
      }else{
        X[1:n_outliers, ] <- X[1:n_outliers, ] + rep(mean_fns[[j]](argvals[[j]]),
                                                     each = n_outliers)
      }
    }
    
    if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  simData = multiFunData(simData)
  simData <- addErrorMult(simData,
                          sd = runif(length(argvals),
                                     error_params[1],
                                     error_params[2]))
  
  return(simData)
}


multimodel9 <- function(no_functions = 100,
                        argvals = list(seq(0,1, length.out = 50),
                                       seq(0,1, length.out = 50),
                                       seq(0, 1, length.out = 50)),
                        basisfunction = "Fourier",
                        eigenfntype = "split",
                        nbasis = 9,
                        eigenvalue_decay = "linear",
                        mean_fns = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                        mu2 = function(t){return(5 * cos(2*pi*t))},
                                        mu3 = function(t){return(5 * ((t-1)^2) )}),
                        mean_fns2 = list(u1 = function(t){return(runif(1,10,10)* t * sin(pi*t))},
                                         u2 = function(t){return(runif(1,11,11)* t * cos(pi*t))},
                                         u3 = function(t){return((runif(1,10,10)*t * sin(2*pi*t))-6)}),
                        mean_fns3 = list(u1 = function(t, z4){return(8*t * sin(pi*t))},
                                         u2 = function(t, z4){return((8-7)*t * cos(pi*t))},
                                         u3 = function(t, z4 ){return(((8-2) * sin(2*pi*t))-3)}),
                        error_params = c(.1, .3),
                        n_outliers = 10,
                        specialp = c())


{
  
  # generate eigenfunctions
  trueFuns <- switch(eigenfntype,
                     split = simMultiSplit(argvals = argvals,
                                           M = nbasis,
                                           eFunType = basisfunction,
                                           ignoreDeg = NULL,
                                           eValType = eigenvalue_decay,
                                           N = no_functions),
                     weighted = simMultiWeight(argvals = argvals,
                                               M = nbasis,
                                               eFunType = basisfunction,
                                               ignoreDeg = NULL,
                                               eValType = eigenvalue_decay,
                                               N = no_functions),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data."))
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(M = Mtotal, type = eigenvalue_decay)
  scores <- t(replicate(no_functions,
                        stats::rnorm(Mtotal, sd = sqrt(eVal(M = Mtotal,
                                                            type = eigenvalue_decay)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  Z4 <- runif(1, 2, 8)
  for(j in seq_len(p))
  {
    #X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    X <- scores %*% trueFuns[[j]]@X
    # add the mean here:
    if(!is.null(mean_fns)){
      #X = X + rep(mean_fns[[j]](argvals[[j]]), rep(N, ncol(X)))
      # X[(n_outliers+1):N, ] <- X[(n_outliers+1):N, ] + rep(mufs[[j]](argvals[[j]]),
      #                                                      rep(N-n_outliers, ncol(X)))
      X[(n_outliers+1):no_functions, ] = X[(n_outliers+1):no_functions, ] +  
        rep(mean_fns[[j]](argvals[[j]]),
            each = no_functions - n_outliers) +
        rep(mean_fns3[[j]](argvals[[j]], Z4),
            each = no_functions - n_outliers) 
    }
    if(length(specialp) == 0){
      if(n_outliers > 0){
        # X[1:n_outliers, ] <- X[1:n_outliers, ] + rep(mean_fns2[[j]](argvals[[j]]),
        #                                              each = n_outliers)
        X[1:n_outliers, ] <- X[1:n_outliers, ] + 
          rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) +
          rep(mean_fns2[[j]](argvals[[j]]), each = n_outliers) 
      }
    }else{
      if(!(j %in% specialp)){
        if(n_outliers > 0){
          X[1:n_outliers, ] <- X[1:n_outliers, ] + 
            rep(mean_fns[[j]](argvals[[j]]), each = n_outliers) +
            rep(mean_fns2[[j]](argvals[[j]]), each = n_outliers) 
        }
      }else{
        X[1:n_outliers, ] <- X[1:n_outliers, ] + 
          rep(mean_fns[[j]](argvals[[j]]),
              each = n_outliers) + 
          rep(mean_fns3[[j]](argvals[[j]], Z4),
              each = n_outliers) 
      }
    }
    
    if(no_functions == 1) dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  simData = multiFunData(simData)
  simData <- addErrorMult(simData,
                          sd = runif(length(argvals),
                                     error_params[1],
                                     error_params[2]))
  
  return(simData)
}
# simMultiFunData2 <- function(type = "split",
#                              argvals = list(seq(0,1, length.out = 50),
#                                             seq(0,1, length.out = 50),
#                                             seq(0, 1, length.out = 50)),
#                              M = 9, eFunType = "Fourier", ignoreDeg = NULL,
#                              eValType = "linear", N = 100,
#                              mufs = list(mu1 = function(t){return(5 * sin(2*pi*t))},
#                                          mu2 = function(t){return(5 * cos(2*pi*t))},
#                                          mu3 = function(t){return(5 * ((t-1)^2) )}), 
#                              kprob = 0.5, n_outliers = 20, 
#                              a = 0, b = 0.9, l = 0.1)
# {
#   
#   # generate eigenfunctions
#   trueFuns <- switch(type,
#                      split = simMultiSplit(argvals, M, eFunType, ignoreDeg, eValType, N),
#                      weighted = simMultiWeight(argvals, M, eFunType, ignoreDeg, eValType, N),
#                      stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data.")
#   )
#   
#   # number of eigenfunctions generated
#   Mtotal <- nObs(trueFuns)
#   
#   # number of elements in multivariate functional basis
#   p <- length(trueFuns)
#   
#   # generate eigenvalues and scores
#   trueVals <- eVal(Mtotal, eValType)
#   scores <- t(replicate(N, stats::rnorm(Mtotal, sd = sqrt(eVal(Mtotal, eValType)))))
#   
#   # generate individual observations
#   simData  <- vector("list", p)
#   
#   for(j in seq_len(p))
#   {
#     X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
#     # add the mean here:
#     if(!is.null(mufs)){
#       X = X + rep(mufs[[j]](argvals[[j]]), rep(N, ncol(X)))
#     }
#     
#     if(n_outliers > 0){
#       qcoeff <- rbinom(n_outliers, 1, kprob) 
#       qcoeff[qcoeff == 0] <- -1 
#       indicator <- sapply(runif(n_outliers, a, b),
#                           function(x) (argvals[[j]] >= x)*(argvals[[j]] <= x + l) )
#       
#       X[1:n_outliers, ] <- X[1:n_outliers, ] + t(indicator*rep(qcoeff*8,
#                                                                each = length(argvals[[j]])))
#     }
#     if(N == 1)
#       dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
#     
#     simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
#   } 
#   
#   return(list(simData = multiFunData(simData),
#               trueFuns = trueFuns,
#               trueVals = trueVals))
# }

# simMultiFunData1 <- function(type = "split",
#                              argvals = list(seq(0,1, length.out = 50),
#                                             seq(0,1, length.out = 50),
#                                             seq(0, 1, length.out = 50)),
#                              M = 9, eFunType = "Fourier", ignoreDeg = NULL,
#                              eValType = "linear", N = 100,
#                              mufs = list(mu1 = function(t){return(5 * sin(2*pi*t))},
#                                          mu2 = function(t){return(5 * cos(2*pi*t))},
#                                          mu3 = function(t){return(5 * ((t-1)^2) )}), 
#                              kprob = 0.5, n_outliers = 20)

# ssm1 <- simMultiFunData1()
# funData::plot(ssm1$simData)
# ssm1$simData <- addErrorMult(ssm1$simData, sd = runif(3, 0.3, 0.5)) 
# funData::plot(ssm1$simData)


# Model 2 isolated magnitude 


# ssm2 <- simMultiFunData2()
# funData::plot(ssm2$simData)
# ssm2$simData <- addErrorMult(ssm2$simData, sd = runif(3, 0.3, 0.5)) 
# funData::plot(ssm2$simData)
# Model 3 shape outleir I
simMultiFunData3 <- function(type = "split",
                             argvals = list(seq(0,1, length.out = 50),
                                            seq(0,1, length.out = 50),
                                            seq(0, 1, length.out = 50)),
                             M = 9, eFunType = "Fourier", ignoreDeg = NULL,
                             eValType = "linear", N = 100,
                             mufs = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                         mu2 = function(t){return(5 * cos(2*pi*t))},
                                         mu3 = function(t){return(5 * ((t-1)^2) )}),
                             mufs2 = list(mu1 = function(t){return(5 * sin(2*pi*(t-0.3)))},
                                          mu2 = function(t){return(5 * cos(2*pi*(t-0.2)))},
                                          mu3 = function(t){return(5 * ((t - 0.5)^2))}),
                             n_outliers = 10)
{
  # generate eigenfunctions
  trueFuns <- switch(type,
                     split = simMultiSplit(argvals, M, eFunType, ignoreDeg, eValType, N),
                     weighted = simMultiWeight(argvals, M, eFunType, ignoreDeg, eValType, N),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data.")
  )
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(Mtotal, eValType)
  scores <- t(replicate(N, stats::rnorm(Mtotal, sd = sqrt(eVal(Mtotal, eValType)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    # add the mean here:
    if(!is.null(mufs)){
      X[(n_outliers+1):N, ] = X[(n_outliers+1):N, ] + rep(mufs[[j]](argvals[[j]]),
                                                          rep(N-n_outliers, ncol(X)))
    }
    
    if(n_outliers > 0){
      X[1:n_outliers, ] <- X[1:n_outliers, ] + rep(mufs2[[j]](argvals[[j]]),
                                                   rep(n_outliers, ncol(X)))
    }
    if(N == 1)
      dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  
  return(list(simData = multiFunData(simData),
              trueFuns = trueFuns,
              trueVals = trueVals))
}
# ssm3 <- simMultiFunData3()
# funData::plot(ssm3$simData)
# ssm3$simData <- addErrorMult(ssm3$simData, sd = runif(3, 0.1, 0.2)) 
# funData::plot(ssm3$simData)
# model 4 shape outlier 2
simMultiFunData4 <- function(type = "split",
                             argvals = list(seq(0,1, length.out = 50),
                                            seq(0,1, length.out = 50),
                                            seq(0, 1, length.out = 50)),
                             M = 9, eFunType = "Fourier", ignoreDeg = NULL,
                             eValType = "linear", N = 100,
                             mufs = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                         mu2 = function(t){return(5 * cos(2*pi*t))},
                                         mu3 = function(t){return(5 * ((t-1)^2) )}),
                             ufs = list(u1 = function(t){return(2 * sin(4*pi*t))},
                                        u2 = function(t){return(2 * cos(4*pi*t))},
                                        u3 = function(t){return(2 * cos(8*pi*t))}),
                             n_outliers = 10)
{
  # generate eigenfunctions
  trueFuns <- switch(type,
                     split = simMultiSplit(argvals, M, eFunType, ignoreDeg, eValType, N),
                     weighted = simMultiWeight(argvals, M, eFunType, ignoreDeg, eValType, N),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data.")
  )
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(Mtotal, eValType)
  scores <- t(replicate(N, stats::rnorm(Mtotal, sd = sqrt(eVal(Mtotal, eValType)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    # add the mean here:
    if(!is.null(mufs)){
      X[(n_outliers+1):N, ] = X[(n_outliers+1):N, ] +
        rep(mufs[[j]](argvals[[j]]), rep(N-n_outliers, ncol(X))) +
        runif(1,-2.1,2.1)
    }
    
    if(n_outliers > 0){
      X[1:n_outliers, ] <- X[1:n_outliers, ] +
        rep(mufs[[j]](argvals[[j]]), rep(n_outliers, ncol(X))) +
        rep(ufs[[j]](argvals[[j]]), each = n_outliers) 
      
    }
    if(N == 1)
      dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  
  return(list(simData = multiFunData(simData),
              trueFuns = trueFuns,
              trueVals = trueVals))
}
# ssm4 <- simMultiFunData4()
# funData::plot(ssm4$simData)
# ssm4$simData <- addErrorMult(ssm4$simData, sd = runif(3, 0.1, 0.2)) 
# funData::plot(ssm4$simData)
# model 5 mixed outleir
simMultiFunData5 <- function(type = "split",
                             argvals = list(seq(0,1, length.out = 50),
                                            seq(0,1, length.out = 50),
                                            seq(0, 1, length.out = 50)),
                             M = 9, eFunType = "Fourier", ignoreDeg = NULL,
                             eValType = "linear", N = 100,
                             mufs = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                         mu2 = function(t){return(5 * cos(2*pi*t))},
                                         mu3 = function(t){return(5 * ((t-1)^2) )}),
                             ufs = list(u1 = function(t, n){return(rep(5 * sin(2*pi*t), each = n) *
                                                                     (2+rexp(n, 2)))},
                                        u2 = function(t, n){return(rep(5 * cos(2*pi*t), each = n) * 
                                                                     (2+rexp(n, 2)))},
                                        u3 = function(t, n ){return(rep((5 * (t-1)^2), each = n) * 
                                                                      (2+rexp(n, 2)) -6)}),
                             n_outliers = 10)
{
  # generate eigenfunctions
  trueFuns <- switch(type,
                     split = simMultiSplit(argvals, M, eFunType, ignoreDeg, eValType, N),
                     weighted = simMultiWeight(argvals, M, eFunType, ignoreDeg, eValType, N),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data.")
  )
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(Mtotal, eValType)
  scores <- t(replicate(N, stats::rnorm(Mtotal, sd = sqrt(eVal(Mtotal, eValType)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    # add the mean here:
    if(!is.null(mufs)){
      X[(n_outliers+1):N, ] = X[(n_outliers+1):N, ] +
        rep(mufs[[j]](argvals[[j]]), rep(N-n_outliers, ncol(X)))
    }
    
    if(n_outliers > 0){
      X[1:n_outliers, ] <- X[1:n_outliers, ] +
        rep(mufs[[j]](argvals[[j]]), rep(n_outliers, ncol(X))) +
        ufs[[j]](argvals[[j]], n_outliers) 
      
    }
    if(N == 1)
      dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  
  return(list(simData = multiFunData(simData),
              trueFuns = trueFuns,
              trueVals = trueVals))
}
# ssm5 <- simMultiFunData5()
# funData::plot(ssm5$simData)
# ssm5$simData <- addErrorMult(ssm5$simData, sd = runif(3, 0.1, 0.2)) 
# funData::plot(ssm5$simData)

# model 6 joint outlier
simMultiFunData6 <- function(type = "split",
                             argvals = list(seq(0,1, length.out = 50),
                                            seq(0,1, length.out = 50),
                                            seq(0, 1, length.out = 50)),
                             M = 9, eFunType = "Fourier", ignoreDeg = NULL,
                             eValType = "linear", N = 100,
                             mufs = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                         mu2 = function(t){return(5 * cos(2*pi*t))},
                                         mu3 = function(t){return(5 * ((t-1)^2) )}),
                             ufs_out = list(u1 = function(t){return(runif(1,2,8)* t * sin(pi*t))},
                                            u2 = function(t, n){return(runif(1,2,8)* t * cos(pi*t))},
                                            u3 = function(t, n ){return((runif(1,2,8)*t * sin(2*pi*t))-6)}),
                             ufs_no = list(u1 = function(t, z4){return(z4*t * sin(pi*t))},
                                           u2 = function(t, z4){return((8-z4)*t * cos(pi*t))},
                                           u3 = function(t, z4 ){return(((z4-2) * sin(2*pi*t))-6)}),
                             n_outliers = 10)
{
  # generate eigenfunctions
  trueFuns <- switch(type,
                     split = simMultiSplit(argvals, M, eFunType, ignoreDeg, eValType, N),
                     weighted = simMultiWeight(argvals, M, eFunType, ignoreDeg, eValType, N),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data.")
  )
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(Mtotal, eValType)
  scores <- t(replicate(N, stats::rnorm(Mtotal, sd = sqrt(eVal(Mtotal, eValType)))))
  
  # generate individual observations
  simData  <- vector("list", p)
  Z4 <- runif(1, 2, 8)
  for(j in seq_len(p))
  {
    X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    # add the mean here:
    if(!is.null(mufs)){
      X[(n_outliers+1):N, ] = X[(n_outliers+1):N, ] +
        rep(mufs[[j]](argvals[[j]]), rep(N-n_outliers, ncol(X))) +
        rep(ufs_no[[j]](argvals[[j]], Z4), rep(N-n_outliers, ncol(X)))
    }
    
    if(n_outliers > 0){
      temp <- ufs_out[[j]](argvals[[j]])
      X[1:n_outliers, ] <- X[1:n_outliers, ] +
        rep(mufs[[j]](argvals[[j]]), rep(n_outliers, ncol(X))) +
        rep(temp, rep(n_outliers, ncol(X)))
      
    }
    if(N == 1)
      dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  
  return(list(simData = multiFunData(simData),
              trueFuns = trueFuns,
              trueVals = trueVals))
}
# ssm6 <- simMultiFunData6()
# funData::plot(ssm6$simData)
# ssm6$simData <- addErrorMult(ssm5$simData, sd = runif(3, 0.1, 0.2)) 
# funData::plot(ssm6$simData)

#model 7 covariance outlier
simMultiFunData7 <- function(type = "split",
                             argvals = list(seq(0,1, length.out = 50),
                                            seq(0,1, length.out = 50),
                                            seq(0, 1, length.out = 50)),
                             M = 9, eFunType = "Fourier", ignoreDeg = NULL,
                             eValType = "linear", N = 100,
                             mufs = list(mu1 = function(t){return(5 * sin(2*pi*t))},
                                         mu2 = function(t){return(5 * cos(2*pi*t))},
                                         mu3 = function(t){return(5 * ((t-1)^2) )}),
                             n_outliers = 10, 
                             nudiag_in = runif(3, 2,3), nudiag_out = runif(3, .3,.5))
{
  # generate eigenfunctions
  trueFuns <- switch(type,
                     split = simMultiSplit(argvals, M, eFunType, ignoreDeg, eValType, N),
                     weighted = simMultiWeight(argvals, M, eFunType, ignoreDeg, eValType, N),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data.")
  )
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eVal(Mtotal, eValType)
  scores <- t(replicate(N, stats::rnorm(Mtotal, sd = sqrt(eVal(Mtotal, eValType)))))
  
  # generate matern errors
  larg <- length(argvals[[1]])
  edata = array(0, dim = c(N, larg, p))
  if(n_outliers > 0){
    for (i in 1:n_outliers)  {
      model <- RandomFields::RMparswm(nudiag = nudiag_out)
      edata[i, , ] = matrix(RandomFields::RFsimulate(model, argvals[[1]]), larg, p)
    }
    for (i in (n_outliers+1):N)  {
      model <- RandomFields::RMparswm(nudiag = nudiag_in)
      edata[i, , ] = matrix(RandomFields::RFsimulate(model, argvals[[1]]), larg, p)
    }
  }else{
    for (i in 1:N)  {
      model <- RandomFields::RMparswm(nudiag = nudiag_in)
      edata[i, , ] = matrix(RandomFields::RFsimulate(model, argvals[[1]]), larg, p)
    }
  }
  
  # generate individual observations
  simData  <- vector("list", p)
  for(j in seq_len(p))
  {
    X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    # add the mean here:
    if(!is.null(mufs)){
      X = X + rep(mufs[[j]](argvals[[j]]), rep(N, ncol(X))) 
    }
    X <- X + edata[,,j]
    if(N == 1)
      dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  return(list(simData = multiFunData(simData),
              trueFuns = trueFuns,
              trueVals = trueVals))
}
# ssm7 <- simMultiFunData7()
# funData::plot(ssm7$simData)
