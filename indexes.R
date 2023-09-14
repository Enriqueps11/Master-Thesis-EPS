###############################################
################## FUNCTIONS ##################
###############################################

## All the functions created for this paper are defined in this R script.

######################
### Epigraph Index ###
######################

## The epigraph index of a curve x (EI) is one minus the proportion of curves in the sample that lie above x.
# This function computes the EI of a bunch of B curves.

EI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.min <- curves[i,]
    for (j in 1:B){
      inpoints <- sum((curves[j,] >= env.min))
      if (inpoints == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (1-index)
}

#######################
### Hypograph Index ###
#######################

## The hypograph index of a curve x (HI) is the proportion of curves in the sample that lie below x.
# This function computes the HI of a bunch of B curves.

HI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.max <- curves[i,]
    for (j in 1:B){
      inpoints <- sum(curves[j,]<=env.max)
      if (inpoints == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (index)
}

##################################
### Generalized Epigraph Index ###
##################################

## The generalized epigraph index of a curve x (MEI) is 
## one minus the "proportion of time" that the curves of the sample are above x.
# This function computes the MEI of a bunch of B curves.

MEI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.min <- curves[i,]
    for (j in 1:B){
      inpoints <- sum(curves[j,] >= env.min)
      index[i] <- index[i]+inpoints
    }
  }
  index <- index/(B*lengthcurves)
  return (1-index)
}

###################################
### Generalized Hypograph Index ###
###################################

## The generalized hypograph index of a curve x (MHI) is 
## the "proportion of time" that the curves of the sample are below x.
# This function computes the MHI of a bunch of B curves.

MHI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.max <- curves[i,]
    for (j in 1:B){
      inpoints <- sum(curves[j,]<=env.max)
      index[i] <- index[i]+inpoints
    }
  }
  index <- index/(B*lengthcurves)
  return (index)
}

######################
### Multivariate Epigraph Index ###
######################

## The epigraph index of a multivariate curve x (mulEI) is one minus the 
# proportion of curves in the sample that lie above x (component by component).
# This function computes the mulEI of a bunch of B curves.

mulEI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  D <- dim(curves)[3]
  for (i in 1:B){
    index[i] <- 0
    inpoints_k <- array(rep(NaN,B*lengthcurves*D),dim=c(B,lengthcurves,D))
    for (k in 1:D){
      env.min <- curves[i,,k]
      inpoints <- c()
      for (j in 1:B){
        inpoints<- rbind(inpoints,curves[j,,k] >= env.min)
      }
      inpoints_k[,,k]<-inpoints
    }
    inpoints_mult<-inpoints_k[,,1];
    for (l in 2 : D){
      inpoints_mult<- inpoints_mult*inpoints_k[,,l]
    }
    for (m in 1:B){
      inpoints_f <- sum(inpoints_mult[m,])
      if (inpoints_f == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (1-index)
}

mulHI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  D <- dim(curves)[3]
  for (i in 1:B){
    index[i] <- 0
    inpoints_k <- array(rep(NaN,B*lengthcurves*D),dim=c(B,lengthcurves,D))
    for (k in 1:D){
      env.min <- curves[i,,k]
      inpoints <- c()
      for (j in 1:B){
        inpoints<- rbind(inpoints,curves[j,,k] <= env.min)
      }
      inpoints_k[,,k]<-inpoints
    }
    inpoints_mult<-inpoints_k[,,1];
    for (l in 2 : D){
      inpoints_mult<- inpoints_mult*inpoints_k[,,l]
    }
    for (m in 1:B){
      inpoints_f <- sum(inpoints_mult[m,])
      if (inpoints_f == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (index)
}

mulMHI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  D <- dim(curves)[3]
  for (i in 1:B){
    index[i] <- 0
    inpoints_k <- array(rep(NaN,B*lengthcurves*D),dim=c(B,lengthcurves,D))
    for (k in 1:D){
      env.min <- curves[i,,k]
      inpoints <- c()
      for (j in 1:B){
        inpoints<- rbind(inpoints,curves[j,,k] <= env.min)
      }
      inpoints_k[,,k]<-inpoints
    }
    inpoints_mult<-inpoints_k[,,1];
    for (l in 2 : D){
      inpoints_mult<- inpoints_mult*inpoints_k[,,l]
    }
    for (m in 1:B){
      inpoints_f <- sum(inpoints_mult[m,])
      index[i] <- index[i]+inpoints_f
    }
  }
  index <- index/(B*lengthcurves)
  return (index)
}

mulMEI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  D <- dim(curves)[3]
  for (i in 1:B){
    index[i] <- 0
    inpoints_k <- array(rep(NaN,B*lengthcurves*D),dim=c(B,lengthcurves,D))
    for (k in 1:D){
      env.min <- curves[i,,k]
      inpoints <- c()
      for (j in 1:B){
        inpoints<- rbind(inpoints,curves[j,,k] >= env.min)
      }
      inpoints_k[,,k]<-inpoints
    }
    inpoints_mult<-inpoints_k[,,1];
    for (l in 2 : D){
      inpoints_mult<- inpoints_mult*inpoints_k[,,l]
    }
    for (m in 1:B){
      inpoints_f <- sum(inpoints_mult[m,])
      index[i] <- index[i]+inpoints_f
    }
  }
  index <- index/(B*lengthcurves)
  return (1-index)
}

##########################
### funspline function ###
##########################

## This function fits a smooth B-spline basis for a bunch of curves X and
## calculates first and second derivatives for each curve.

funspline <- function(X,t,basis){
  
  # INPUT
  # X <- functional data object
  # t <- sequence for creating the basis
  # basis <- create.bspline.basis() object
  ys = smooth.basis(argvals=t, y=t(X), fdParobj=basis) #data obtained when applying a smoothed B-spline basis
  smooth <- t(eval.fd(t,ys$fd,0)) #smoothed data
  deriv <- t(eval.fd(t,ys$fd,1)) #first derivatives
  deriv2 <- t(eval.fd(t,ys$fd,2)) #second derivatives
  
  res <- list("spline"=ys, "smooth"=smooth, "deriv"=deriv, "deriv2"=deriv2) #object containing the data and 1st and 2nd derivatives
  return(res)
}

####################
### ind function ###
####################

## This function generates a dataframe containing EI, HI, MEI and MHI
## for each curve on the original data, first and second derivatives.

ind <- function(X,t,rangeval,nbasis,norder){
  
  # INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the
  #             functional data object can be evaluated. By default, we consider the interval [0,1]
  # nbasis <- number of basis functions 
  # norder <- order of b-splines
  
  basisobj <- create.bspline.basis(rangeval = rangeval, nbasis = nbasis, norder = norder) #b-spline basis (order 3 polynomials when norder = 4)
  
  # dta <- X #original data
  
  if(length(dim(X))==3){
    # Multivariate functional data
    
    N <- dim(X)[1]
    P <- dim(X)[2]
    K <- dim(X)[3]
    
    dta <- array(rep(NaN,N*P),dim=c(N,P,K))
    ddta <- array(rep(NaN,N*P),dim=c(N,P,K))
    d2dta <- array(rep(NaN,N*P),dim=c(N,P,K))
    for(k in 1:K){
      dta_funspline <- funspline(X[,,k], t, basisobj)
      dta[,,k] <- dta_funspline$smooth #smoothed data
      ddta[,,k] <- dta_funspline$deriv #first derivative
      d2dta[,,k] <- dta_funspline$deriv2 #second derivative
    }  
    # Applying the indexes to the original data
    dtaEI <- mulEI(dta) 
    dtaHI <- mulHI(dta)
    dtaMEI <- mulMEI(dta)
    dtaMHI <- mulMHI(dta)
    
    # Applying the indexes to the data first derivatives
    ddtaEI <- mulEI(ddta)
    ddtaHI <- mulHI(ddta)
    ddtaMEI <- mulMEI(ddta)
    ddtaMHI <- mulMHI(ddta)
    
    # Applying the indexes to the data second derivatives
    d2dtaEI <- mulEI(d2dta)
    d2dtaHI <- mulHI(d2dta)
    d2dtaMEI <- mulMEI(d2dta)
    d2dtaMHI <- mulMHI(d2dta)
    
    # New multivariate data set
    ind.data <- data.frame(dtaEI,dtaHI,dtaMEI,dtaMHI,ddtaEI,ddtaHI,ddtaMEI, 
                           ddtaMHI,d2dtaEI,d2dtaHI,d2dtaMEI,d2dtaMHI)
  } else if(length(dim(X))==2){
    # Univariate functional data
    
    dtaX <- funspline(X,t,basisobj)
    dta <- dtaX$smooth # smoothed data coefs
    ddta <- dtaX$deriv #first derivatives
    d2dta <- dtaX$deriv2 #second derivatives
    
    dtaEI <- EI(dta) 
    dtaHI <- HI(dta)
    dtaMEI <- MEI(dta)
    dtaMHI <- MHI(dta)
    
    # Applying the indexes to the data first derivatives
    ddtaEI <- EI(ddta)
    ddtaHI <- HI(ddta)
    ddtaMEI <- MEI(ddta)
    ddtaMHI <- MHI(ddta)
    
    # Applying the indexes to the data second derivatives
    d2dtaEI <- EI(d2dta)
    d2dtaHI <- HI(d2dta)
    d2dtaMEI <- MEI(d2dta)
    d2dtaMHI <- MHI(d2dta)
    
    # New multivariate data set
    ind.data <- data.frame(dtaEI,dtaHI,dtaMEI,dtaMHI,ddtaEI,ddtaHI,ddtaMEI, 
                           ddtaMHI,d2dtaEI,d2dtaHI,d2dtaMEI,d2dtaMHI)
  } else{
    print("Non valid data dimension")
    ind.data <- "Non valid data dimension"
  }
  
  results <- list("ind"=ind.data, "sdata"=dta, "sddata"=ddta, "sd2data"=d2dta)
  
  return(results)
}

##########################################################
### ind function using the medians of the 3 components ###
##########################################################

## This function generates a dataframe containing EI, HI, MEI and MHI
## for each curve on the original data, first and second derivatives and computes de medians 
## the indexes taking into account the 3 components.

indM <- function(X,t,rangeval=c(0,1),nbasis=40,norder=4){
  
  # Applying the indexes calculated in the roahd way
  
  # INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the
  #             functional data object can be evaluated. By default, we consider the interval [0,1]
  # nbasis <- number of basis functions --- PARAMETRO QUE NO ESTOY USANDO AHORA
  # norder <- order of b-splines
  
  basisobj <- create.bspline.basis(rangeval=rangeval,norder=norder) #b-spline basis (order 3 polynomials when norder = 4)
  
  if(length(dim(X))==3){
    # Multivariate functional data
    
    N <- dim(X)[1]
    P <- dim(X)[2]
    K <- dim(X)[3]
    
    dta <- array(rep(NaN,N*P),dim=c(N,P,K))
    ddta <- array(rep(NaN,N*P),dim=c(N,P,K))
    d2dta <- array(rep(NaN,N*P),dim=c(N,P,K))
    for(k in 1:K){
      dta_funspline <- funspline(X[,,k],t,basisobj)
      dta[,,k] <- dta_funspline$smooth
      ddta[,,k] <- dta_funspline$deriv #first derivative
      d2dta[,,k] <- dta_funspline$deriv2 #second derivative
    }
    # Applying the indexes to the origial data
    dtaEI <- rowMeans(apply(dta,3,EI))
    dtaHI <- rowMeans(apply(dta,3,HI))
    dtaMEI <- rowMeans(apply(dta,3,MEI))
    dtaMHI <- rowMeans(apply(dta,3,MHI))
    
    # Applying the indexes to the data first derivatives
    ddtaEI <- rowMeans(apply(ddta,3,EI))
    ddtaHI <- rowMeans(apply(ddta,3,HI))
    ddtaMEI <- rowMeans(apply(ddta,3,MEI))
    ddtaMHI <- rowMeans(apply(ddta,3,MHI))
    
    # Applying the indexes to the data second derivatives
    d2dtaEI <- rowMeans(apply(d2dta,3,EI))
    d2dtaHI <- rowMeans(apply(d2dta,3,HI))
    d2dtaMEI <- rowMeans(apply(d2dta,3,MEI))
    d2dtaMHI <- rowMeans(apply(d2dta,3,MHI))
  }else{
    stop("The given data must be a three dimensional array")
  }
  
  # New multivariate data set
  ind.data <- data.frame(dtaEI,dtaHI,dtaMEI,dtaMHI,ddtaEI,ddtaHI,ddtaMEI,
                         ddtaMHI,d2dtaEI,d2dtaHI,d2dtaMEI,d2dtaMHI)
  
  results <- list("ind"=ind.data, "sdata"=dta, "sddata"=ddta, "sd2data"=d2dta)
  
  return(results)
}