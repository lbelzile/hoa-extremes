# Derivatives of log-likelihood, weighting matrices, canonical parameters and derivatives
dydmut <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  t - t0
}
dydp <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  -sigma*(-log(-pexc + 1))^(-xi - 1)/(pexc - 1)
}
dydsigma <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  -(mut*t - mut*t0 - y[,1] + z)/sigma
}
dydxi <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  sigma*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)*((xi*(sigma*(1/(-log(-pexc + 1))^xi - 1)/xi^2 + sigma*log(-log(-pexc + 1))/(xi*(-log(-pexc + 1))^xi))/sigma + (mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)/sigma)/(((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)*xi) - log(-(mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)/xi^2)
}

dyrdmut <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  -sigma*((t - t0)*(-(mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma - (t - t0)*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma)*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(1/xi + 1)
}
dyrdp <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  sigma*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(1/xi + 1)*((-(mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*(-log(-pexc + 1))^(-xi - 1)/(pexc - 1) - (-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*(-log(-pexc + 1))^(-xi - 1)/(pexc - 1))
}
dyrdsigma <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  sigma*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(1/xi + 1)*((-(mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma^2 + (1/(-log(-pexc + 1))^xi - 1)/sigma)/xi - (-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma^2 + (1/(-log(-pexc + 1))^xi - 1)/sigma)/xi)
}
dyrdxi <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  sigma*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(1/xi + 1)*
    (((xi*(sigma*(1/(-log(-pexc + 1))^xi - 1)/xi^2 + sigma*log(-log(-pexc + 1))/(xi*(-log(-pexc + 1))^xi))/sigma + (mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)/sigma)/
        (((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)*xi) - log(-(mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)/xi^2)/(-(mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) -
       ((xi*(sigma*(1/(-log(-pexc + 1))^xi - 1)/xi^2 + sigma*log(-log(-pexc + 1))/(xi*(-log(-pexc + 1))^xi))/sigma + (mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)/sigma)/(((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)*xi) - log(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)/xi^2)/(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(1/xi))
}

dlldy <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  (xi + 1)/(sigma*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1))

}
dlldyr <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  (-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma + xi*(1/xi + 1)/(sigma*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1))
}

dlldydmut <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  -(t - t0)*xi^2*(1/xi + 1)/(sigma^2*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)^2)
}
dlldydsigma <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  xi*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma^2 + (1/(-log(-pexc + 1))^xi - 1)/sigma)*(1/xi + 1)/(sigma*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)^2) - xi*(1/xi + 1)/(sigma^2*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1))
}
dlldydpex <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  xi^2*(-log(-pexc + 1))^(-xi - 1)*(1/xi + 1)/((pexc - 1)*sigma*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)^2)
}

dlldydxi <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  -xi*(xi*(sigma*(1/(-log(-pexc + 1))^xi - 1)/xi^2 + sigma*log(-log(-pexc + 1))/(xi*(-log(-pexc + 1))^xi))/sigma + (mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)/sigma)*(1/xi + 1)/(sigma*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)^2) + (1/xi + 1)/(sigma*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)) - 1/(sigma*((mut*t - mut*t0 - y[,1] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)*xi)
}

dlldyrdmut <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  -(t - t0)*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*xi*(1/xi - 1)/sigma^2 + 2*(t - t0)*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)/sigma^2 - (t - t0)*xi^2*(1/xi + 1)/(sigma^2*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)^2)
}
dlldyrdsigma <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  (-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma^2 + (1/(-log(-pexc + 1))^xi - 1)/sigma)*(1/xi - 1)/sigma - 2*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma^2 + (1/(-log(-pexc + 1))^xi - 1)/sigma)/(sigma*xi) + xi*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma^2 + (1/(-log(-pexc + 1))^xi - 1)/sigma)*(1/xi + 1)/(sigma*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)^2) - (-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma^2 - xi*(1/xi + 1)/(sigma^2*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1))
}
dlldyrdpex <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  (-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*xi*(-log(-pexc + 1))^(-xi - 1)*(1/xi - 1)/((pexc - 1)*sigma) - 2*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*(-log(-pexc + 1))^(-xi - 1)/((pexc - 1)*sigma) + xi^2*(-log(-pexc + 1))^(-xi - 1)*(1/xi + 1)/((pexc - 1)*sigma*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)^2)
}

dlldyrdxi <- function(y, mut, sigma, pexc, xi, z, t, t0, ...){
  (-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((xi*(sigma*(1/(-log(-pexc + 1))^xi - 1)/xi^2 + sigma*log(-log(-pexc + 1))/(xi*(-log(-pexc + 1))^xi))/sigma + (mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)/sigma)*(1/xi - 1)/((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1) - log(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)/xi^2)/sigma -
    2*(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((xi*(sigma*(1/(-log(-pexc + 1))^xi - 1)/xi^2 + sigma*log(-log(-pexc + 1))/(xi*(-log(-pexc + 1))^xi))/sigma + (mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)/sigma)/(((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)*xi) -
                                                                                                          log(-(mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma + 1)/xi^2)/sigma - xi*(xi*(sigma*(1/(-log(-pexc + 1))^xi - 1)/xi^2 + sigma*log(-log(-pexc + 1))/(xi*(-log(-pexc + 1))^xi))/sigma + (mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)/sigma)*(1/xi + 1)/(sigma*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)^2) +
    (1/xi + 1)/(sigma*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)) - 1/(sigma*((mut*t - mut*t0 - y[,2] + z - sigma*(1/(-log(-pexc + 1))^xi - 1)/xi)*xi/sigma - 1)*xi)
}

pex <- function(sigma, mu, mut, xi, z, t0, ...){1-exp(-(1+xi*(z-mu-mut*t0)/sigma)^(-1/xi))}

rlargp.V <- function(dat, par, t, t0, z, ...){
  r <- ncol(dat)
  mut <- par[1]
  pexc <- par[2]
  sigma <- par[3]
  xi <- par[4]
  Varr <- array(0, dim = c(nrow(dat), 2L, 4L)) #n x 2 x p
  Varr[,1,] <- cbind(dydp(y = dat, mut = mut, sigma = sigma, pexc = pexc, xi = xi, z=z, t = t, t0 = t0, ...),
                     dydmut(y = dat, mut = mut, sigma = sigma, pexc = pexc, xi = xi, z=z, t = t, t0 = t0, ...),
                     dydsigma(y = dat, mut = mut, sigma = sigma, pexc = pexc, xi = xi, z=z, t = t, t0 = t0, ...),
                     dydxi(y = dat, mut = mut, sigma = sigma, pexc = pexc, xi = xi, z=z, t = t, t0 = t0, ...))
  Varr[,2,] <- cbind(dyrdp(y = dat, mut = mut, sigma = sigma, pexc = pexc, xi = xi, z=z, t = t, t0 = t0, ...),
                     dyrdmut(y = dat, mut = mut, sigma = sigma, pexc = pexc, xi = xi, z=z, t = t, t0 = t0, ...),
                     dyrdsigma(y = dat, mut = mut, sigma = sigma, pexc = pexc, xi = xi, z=z, t = t, t0 = t0, ...),
                     dyrdxi(y = dat, mut = mut, sigma = sigma, pexc = pexc, xi = xi, z=z, t = t, t0 = t0, ...))
  return(Varr)
}

rlargp.phi <- function(dat, par, z, V, t, t0, ...){
  mut <- par[1]
  pexc <- par[2]
  sigma <- par[3]
  xi <- par[4]
  colSums(V[,1,] * dlldy(y = dat, mut, sigma, pexc, xi, z, t, t0, ...)) + colSums(V[,2,] * dlldyr(y = dat, mut, sigma, pexc, xi, z, t, t0, ...))
}

rlargp.dphi <- function(dat, par, z, V, t, t0, ...){
  mut <- par[1]
  pexc <- par[2]
  sigma <- par[3]
  xi <- par[4]
  rbind(
    colSums(V[,1,] * dlldydmut(y = dat, mut, sigma, pexc, xi, z, t, t0, ...)) + colSums(V[,2,] * dlldyrdmut(y = dat, mut, sigma, pexc, xi, z, t, t0, ...)),
    colSums(V[,1,] * dlldydpex(y = dat, mut, sigma, pexc, xi, z, t, t0, ...)) + colSums(V[,2,] * dlldyrdpex(y = dat, mut, sigma, pexc, xi, z, t, t0, ...)),
    colSums(V[,1,] * dlldydsigma(y = dat, mut, sigma, pexc, xi, z, t, t0, ...)) + colSums(V[,2,] * dlldyrdsigma(y = dat, mut, sigma, pexc, xi, z, t, t0, ...)),
    colSums(V[,1,] * dlldydxi(y = dat, mut, sigma, pexc, xi, z, t, t0, ...)) + colSums(V[,2,] * dlldyrdxi(y = dat, mut, sigma, pexc, xi, z, t, t0, ...)))
}

rlargp.ll <- function(par, dat, z, t, t0){
  dat <- as.matrix(dat)
  if(!isTRUE(all(dat[1,1] >= dat[1,]))){
    stop("Observations in `dat` must be ordered from largest to smallest in each row.")
  }
  r <- ncol(dat)
  n <- nrow(dat)
  xmax <- max(dat[,1])
  xmin <- min(dat[,r])
  mut <- par[1]
  sigma <- par[3]
  xi <- par[4]
  if(par[2] < 0 || par[2] > 1 || par[3] < 0){
    return(-1e10)
  }
  mu <-  z - qgev(p = 1-par[2], loc = 0, scale = sigma, shape = xi) - mut*t0
  loc <- mu + par[1]*t
  if((xi < 0)&(xmax > max(loc) - sigma/xi) || (xi > 0)&(xmin < min(loc) - sigma/xi)){
    return(-1e10)
  }

  if(isTRUE(all.equal(xi,-1, check.attributes = FALSE))){
    - sum((1+xi*(dat[,r]-loc)/sigma)^(-1/xi)) - n*r*log(sigma)
  } else{
    if(abs(xi) > 1e-7){
      - sum((1+xi*(dat[,r]-loc)/sigma)^(-1/xi)) - n*r*log(sigma) -
        (1/xi+1)*sum(log1p(xi*(dat-loc)/sigma))
    } else{
      -sum(exp((loc-dat[,r])/sigma)) - n*r*log(sigma) - sum((dat-loc)/sigma)
    }
  }
}

rlargp.score <- function(par, dat, z, t, t0){
  # return(numDeriv::grad(function(x){rlargp.ll(par = x, dat = dat, z = z, t = t, t0 = t0)}, x = par))
  mut <- par[1]; pex <- par[2]; sigma <- par[3]; xi <- par[4]; r <- 2
  c(sum(-(t - t0)*xi*(1/xi + 1)/(sigma*((mut*t - mut*t0 - dat[,1] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma - 1)) - (t - t0)*xi*(1/xi + 1)/(sigma*((mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma - 1)) - (t - t0)*(-(mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma + 1)^(1/xi - 1)/(sigma*(-(mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma + 1)^(2/xi))),
    sum((-(mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*(-log(-pex + 1))^(-xi - 1)/(pex - 1) + xi*(-log(-pex + 1))^(-xi - 1)*(1/xi + 1)/((pex - 1)*((mut*t - mut*t0 - dat[,1] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma - 1)) + xi*(-log(-pex + 1))^(-xi - 1)*(1/xi + 1)/((pex - 1)*((mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma - 1))),
    sum((-(mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma^2 + (1/(-log(-pex + 1))^xi - 1)/sigma)/xi + ((mut*t - mut*t0 - dat[,1] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma^2 + (1/(-log(-pex + 1))^xi - 1)/sigma)*(1/xi + 1)/((mut*t - mut*t0 - dat[,1] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma - 1) + ((mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma^2 + (1/(-log(-pex + 1))^xi - 1)/sigma)*(1/xi + 1)/((mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma - 1) - r/sigma),
    sum(-(xi*(sigma*(1/(-log(-pex + 1))^xi - 1)/xi^2 + sigma*log(-log(-pex + 1))/(xi*(-log(-pex + 1))^xi))/sigma + (mut*t - mut*t0 - dat[,1] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)/sigma)*(1/xi + 1)/((mut*t - mut*t0 - dat[,1] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma - 1) -
          (xi*(sigma*(1/(-log(-pex + 1))^xi - 1)/xi^2 + sigma*log(-log(-pex + 1))/(xi*(-log(-pex + 1))^xi))/sigma + (mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)/sigma)*(1/xi + 1)/((mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma - 1) + ((xi*(sigma*(1/(-log(-pex + 1))^xi - 1)/xi^2 + sigma*log(-log(-pex + 1))/(xi*(-log(-pex + 1))^xi))/sigma + (mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)/sigma)/(((mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma - 1)*xi) - log(-(mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma + 1)/xi^2)/(-(mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) +
          log(-(mut*t - mut*t0 - dat[,1] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma + 1)/xi^2 + log(-(mut*t - mut*t0 - dat[,2] + z - sigma*(1/(-log(-pex + 1))^xi - 1)/xi)*xi/sigma + 1)/xi^2))

}

rlargp.infomat <- function(par, dat, z, t, t0){
  -numDeriv::hessian(function(x){
    rlargp.ll(par = x, dat = dat, z = z, t = t, t0 = t0)}, x = par)
}

# #Check log-likelihood
# ismev::rlarg.fit(show = FALSE, r=2, xdat = dat, ydat = data.frame(yrs), mul = 1)$nllh +
#   mev::rlarg.ll(par = c(0,1,mlev[4]), (dat-mlev[1]-mlev[2]*yrs)/mlev[3])-prod(dim(dat))*log(mlev[3])
# mev::rlarg.ll(par = c(0,1,mlev[4]), (dat-mlev[1]-mlev[2]*yrs)/mlev[3])-prod(dim(dat))*log(mlev[3])
# rlargp.ll(par = pars, dat = dat, z = z, t = yrs, t0 = t0)
