# Nonhomogeneous Poisson process likelihood and HOA approx
# Parametrized in terms of pth quantile of the T-year maximum

#' Matrix of sufficient directions for the nonhomogeneous
#' Poisson process likelihood parametrized in terms of the
#' p quantile of the T observation maximum
#' @param par vector of parameters estimates, with quantiles of the T\code{ny}-observation, scale and shape parameters
#' @param dat vector of threshold exceedances
#' @param ql probability level of the quantile
#' @param N number of observations for the maximum
nhppq.Vfun <- function(par, dat, ql, N, ny){
  y <- dat
  stopifnot(length(par) == 3L, par[1:2] > 0)
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  V_qp <- ny*(-(qp - y - sigma*((-N/log(ql))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*xi*(1/xi + 1)/sigma^2
  V_sigma <- -ny*(-(qp - y - sigma*((-N/log(ql))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*((qp - y - sigma*((-N/log(ql))^xi - 1)/xi)*xi/sigma^2 + ((-N/log(ql))^xi - 1)/sigma)*(1/xi + 1)/sigma - ny*(-(qp - y - sigma*((-N/log(ql))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma^2
  V_xi <- ny*(-(qp - y - sigma*((-N/log(ql))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*(((sigma*(-N/log(ql))^xi*log(-N/log(ql))/xi - sigma*((-N/log(ql))^xi - 1)/xi^2)*xi/sigma - (qp - y - sigma*((-N/log(ql))^xi - 1)/xi)/sigma)*(1/xi + 1)/((qp - y - sigma*((-N/log(ql))^xi - 1)/xi)*xi/sigma - 1) + log(-(qp - y - sigma*((-N/log(ql))^xi - 1)/xi)*xi/sigma + 1)/xi^2)/sigma
  cbind(V_qp, V_sigma, V_xi)
}

#' Directional derivative of log likelihood
#'
#' This term is \code{V} times the sample
#' space derivative of the log-likelihood
#' @inheritParams nhppq.Vfun
#' @param V an \code{n} by \code{p} matrix whose columns span the sufficient directions
nhppq.dphi <- function(par, dat, ql, N, ny, V){
  y <- dat
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  dphi_qp <- -((xi + 1)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi)*log(-ny/((qp*xi - xi*y - sigma*(-N/log(ql))^xi)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi))) - (xi + 1)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi))/ny
  dphi_sigma <- ((sigma*(-N/log(ql))^xi + qp - y)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi)*log(-ny/((qp*xi - xi*y - sigma*(-N/log(ql))^xi)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi))) - (sigma*(-N/log(ql))^xi + qp - y)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi))/(ny*sigma)
  dphi_xi <- -((qp*xi - xi*y - sigma*(-N/log(ql))^xi)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi)*log(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma) + ((sigma*(-N/log(ql))^xi*log(-N/log(ql)) - qp)*xi^2 + (sigma*(-N/log(ql))^xi*log(-N/log(ql)) - qp)*xi + (xi^2 + xi)*y)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi) - ((qp*xi - xi*y - sigma*(-N/log(ql))^xi)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi)*log(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma) + ((sigma*(-N/log(ql))^xi*log(-N/log(ql)) - qp)*xi^2 + (sigma*(-N/log(ql))^xi*log(-N/log(ql)) - qp)*xi + (xi^2 + xi)*y)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi))*log(-ny/((qp*xi - xi*y - sigma*(-N/log(ql))^xi)*(-(qp*xi - xi*y - sigma*(-N/log(ql))^xi)/sigma)^(1/xi))))/(ny*xi^2)
  rbind(dphi_qp, dphi_sigma, dphi_xi) %*% V
}

#' Directional derivative of the log likelihood
#'
#' This term is \code{V} times the sample
#' space derivative of the log-likelihood
#'
#' @inheritParams nhppq.Vfun
#' @inheritParams nhppq.dphi
nhppq.phi <- function(par, dat, ql, N, ny, V){
  y <- dat
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  phi_a <- sigma*log(ny*(-(qp - y - sigma*((-N/log(ql))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma)/(ny*(-(qp - y - sigma*((-N/log(ql))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1))
  rbind(phi_a) %*% V
}

#' Matrix of sufficient directions for the (truncated)
#' generalized Pareto components likelihood parametrized in terms of the
#' p quantile of the T observation maximum
#' @param par vector of parameters estimates, with quantiles of the T\code{ny}-observation, scale and shape parameters
#' @param dat vector of threshold exceedances
#' @param s value of the stopping rule
#' @param p probability
#' @param N number of observations for the maximum
#' @param ny number of observations per year
maiq_fc.Vfun <- function(par, dat, s, p, N, ny){
  y <- dat[-length(dat)]
  ynu <- dat[length(dat)]
  stopifnot(length(par) == 3L, par[1:2] > 0)
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  V_qp_nu <- -(s*xi - xi*ynu)/((qp - s)*xi - sigma*(-N/log(p))^xi)
  V_sigma_nu <- (s*(-N/log(p))^xi - ynu*(-N/log(p))^xi)/((qp - s)*xi - sigma*(-N/log(p))^xi)
  V_xi_nu <- (s*sigma*xi^2*(-N/log(p))^xi*log(-N/log(p)) - s*sigma*xi*(-N/log(p))^xi - (sigma*xi^2*(-N/log(p))^xi*log(-N/log(p)) - sigma*xi*(-N/log(p))^xi + ((qp - s)*xi^2 - sigma*xi*(-N/log(p))^xi)*log(-((qp - s)*xi - sigma*(-N/log(p))^xi)/sigma))*ynu - ((2*qp*(-N/log(p))^xi - s*(-N/log(p))^xi)*sigma*xi - (qp^2 - qp*s)*xi^2 - sigma^2*(-N/log(p))^(2*xi))*log(-((qp - s)*xi - sigma*(-N/log(p))^xi)/sigma) + ((2*qp*(-N/log(p))^xi - s*(-N/log(p))^xi)*sigma*xi - (qp^2 - qp*s)*xi^2 - sigma^2*(-N/log(p))^(2*xi) + ((qp - s)*xi^2 - sigma*xi*(-N/log(p))^xi)*ynu)*log(-(qp*xi - xi*ynu - sigma*(-N/log(p))^xi)/sigma))/((qp - s)*xi^3 - sigma*xi^2*(-N/log(p))^xi)


  V_qp <- -sigma*(-(qp - y - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi + 1)*(ny/(-(qp - s - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) - ny/(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi))*((ny*(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma - ny*(-(qp - y - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma)/(ny/(-(qp - s - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) - ny/(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi)) - (ny*(-(qp - s - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma - ny*(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma)*(ny/(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) - ny/(-(qp - y - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi))/(ny/(-(qp - s - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) - ny/(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi))^2)/ny
  V_sigma <- sigma*(-(qp - y - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi + 1)*(ny/(-(qp - s - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) - ny/(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi))*((ny*(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma^2 + ((-N/log(p))^xi - 1)/sigma)/xi - ny*(-(qp - y - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((qp - y - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma^2 + ((-N/log(p))^xi - 1)/sigma)/xi)/(ny/(-(qp - s - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) - ny/(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi)) - (ny*(-(qp - s - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((qp - s - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma^2 + ((-N/log(p))^xi - 1)/sigma)/xi - ny*(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma^2 + ((-N/log(p))^xi - 1)/sigma)/xi)*(ny/(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) - ny/(-(qp - y - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi))/(ny/(-(qp - s - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi) - ny/(-(qp - u - sigma*((-N/log(p))^xi - 1)/xi)*xi/sigma + 1)^(1/xi))^2)/ny
  V_xi <- -(((qp - s)*xi - sigma*(-N/log(p))^xi)*u*(-((qp - s)*xi - sigma*(-N/log(p))^xi)/sigma)^(1/xi) - ((qp - u)*xi - sigma*(-N/log(p))^xi)*s*(-((qp - u)*xi - sigma*(-N/log(p))^xi)/sigma)^(1/xi) + (qp*xi - xi*y - sigma*(-N/log(p))^xi)*(s - u)*(-(qp*xi - xi*y - sigma*(-N/log(p))^xi)/sigma)^(1/xi) - (((qp - s)*xi - sigma*(-N/log(p))^xi)*(-((qp - s)*xi - sigma*(-N/log(p))^xi)/sigma)^(1/xi) - ((qp - u)*xi - sigma*(-N/log(p))^xi)*(-((qp - u)*xi - sigma*(-N/log(p))^xi)/sigma)^(1/xi))*y)*xi/((qp^2*xi^2 + s*u*xi^2 + sigma^2*(-N/log(p))^(2*xi) + (s*xi + u*xi)*sigma*(-N/log(p))^xi - (s*xi^2 + u*xi^2 + 2*sigma*xi*(-N/log(p))^xi)*qp)*((-((qp - s)*xi - sigma*(-N/log(p))^xi)/sigma)^(1/xi) - (-((qp - u)*xi - sigma*(-N/log(p))^xi)/sigma)^(1/xi)))


  cbind(c(V_qp, V_qp_nu), c(V_sigma, V_sigma_nu), c(V_xi, V_xi_nu))
}

maiq_fc.phi <- function(par, dat, s, p, N, ny, V){
  stopifnot(length(par) == 3L, par[1:2] > 0)
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  phi_a <- (xi + 1)/(qp*xi - xi*dat - sigma*(-N/log(p))^xi)
  rbind(phi_a) %*% V
}

maiq_fc.dphi <- function(par, dat, s, p, N, ny, V){
  stopifnot(length(par) == 3L, par[1:2] > 0)
  y <- dat #pointer
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  dphi_qp <- -(xi^2 + xi)/(qp^2*xi^2 + xi^2*y^2 - 2*qp*sigma*xi*(-N/log(p))^xi + sigma^2*(-N/log(p))^(2*xi) - 2*(qp*xi^2 - sigma*xi*(-N/log(p))^xi)*y)
  dphi_sigma <- (xi*(-N/log(p))^xi + (-N/log(p))^xi)/(qp^2*xi^2 + xi^2*y^2 - 2*qp*sigma*xi*(-N/log(p))^xi + sigma^2*(-N/log(p))^(2*xi) - 2*(qp*xi^2 - sigma*xi*(-N/log(p))^xi)*y)
  dphi_xi <- (sigma*xi*(-N/log(p))^xi*log(-N/log(p)) + sigma*(-N/log(p))^xi*(log(-N/log(p)) - 1) - qp + y)/(qp^2*xi^2 + xi^2*y^2 - 2*qp*sigma*xi*(-N/log(p))^xi + sigma^2*(-N/log(p))^(2*xi) - 2*(qp*xi^2 - sigma*xi*(-N/log(p))^xi)*y)
  rbind(dphi_qp, dphi_sigma, dphi_xi) %*% V
}


#' Standard log likelihood for the Maiquetia data
#'
#' We consider a nonhomogeneous Poisson process for threshold exceedances \code{y} above \code{u},
#' with intensity function \deqn{\Lambda(y) = n_y(1+\xi(y-\mu)/\sigma)^(-1/\xi),}
#' where \eqn{n_y} is the average number of exceedances per
#' year and a generalized extreme value likelihood for the
#' yearly maximum \code{z}.
#'
#' @param par vector of parameters, \code{p} quantile of the \code{T} year maximum distribution, scale and shape parameters
#' @param u threshold
#' @param y threshold exceedances
#' @param z yearly maximum observations
#' @param ny average number of exceedances per year
#' @param N total number of years for N-year maximum
#' @param ql quantile level
maiq_std.nll <- function(par, u, y, z, ny, N, ql){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  mu <- qp + sigma/xi*(1-(-N/log(ql))^xi)
  - suppressWarnings(sum(mev::pp.ll(par = c(mu, sigma, xi), dat = y, u = u, np = ny)) -
                       sum(mev::gev.ll(par = c(mu, sigma, xi), dat = z)))
}

#' Full conditional log likelihood for the Maiquetia data
#'
#' We consider a nonhomogeneous Poisson process for threshold exceedances \code{y} above \code{u},
#' with intensity function \deqn{\Lambda(y) = n_y(1+\xi(y-\mu)/\sigma)^(-1/\xi),}
#' where \eqn{n_y} is the average number of exceedances per
#' year and a generalized extreme value likelihood for the
#' yearly maximum \code{z}. This log likelihood function is obtained by
#' conditioning on the sets in which observations fall and the total counts.
#'
#' @inheritParams llmaiq_std
#' @param s value triggering the stopping rule
maiq_fc.nll <- function(par, u, s, y, z, ny, N, ql){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  mu <- qp + sigma/xi*(1-(-N/log(ql))^xi)
  Lambda_f <- function(x){ny*pmax(1+xi*(x-mu)/sigma,0)^(-1/xi)}
  Lambda_s <- Lambda_f(s)
  Lambda_u <- Lambda_f(u)
  if(any(isTRUE(all.equal(Lambda_s,0, check.attributes = FALSE)),
         isTRUE(all.equal(Lambda_s,0, check.attributes = FALSE)))){
    return(1e10)
  }
  ll <- -length(y)*log(sigma) - (1+1/xi)*sum(log(pmax(1+xi*(y-mu)/sigma, 0))) -
    (length(y)-1)*log(Lambda_u-Lambda_s) - log(Lambda_s) +
    sum(mev::gev.ll(par = c(mu, sigma, xi), dat = z))
  if(!is.numeric(ll)){
    ll <- 1e10
  }
  return(as.numeric(-ll))
}

#' Matrix of sufficient directions for the
#' generalized extreme value likelihood parametrized in terms of the
#' p quantile of the T observation maximum
#' @param par vector of parameters estimates, with quantiles of the T\code{ny}-observation, scale and shape parameters
#' @param z vector of yearly maximum
#' @param ql probability level for the quantile
#' @param N number of block for the maximum
gevq.Vfun <- function(par, z, N, ql){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  Nr <- -N/log(ql)
cbind(-1,
      (qp - z)/sigma,
      ((sigma*(Nr)^xi*log(Nr) - qp)*xi + xi*z + (qp*xi - xi*z - sigma*(Nr)^xi)*log(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma))/xi^2
)
}

gevq.phi <- function(par, z, N, ql, V){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
 phi_a <- as.numeric(((xi + 1)*(-(qp*xi - xi*z - sigma*(-N/log(ql))^xi)/sigma)^(1/xi) - 1)/((qp*xi - xi*z - sigma*(-N/log(ql))^xi)*(-(qp*xi - xi*z - sigma*(-N/log(ql))^xi)/sigma)^(1/xi)))
 rbind(phi_a) %*% V
}

gevq.dphi <- function(par, z, N, ql, V){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  Nr <- -N/log(ql)
  dphi_qp <- -((xi^2 + xi)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi) - xi - 1)/((qp^2*xi^2 + xi^2*z^2 - 2*qp*sigma*xi*(Nr)^xi + sigma^2*(Nr)^(2*xi) - 2*(qp*xi^2 - sigma*xi*(Nr)^xi)*z)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi))
  dphi_sigma <- ((sigma*xi*(Nr)^xi + sigma*(Nr)^xi)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi) - sigma*(Nr)^xi - qp + z)/((qp^2*sigma*xi^2 + sigma*xi^2*z^2 - 2*qp*sigma^2*xi*(Nr)^xi + sigma^3*(Nr)^(2*xi) - 2*(qp*sigma*xi^2 - sigma^2*xi*(Nr)^xi)*z)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi))
  dphi_xi <- -((sigma*(Nr)^xi*log(Nr) - qp)*xi^2 + (sigma*(Nr)^xi*log(Nr) - qp)*xi + (xi^2 + xi)*z - (sigma*xi^3*(Nr)^xi*log(Nr) + (sigma*(Nr)^xi*(log(Nr) - 1) - qp)*xi^2 + xi^2*z)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi) + (qp*xi - xi*z - sigma*(Nr)^xi)*log(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma))/((qp^2*xi^4 + xi^4*z^2 - 2*qp*sigma*xi^3*(Nr)^xi + sigma^2*xi^2*(Nr)^(2*xi) - 2*(qp*xi^4 - sigma*xi^3*(Nr)^xi)*z)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi))
  rbind(dphi_qp, dphi_sigma, dphi_xi) %*% V
}

#' Partial conditional log likelihood for the Maiquetia data
#'
#' We consider a non homogeneous Poisson process for threshold exceedances \code{y} above \code{u},
#' with intensity function \deqn{\Lambda(y) = n_y(1+\xi(y-\mu)/\sigma)^(-1/\xi),}
#' where \eqn{n_y} is the average number of exceedances per
#' year and a generalized extreme value likelihood for the
#' yearly maximum \code{z}. This log likelihood function is obtained by
#' conditioning on the number of observations and the fact that the last observation exceeds \code{s}
#'
#' @inheritParams llmaiq_std
#' @param s value triggering the stopping rule
# nllmaiq_pc <- function(par, u, s, y, z, ny, N, ql){
#   qp <- par[1]; sigma <- par[2]; xi <- par[3]
#   mu <- qp + sigma/xi*(1-(-N/log(ql))^xi)
#   Lambda_f <- function(x){pmax(0, ny*(1+xi*(x-mu)/sigma)^(-1/xi))}
#   Lambda_s <- Lambda_f(s)
#   if(isTRUE(all.equal(Lambda_s,0, check.attributes = FALSE))){
#     return(1e10)
#   }
#   ll <- -length(y)*log(sigma) - (1+1/xi)*sum(log(pmax(1+xi*(y-mu)/sigma, 0))) -
#     - log(Lambda_s) +  sum(mev::gev.ll(par = c(mu, sigma, xi), dat = z))
#   return(as.numeric(-ll))
# }

#' Inequality constraints for optimization
#'
#' This function is used for the Maiquetia
#' full conditional likelihood with stopping rule
#' @inheritParams llmaiq_std
#' @inheritParams llmaiq_fc
#' @keywords internal
constraints <- function(par, u, s, y, z, ny, N, ql){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  mu <- qp + sigma/xi*(1-(-N/log(ql))^xi)
  xv <- ifelse(xi > 0, min(c(u,s,y,z)), max(c(u,s,y,z)))
  as.numeric(c(sigma, 1 - xi, xi + 1, sigma + xi*(xv-mu)))
}


# Profile log likelihood for pth quantile of T max
profile_qp_fc_maiq <- function(psi, start, mod = c("profile","tem"),
                               y, z, s, ny, N = 50, ql = 0.5, correction = TRUE,
                               threshold = 0, plot = TRUE, ...){
  names_mle <- c("Nquant","scale","shape")
  param <- "Nquant"
  mod <- match.arg(mod, c("profile", "tem"), several.ok = TRUE)
  ymin <- min(y)
  ymax <- max(y)
  zmin <- min(z)
  zmax <- max(z)
  #' Wrapper for inequality constraints
  #'
  #' Wrapper around \code{constraints} function
  #' @keywords internal
  constraints_fn <- function(par){
    constraints(par = par,
                z = z,
                s = s,
                y = y,
                u = threshold, ny = ny, N = N, ql = ql)}

  #' Full conditional log likelihood for Maiquetia data
  #'
  #' Wrapper around \code{maiq_fc.nll} for \code{alabama::auglag},
  #' which complains about parameters in \code{hin} or
  #' \code{fn} without the same arguments.
  #'
  #' @keywords internal
  maiq_fc_nll_fn <- function(par){
    maiq_fc.nll(par = par, z = z,
               s = s,
               y = y,
               u = threshold, ny = ny, N = N, ql = ql)
  }
  mle_opt <- Rsolnp::solnp(pars = start,
                           fun = maiq_fc_nll_fn,
                           ineqfun = constraints_fn,
                           ineqLB = rep(0,4),
                           ineqUB = rep(Inf, 4),
                           control = list(trace = 0))
  mle <- mle_opt$par
  names(mle) <- names_mle
  mle_opt$hessian <- numDeriv::hessian(fun = maiq_fc_nll_fn, x = mle, method.args=list(eps=1e-4, d=1e-5, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))
  std_error <- sqrt(diag(solve(mle_opt$hessian)))
  # Extract the components, notably V for model `tem`. Keep other components for optimization
  if("tem" %in% mod){
    V1 <- maiq_fc.Vfun(par = mle, dat = y, s = s, N = N, p = ql, ny = ny)
    V2 <- gevq.Vfun(par = mle, z = z, N = N, ql = ql)
    V <- rbind(V1, V2)
  } else{
    V <- NULL
  }
  maxll <- -mle_opt$value[length(mle_opt$value)]
  constr.mle.Nquant <- function(psi, start = mle[-1]) {
      suppressWarnings(Rsolnp::solnp(
        pars = start,
        fun = function(lambda){maiq_fc.nll(par = c(psi, lambda), y = y, z = z, s = s, u = threshold, ny = ny, N = N, ql = ql) },
        ineqfun = function(lambda){constraints_fn(c(psi, lambda))},
        ineqUB = rep(Inf, 4),
        ineqLB = rep(0, 4),
        control = list(trace=0))$pars)
    # alabama::auglag(
    #     par = start,
    #     fn = function(lambda, psi) {maiq_std.nll(par = c(psi, lambda), y = y, z = z, u = threshold, ny = ny, N = N, ql = ql) },
    #     hin = function(lambda, psi){constraints_fn(c(psi, lambda))},
    #     control.outer = list(method = "BFGS", trace = FALSE),
    #     control.optim = list(parscale = abs(mle[-1]),reltol = 1e-10, maxit = 1e4),
    #     psi = psi)$par
    }
  # Missing psi vector
  if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
    #compute profile log-lik on a grid left and right of the MLE
    psi <- unique(seq(from = mle[1]-2*std_error, to = mle[1]+3*std_error, length.out = 101L))
    psi <- psi[psi > 0]
  }
  pars <- matrix(nrow = length(psi), ncol = 3)
  pars[,1] <- psi
  mid <- which.min(sapply(psi, function(p){abs(p-mle[1])}))
  pars[mid, 2:3] <- constr.mle.Nquant(psi = psi[mid])
  for(i in (mid-1):1){
    pars[i,2:3] <- constr.mle.Nquant(psi = psi[i], start = pars[i+1,2:3])
  }
  for(i in (mid+1):length(psi)){
    pars[i,2:3] <- constr.mle.Nquant(psi = psi[i], start = pars[i-1,2:3])
  }
  # Profile log likelihood values for psi
  profll <- apply(pars, 1, function(par) {
    -maiq_fc.nll(par = par, y = y, z = z, s = s, u = threshold, ny = ny, N = N, ql = ql)
  })
  r <- sign(mle[param] - psi) * sqrt(2 * pmax(0, maxll - profll))
  #Sometimes mle difference is -1e-6, which gives NA...
  if ("tem" %in% mod) {
    phi_mle <- maiq_fc.phi(par = mle, dat = y, N = N, V = V1, p = ql, s = s) +
                 gevq.phi(par = mle, z = z, N = N, V = V2, ql = ql)
    dphi_mle <- maiq_fc.dphi(par = mle, dat = y, N = N, V = V1, p = ql, s = s) +
      gevq.dphi(par = mle, z = z, N = N, V = V2, ql = ql)
    det_dphi_mle <- det(dphi_mle)
    phi_par <- function(par){maiq_fc.phi(par = par, dat = y, N = N, V = V1, p = ql, s = s) +
      gevq.phi(par = par, z = z, N = N, V = V2, ql = ql)}
    dphi_par <- function(par){maiq_fc.dphi(par = par, dat = y, N = N, V = V1, p = ql, s = s) +
        gevq.dphi(par = par, z = z, N = N, V = V2, ql = ql)}

    q2num <- apply(pars, 1, function(par) {
      det(rbind(c(phi_mle - phi_par(par)),
                dphi_par(par)[-1, ]))
    })
    if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
      warning("Correction factor and likelihood root are of opposite sign - check output")
    }
    # Obtain information matrix through numerical differentiation
    maiq.infomat <- function(par){ #negative of hessian of ll
      numDeriv::hessian(func = function(x){ #loglik already negated
        maiq_fc.nll(par = x, u = threshold, s = s, y = y, z = z, ny = ny, N = N, ql = ql)
        }, x = par, method.args=list(eps=1e-4, d=1e-5, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))}
    logq <- apply(pars, 1, function(par) {
      -0.5 * log(det(maiq.infomat(par)[-1, -1]))
    }) + log(abs(q2num))
    qmlecontrib <- -log(abs(det_dphi_mle)) + 0.5 *
      log(det(mle_opt$hessian))
    logq <- logq + qmlecontrib
    qcor <- sign(q2num) * exp(logq)
    rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r)))/r)
    # Not working so far...
    tem.max.opt <- function(psi){
      para <- c(psi, constr.mle.Nquant(psi, start = pars[which.min(abs(psi - pars[,1]))[1],2:3])) #better starting values?
      pll <- -maiq_fc.nll(par = para, z = z, y = y, ny = ny, u = threshold, s = s, ql = ql, N = N)
      rs <- 2 * pmax(0, maxll - pll)
      logq <- -0.5 * log(det(maiq.infomat(par = para)[-1,-1])) + qmlecontrib +
        log(abs(det(rbind(c(phi_mle) -c(phi_par(par = para)),dphi_par(par = para)[-1, ]))))
      rs + (logq - log(sqrt(abs(rs))))
    }
    opt.tem <- optim(par = mle[1]+2, fn = tem.max.opt, method = "Brent", lower = max(1e-05, mle[1] - std_error[1]), upper = mle[1] + std_error[1], control = list(abstol = 1e-10))
    tem.max <- opt.tem$par
    tem.rel <- opt.tem$value
  }
  # Return profile likelihood and quantities of interest (modified likelihoods)
  colnames(pars) <- names(mle)
  ans <- list(mle = mle, pars = pars, psi.max = as.vector(mle[param]),
              param = param, std_error = std_error,
              psi = psi, pll = profll, maxpll = maxll, r = r, V = V)
  ans$normal <- c(ans$psi.max, ans$std_error)
  if ("tem" %in% mod) {
    ans$q <- qcor
    ans$rstar <- rstar
    ans$tem.psimax <- as.vector(tem.max)
    ans$tem.reliability <- tem.rel
    if (correction && length(psi) > 10) {
      ans <- spline.corr(ans)
    }
  }
  if ("tem" %in% mod) {
    class(ans) <- c("eprof", "fr")
  } else {
    class(ans) <- "eprof"
  }
  if(plot){
    plot(ans)
  }
  return(invisible(ans))

}



# Profile log likelihood for pth quantile of T max
profile_qp_std_maiq <- function(psi, start, mod = c("profile","tem"),
                                y, z, ny, N = 50, ql = 0.5,
                                correction = TRUE,
                                threshold = 0, plot = TRUE, ...){
  names_mle <- c("Nquant","scale","shape")
  param <- "Nquant"
  mod <- match.arg(mod, c("profile", "tem"), several.ok = TRUE)
  ymin <- min(y)
  ymax <- max(y)
  zmin <- min(z)
  zmax <- max(z)

  #' Wrapper for inequality constraints
  #'
  #' Wrapper around \code{constraints} function
  #' @keywords internal
  constraints_fn <- function(par){
    qp <- par[1]; sigma <- par[2]; xi <- par[3]
    mu <- qp + sigma/xi*(1-(-N/log(ql))^xi)
    xv <- ifelse(xi > 0, min(c(y, z, threshold)), max(c(threshold, y, z)))
    as.numeric(c(sigma/1000, 1 - xi, xi + 1, sigma + xi*(xv-mu)))
  }

  #' Full conditional log likelihood for Maiquetia data
  #'
  #' Wrapper around \code{maiq_std.nll} for \code{alabama::auglag},
  #' which complains about parameters in \code{hin} or
  #' \code{fn} without the same arguments.
  #'
  #' @keywords internal
  maiq_std_nll_fn <- function(par){
    maiq_std.nll(par = par, z = z, y = y,
                u = threshold, ny = ny, N = N, ql = ql)
  }

  mle_opt <- Rsolnp::solnp(pars = start,
                           fun = maiq_std_nll_fn,
                           ineqfun = constraints_fn,
                           ineqLB = rep(0, 4),
                           ineqUB = rep(Inf, 4),
                           control = list(trace = 0))
  mle <- mle_opt$pars
  names(mle) <- names_mle
  mle_opt$hessian <- numDeriv::hessian(fun = maiq_std_nll_fn, x = mle, method.args=list(eps=1e-4, d=1e-5, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))
  std_error <- sqrt(diag(solve(mle_opt$hessian)))
  # Extract the components, notably V for model `tem`. Keep other components for optimization
  if("tem" %in% mod){
     V1 <- nhppq.Vfun(par = mle, dat = y, N = N, ql = ql, ny = ny)
     V2 <- gevq.Vfun(par = mle, z = z, N = N, ql = ql)
     V <- rbind(V1, V2)
  } else{
    V <- NULL
  }
  maxll <- -mle_opt$value[length(mle_opt$value)]

  constr.mle.Nquant <- function(psi, start = mle[-1]) {
    suppressWarnings(
      Rsolnp::solnp(pars = start,
                     fun = function(lambda) {maiq_std.nll(par = c(psi, lambda), y = y, z = z, u = threshold, ny = ny, N = N, ql = ql) },
                     ineqfun = function(lambda){constraints_fn(c(psi, lambda))},
                    ineqLB = rep(0,4), ineqUB = rep(Inf, 4), control =list(trace=0))$pars
      # auglag(x0 = start,
      #                fn = function(lambda) {maiq_std.nll(par = c(psi, lambda), y = y, z = z, u = threshold, ny = ny, N = N, ql = ql) },
      #                hin = function(lambda){constraints_fn(c(psi, lambda))},
      #                localsolver = "BOBYQA")$par
      # alabama::auglag(
      #   par = start,
      #   fn = function(lambda, psi) {maiq_std.nll(par = c(psi, lambda), y = y, z = z, u = threshold, ny = ny, N = N, ql = ql) },
      #   hin = function(lambda, psi){constraints_fn(c(psi, lambda))},
      #   control.outer = list(method = "BFGS", trace = FALSE),
      #   control.optim = list(parscale = abs(mle[-1]),reltol = 1e-10, maxit = 1e4),
      #   psi = psi)$par
    )
  }
  # Missing psi vector
  if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
    #compute profile log-lik on a grid left and right of the MLE
    psi <- unique(seq(from = mle[1]-2.5*std_error, to = mle[1]+3*std_error, length.out = 101L))
    psi <- psi[psi > 0]
  }
  pars <- matrix(nrow = length(psi), ncol = 3)
  pars[,1] <- psi
  mid <- which.min(sapply(psi, function(p){abs(p-mle[1])}))
  pars[mid, 2:3] <- constr.mle.Nquant(psi = psi[mid])
  for(i in (mid-1):1){
    pars[i,2:3] <- constr.mle.Nquant(psi = psi[i], start = pars[i+1,2:3])
  }
  for(i in mid:length(psi)){
    pars[i,2:3] <- constr.mle.Nquant(psi = psi[i], start = pars[i-1,2:3])
  }
  # par(mfrow = c(1,2), mar = c(4,4,1,1), bty = "l")
  # plot(pars[,1], pars[,2], type = "l")
  # plot(pars[,1], pars[,3], type = "l")
  # Profile log likelihood values for psi
  profll <- apply(pars, 1, function(par) {
    -maiq_std.nll(par = par, y = y, z = z, u = threshold, ny = ny, N = N, ql = ql)
  })
  r <- sign(mle[param] - psi) * sqrt(2 * pmax(0, maxll - profll))
  #Sometimes mle difference is -1e-6, which gives NA...
  if ("tem" %in% mod) {
    maiq_std_phi <- function(para){
      nhppq.phi(par = para, dat = y, ql = ql, N = N, ny = ny, V = V1) +
        gevq.phi(par = para, z = z, N= N, ql = ql, V = V2)
    }
    maiq_std_dphi <- function(para){
      nhppq.dphi(par = para, dat = y, ql = ql, N = N, ny = ny, V = V1) +
        gevq.dphi(par = para, z = z, N= N, ql = ql, V = V2)
    }
    phi.mle <- maiq_std_phi(para = mle)
    q2num <- apply(pars, 1, function(par) {
      det(rbind(c(phi.mle - maiq_std_phi(para = par)),
                maiq_std_dphi(para = par)[-1, ]))
    })
    det_dphi_mle <- det(maiq_std_dphi(para = mle))
    if (isTRUE(any(sign(q2num)/sign(det_dphi_mle) * sign(r) < 0, na.rm = TRUE))) {
      warning("Correction factor and likelihood root are of opposite sign - check output")
    }
    # Obtain information matrix through numerical differentiation
    maiq_std_infomat <- function(par){ #negative of hessian of ll
      numDeriv::hessian(func = function(x){ #loglik already negated
        maiq_std.nll(par = x, u = threshold, y = y, z = z, ny = ny, N = N, ql = ql)
      }, x = par, method.args=list(eps=1e-4, d=1e-5, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))}

    logq <- apply(pars, 1, function(par) {
      -0.5 * log(det(maiq_std_infomat(par = par)[-1, -1]))
    }) + log(abs(q2num))
    qmlecontrib <- -log(abs(det_dphi_mle)) + 0.5 *
      log(det(maiq_std_infomat(par = mle)))
    logq <- logq + qmlecontrib
    qcor <- sign(q2num) * exp(logq)
    rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r)))/r)
    tem.max.opt <- function(psi) {
      para <- c(psi, constr.mle.Nquant(psi, start = pars[which.min(abs(psi - pars[,1]))[1],2:3]))
      pll <- -maiq_std.nll(par = para, u = threshold, y = y, z = z, ny = ny, N = N, ql = ql)
      rs <- 2 * pmax(0, maxll - pll)
      logq <- -0.5 * log(det(maiq_std_infomat(par = para)[-1,-1])) +
        qmlecontrib +
        log(abs(det(rbind(c(phi.mle - maiq_std_phi(para = para)), maiq_std_dphi(para = para)[-1, ]))))
      rs + (logq - log(sqrt(abs(rs))))
    }
    # Technically, this equation should solve rstar=0, so solving r*rstar=0 gives the second root if rstar != r
    # #Bespoke implementation: check only to the right of the MLE
    opt.tem <- optim(par = mle[1]+10, fn = tem.max.opt, method = "Brent",  lower = max(1e-05, mle[1] + 0.01* std_error[1]), upper = mle[1] +1.5* std_error[1], control = list(abstol = 1e-10))
    tem.max <- opt.tem$par
    tem.rel <- opt.tem$value
  }
  # Return profile likelihood and quantities of interest (modified likelihoods)
  colnames(pars) <- names(mle)
  ans <- list(mle = mle, pars = pars, psi.max = as.vector(mle[param]),
              param = param, std_error = std_error,
              psi = psi, pll = profll, maxpll = maxll, r = r, V = V)
  if ("tem" %in% mod) {
    ans$q <- qcor
    ans$rstar <- rstar
    ans$tem.psimax <- as.vector(tem.max)
    ans$normal <- c(ans$psi.max, ans$std_error)
    ans$tem.reliability <- tem.rel
    if (correction && length(psi) > 10) {
      ans <- spline.corr(ans)
    }
  }
  if ("tem" %in% mod) {
    class(ans) <- c("eprof", "fr")
  } else {
    class(ans) <- "eprof"
  }
  if(plot){
    plot(ans)
  }
  return(invisible(ans))

}
