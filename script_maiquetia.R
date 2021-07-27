#####################################################################
## Nonhomogeneous Poisson process likelihood and HOA approximation ##
##  parametrized in terms of pth quantile of the T-year maximum    ##
#####################################################################



#' Observed information matrix of the Poisson process (quantile)
#'
#' The inhomogeneous Poisson process likelihood is parametrized so that its first argument is the \code{q} quantile of the \code{N} year maximum, assuming that there are on average \code{np} observations by year.
#' @note Currently, the function does not support the case when the shape is zero.
#'
#' @param par vector of parameters (quantile, scale, shape)
#' @param dat vector of data
#' @param method string; only observed information (\code{obs})
#' @param u threshold
#' @param np number of periods of observations. This is a \emph{post hoc} adjustment for the intensity so that the parameters of the model coincide with those of a generalized extreme value distribution with block size \code{length(dat)/np}.
#' @param q probability, corresponding to \code{q}th quantile of the N-block maximum
#' @param N the first parameter is the \code{q} quantile of the \code{N}-block maximum
nhpp.qp.infomat <- function(par, dat, method = "obs", N, u, np, q = 0.5){
  qp <- par[1]
  sigma = par[2]
  xi = par[3]
  if(isTRUE(all.equal(xi, 0, check.attributes = FALSE))){
    stop("Observed information matrix not implemented for the case xi=0.")
  }
  mNlogq = -N/log(q)
  logmNlogq = log(mNlogq)
  infomat <- matrix(0, nrow = 3, ncol = 3)
  infomat[1,1] <- np*(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)^(1/xi - 2)*xi*(1/xi - 1)/(sigma^2*(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)^(2/xi)) -
    2*np*(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)^(2/xi - 2)/(sigma^2*(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)^(3/xi)) + xi*(1 + xi)*sum(1/(sigma^2*((qp - dat - sigma*(mNlogq^xi - 1)/xi)*xi/sigma - 1)^2))


  infomat[1,2] <- infomat[2,1] <-
    (np*sigma*mNlogq^xi + np*qp - np*u)/((sigma^3*mNlogq^(2*xi) + (qp^2*sigma - 2*qp*sigma*u + sigma*u^2)*xi^2 - 2*(qp*sigma^2*mNlogq^xi - sigma^2*u*mNlogq^xi)*xi)*(-((qp - u)*xi - sigma*mNlogq^xi)/sigma)^(1/xi)) -
    sum((xi*mNlogq^xi + mNlogq^xi)/((qp^2 - 2*qp*dat + dat^2)*xi^2 + sigma^2*mNlogq^(2*xi) - 2*(qp*sigma*mNlogq^xi - sigma*dat*mNlogq^xi)*xi)
    )

  infomat[1,3] <- infomat[3,1] <- ((np*sigma*mNlogq^xi*logmNlogq - np*qp + np*u)*xi^2 + (np*sigma*mNlogq^xi*logmNlogq - np*qp + np*u)*xi -
(np*sigma*mNlogq^xi - (np*qp - np*u)*xi)*log(-((qp - u)*xi - sigma*mNlogq^xi)/sigma))/(((qp^2 - 2*qp*u + u^2)*xi^4 + sigma^2*xi^2*mNlogq^(2*xi) - 2*(qp*sigma*mNlogq^xi -
 sigma*u*mNlogq^xi)*xi^3)*(-((qp - u)*xi - sigma*mNlogq^xi)/sigma)^(1/xi)) + sum(-(sigma*xi*mNlogq^xi*logmNlogq + sigma*mNlogq^xi*(logmNlogq - 1) - qp + dat)/
((qp^2 - 2*qp*dat + dat^2)*xi^2 + sigma^2*mNlogq^(2*xi) - 2*(qp*sigma*mNlogq^xi - sigma*dat*mNlogq^xi)*xi))
  infomat[2,3] <- infomat[3,2] <- sum((qp*sigma*mNlogq^xi*(logmNlogq - 1) - qp^2 - (sigma*mNlogq^xi*(logmNlogq - 1) - 2*qp)*dat - dat^2 + (qp*sigma*mNlogq^xi*logmNlogq - sigma*dat*mNlogq^xi*logmNlogq)*xi)/(sigma^3*mNlogq^(2*xi) + (qp^2*sigma - 2*qp*sigma*dat + sigma*dat^2)*xi^2 - 2*(qp*sigma^2*mNlogq^xi - sigma^2*dat*mNlogq^xi)*xi)) -
    ((np*qp*sigma*mNlogq^xi*logmNlogq - np*qp^2 - np*u^2 - (np*sigma*mNlogq^xi*logmNlogq - 2*np*qp)*u)*xi^2 + (np*qp*sigma*mNlogq^xi*logmNlogq - np*qp^2 - np*u^2 - (np*sigma*mNlogq^xi*logmNlogq - 2*np*qp)*u)*xi - (np*qp*sigma*mNlogq^xi - np*sigma*u*mNlogq^xi - (np*qp^2 - 2*np*qp*u + np*u^2)*xi)*log(-((qp - u)*xi -
sigma*mNlogq^xi)/sigma))/((sigma^3*xi^2*mNlogq^(2*xi) + (qp^2*sigma - 2*qp*sigma*u + sigma*u^2)*xi^4 - 2*(qp*sigma^2*mNlogq^xi - sigma^2*u*mNlogq^xi)*xi^3)*(-((qp - u)*xi - sigma*mNlogq^xi)/sigma)^(1/xi))

  infomat[2,2] <- sum(-((qp^2 - 2*qp*dat + dat^2)*xi^2 - 2*qp*sigma*mNlogq^xi + 2*sigma*dat*mNlogq^xi - (2*qp*sigma*mNlogq^xi - qp^2 - 2*(sigma*mNlogq^xi - qp)*dat - dat^2)*xi)/(sigma^4*mNlogq^(2*xi) + (qp^2*sigma^2 - 2*qp*sigma^2*dat +
sigma^2*dat^2)*xi^2 - 2*(qp*sigma^3*mNlogq^xi - sigma^3*dat*mNlogq^xi)*xi)) + np*(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*((qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma^2 + (mNlogq^xi - 1)/sigma)^2*(1/xi - 1)/xi - 2*np*(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*((qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma^2 + (mNlogq^xi - 1)/sigma)^2/xi^2 -
    2*np*(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*((qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma^3 + (mNlogq^xi - 1)/sigma^2)/xi + length(dat)/sigma^2

  infomat[3,3] <- sum(((sigma*mNlogq^xi*logmNlogq/xi - sigma*(mNlogq^xi - 1)/xi^2)*xi/sigma - (qp - dat - sigma*(mNlogq^xi - 1)/xi)/sigma)^2*(1/xi + 1)/((qp - dat - sigma*(mNlogq^xi - 1)/xi)*xi/sigma - 1)^2 + ((sigma*mNlogq^xi*logmNlogq^2/xi - 2*sigma*mNlogq^xi*logmNlogq/xi^2 + 2*sigma*(mNlogq^xi - 1)/xi^3)*xi/sigma +
2*(sigma*mNlogq^xi*logmNlogq/xi - sigma*(mNlogq^xi - 1)/xi^2)/sigma)*(1/xi + 1)/((qp - dat - sigma*(mNlogq^xi - 1)/xi)*xi/sigma - 1) - 2*((sigma*mNlogq^xi*logmNlogq/xi - sigma*(mNlogq^xi - 1)/xi^2)*xi/sigma -
(qp - dat - sigma*(mNlogq^xi - 1)/xi)/sigma)/(((qp - dat - sigma*(mNlogq^xi - 1)/xi)*xi/sigma - 1)*xi^2) - 2*log(-(qp - dat - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)/xi^3) -np*(((sigma*mNlogq^xi*logmNlogq/xi - sigma*(mNlogq^xi - 1)/xi^2)*xi/sigma -
(qp - u - sigma*(mNlogq^xi - 1)/xi)/sigma)/(((qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma - 1)*xi) +
log(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)/xi^2)^2/(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)^(1/xi) - np*(((sigma*mNlogq^xi*logmNlogq/xi - sigma*(mNlogq^xi - 1)/xi^2)*xi/sigma -
(qp - u - sigma*(mNlogq^xi - 1)/xi)/sigma)^2/(((qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma - 1)^2*xi) + ((sigma*mNlogq^xi*logmNlogq^2/xi - 2*sigma*mNlogq^xi*logmNlogq/xi^2 + 2*sigma*(mNlogq^xi - 1)/xi^3)*xi/sigma + 2*(sigma*mNlogq^xi*logmNlogq/xi - sigma*(mNlogq^xi - 1)/xi^2)/sigma)/(((qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma - 1)*xi) -
2*((sigma*mNlogq^xi*logmNlogq/xi - sigma*(mNlogq^xi - 1)/xi^2)*xi/sigma - (qp - u - sigma*(mNlogq^xi - 1)/xi)/sigma)/(((qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma - 1)*xi^2) -
2*log(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)/xi^3)/(-(qp - u - sigma*(mNlogq^xi - 1)/xi)*xi/sigma + 1)^(1/xi)

  return(-infomat)
}


#' Sufficient directions for nonhomogeneous Poisson process
#'
#' Matrix of sufficient directions for the nonhomogeneous
#' Poisson process likelihood parametrized in terms of the
#' p quantile of the T observation maximum
#' @param par vector of parameters estimates, with quantiles of the T\code{ny}-observation, scale and shape parameters
#' @param dat vector of threshold exceedances
#' @param q probability level of the quantile
#' @param N number of observations for the maximum
nhpp.qp.Vfun <- function(par, dat, q, N, ny){
  y <- dat
  stopifnot(length(par) == 3L, par[1:2] > 0)
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  cst <- sigma*((-N/log(q))^xi - 1)/xi
  V_qp <- ny*(-(qp - y - cst)*xi/sigma + 1)^(-1/xi - 2)*(1+xi)/sigma^2
  V_sigma <- -ny*(-(qp - y - cst)*xi/sigma + 1)^(-1/xi - 2)*((qp - y - cst)*xi/sigma^2 + ((-N/log(q))^xi - 1)/sigma)*(1/xi + 1)/sigma - ny*(-(qp - y - cst)*xi/sigma + 1)^(-1/xi - 1)/sigma^2
  V_xi <- ny*(-(qp - y - cst)*xi/sigma + 1)^(-1/xi - 1)*(((sigma*(-N/log(q))^xi*log(-N/log(q))/xi - cst/xi)*xi/sigma - (qp - y - cst)/sigma)*(1/xi + 1)/((qp - y - cst)*xi/sigma - 1) + log(-(qp - y - cst)*xi/sigma + 1)/xi^2)/sigma
  cbind(V_qp, V_sigma, V_xi)
}

#' Directional derivative of the log likelihood
#'
#' This term is \code{V} times the sample
#' space derivative of the log-likelihood
#'
#' @inheritParams nhpp.qp.Vfun
#' @inheritParams nhpp.qp.dphi
nhpp.qp.phi <- function(par, dat, q, N, ny, V){
  y <- dat
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  phi_a <- sigma*log(ny*(-(qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma)/(ny*(-(qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1))
  rbind(phi_a) %*% V
}

#' Mixed derivative of log likelihood
#'
#' This term is \code{V} times the sample
#' space derivative of the log-likelihood
#' @inheritParams nhpp.qp.Vfun
#' @param V an \code{n} by \code{p} matrix whose columns span the sufficient directions
nhpp.qp.dphi <- function(par, dat, q, N, ny, V){
  y <- dat
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  cst <- sigma*((-N/log(q))^xi - 1)/xi
  cst2 <- sigma*(-N/log(q))^xi
  dphi_qp <- (-(qp - y -cst)*xi/sigma + 1)^(1/xi)*(xi + 1)/ny

  dphi_sigma <- -(cst2 + qp - y)*(-(qp*xi - xi*y - cst2)/sigma)^(1/xi)/(ny*sigma)
  dphi_xi <- sigma*(-(qp - y - cst)*xi/sigma + 1)^(1/xi + 1)*
  (((cst2*log(-N/log(q))/xi - cst/xi)*xi/sigma - (qp - y -cst)/sigma)*(1/xi + 1)/((qp - y - cst)*xi/sigma - 1) + log(-(qp - y -cst)*xi/sigma + 1)/xi^2)/ny

  rbind(dphi_qp, dphi_sigma, dphi_xi) %*% V
}

#' Negative log likelihood of the inhomogeneous Poisson process
#' @inheritParams nhpp.qp.Vfun
nhpp.nll <- function(par, u, y, ny, N, q){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  mu <- qp + sigma/xi*(1-(-N/log(q))^xi)
  - suppressWarnings(sum(mev::pp.ll(par = c(mu, sigma, xi), dat = y, u = u, np = ny)))
}


# Profile log likelihood for pth quantile of T max
nhpp.qp.profile <- function(psi,
                            start,
                            mod = c("profile","tem"),
                            y,
                            ny = 1,
                            N = 50,
                            q = 0.5,
                            correction = FALSE,
                            threshold = 0,
                            plot = FALSE, ...){
  names_mle <- c("Nquant","scale","shape")
  param <- "Nquant"
  mod <- match.arg(mod,
                   c("profile", "tem"),
                   several.ok = TRUE)
  #' Wrapper for inequality constraints
  #'
  #' Wrapper around \code{constraints} function
  #' @keywords internal
  constraints.fn <- function(par){
    qp <- par[1];
    sigma <- par[2];
    xi <- par[3]
    mu <- qp + sigma/xi*(1-(-N/log(q))^xi)
    xv <- ifelse(xi > 0,
                 min(c(y, threshold)),
                 max(c(threshold, y)))
    as.numeric(c(sigma/1000,
                 1 - xi,
                 xi + 1,
                 sigma + xi*(xv-mu)))
  }

  nhpp.nlls <- function(par){
    nhpp.nll(par = par,
             y = y,
             u = threshold,
             ny = ny,
             N = N,
             q = q)
  }

  mle_opt <- suppressWarnings(
    Rsolnp::solnp(pars = start,
                  fun = nhpp.nlls,
                  ineqfun = constraints.fn,
                  ineqLB = rep(0, 4),
                  ineqUB = rep(Inf, 4),
                  control = list(trace = 0))
  )
  mle <- mle_opt$pars
  names(mle) <- names_mle
  mle_opt$hessian <-
    nhpp.qp.infomat(par = mle,
                    dat = y,
                    u = threshold,
                    np = ny,
                    N = N,
                    q = q)
  std_error <- sqrt(diag(solve(mle_opt$hessian)))
  # Extract the components, notably V for model `tem`. Keep other components for optimization
  if("tem" %in% mod){
    V <- nhpp.qp.Vfun(par = mle, dat = y, N = N, q = q, ny = ny)
#     # Check the sufficient directions are correct.
    intensQ <- function(par, dat, q, N, ny){
      stopifnot(length(par) == 3L, par[1:2] > 0)
      qp <- par[1]
      sigma <- par[2]
      xi <- par[3]
      ny*(-(qp - dat - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma
    }
# numDeriv::jacobian(intensQ,
#                    x = mle,
#                    dat = y,
#                    N = N,
#                    q = q,
#                    ny = ny)
  } else{
    V <- NULL
  }
  maxll <- -mle_opt$value[length(mle_opt$value)]

  constr.mle.Nquant <- function(psi, start = mle[-1]) {
    suppressWarnings(
      Rsolnp::solnp(pars = start,
                    fun = function(lambda) {
                      nhpp.nll(par = c(psi, lambda),
                               y = y,
                               u = threshold,
                               ny = ny,
                               N = N,
                               q = q)
                      },
                    ineqfun = function(lambda){
                      constraints.fn(c(psi, lambda))},
                    ineqLB = rep(0,4),
                    ineqUB = rep(Inf, 4),
                    control = list(trace = 0))$pars
    )
  }
  # Missing psi vector
  if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
    #compute profile log-lik on a grid left and right of the MLE
    psi <- unique(seq(from = mle[1]-2.5*std_error,
                      to = mle[1]+3*std_error,
                      length.out = 101L))
    psi <- psi[psi > 0]
  }
  pars <- matrix(nrow = length(psi), ncol = 3)
  pars[,1] <- psi
  mid <- which.min(sapply(psi, function(p){
    abs(p-mle[1])}))
  pars[mid, 2:3] <- constr.mle.Nquant(psi = psi[mid])
  for(i in (mid-1):1){
    pars[i,2:3] <- constr.mle.Nquant(psi = psi[i],
                                     start = pars[i+1,2:3])
  }
  for(i in mid:length(psi)){
    pars[i,2:3] <- constr.mle.Nquant(psi = psi[i],
                                     start = pars[i-1,2:3])
  }
  profll <- apply(pars, 1, function(par) {
    -nhpp.nll(par = par,
              y = y,
              u = threshold,
              ny = ny,
              N = N,
              q = q)
  })
  r <- sign(mle[param] - psi) * sqrt(2 * pmax(0, maxll - profll))
  #Sometimes mle difference is -1e-6, which gives NA...
  if ("tem" %in% mod) {
    nhpp_phi <- function(para){
      nhpp.qp.phi(par = para,
                  dat = y,
                  q = q,
                  N = N,
                  ny = ny,
                  V = V)
    }
     nhpp_dphi <- function(para){
       # sapply(y, function(yi) {
       #   numDeriv::grad(
       #     x = mle,
       #     func = function(par) {
       #       log(intensQ(
       #         par = par,
       #         dat = yi,
       #         q = q,
       #         N = N,
       #         ny = ny
       #       ))
       #     }
       #   ) / intensQ(
       #     par = mle,
       #     dat = yi,
       #     q = q,
       #     N = N,
       #     ny = ny
       #   )
       # }) %*% V
      nhpp.qp.dphi2(par = para,
                   dat = y,
                   q = q,
                   N = N,
                   ny = ny,
                   V = V)
    }
    phi.mle <- nhpp_phi(para = mle)
    # sapply(y, function(yi) {
    #   log(intensQ(
    #     par = mle,
    #     dat = yi,
    #     q = q,
    #     N = N,
    #     ny = ny)) /
    #     intensQ(
    #       par = mle,
    #       dat = yi,
    #       q = q,
    #       N = N,
    #       ny = ny)}) %*% V -
    #   phi.mle
    # nhpp.qp.dphi(
    #   par = mle,
    #   dat = y,
    #   q = q,
    #   N = N,
    #   ny = ny,
    #   V = V) -
    #   sapply(y, function(yi) {
    #     numDeriv::grad(
    #       x = mle,
    #       func = function(par) {
    #         log(intensQ(
    #           par = par,
    #           dat = yi,
    #           q = q,
    #           N = N,
    #           ny = ny
    #         ))
    #       }
    #     ) / intensQ(
    #       par = mle,
    #       dat = yi,
    #       q = q,
    #       N = N,
    #       ny = ny
    #     )
    #   }) %*% V
    q2num <- apply(pars, 1, function(par) {
      det(rbind(c(phi.mle - nhpp_phi(para = par)),
                nhpp_dphi(para = par)[-1, ]))
    })
    det_dphi_mle <- determinant(nhpp_dphi(para = mle))
    if (isTRUE(any(sign(q2num)/det_dphi_mle$sign * sign(r) < 0,
                   na.rm = TRUE))) {
      warning("Correction factor and likelihood root are of opposite sign - check output")
    }
    logq <- apply(pars, 1, function(par) {
      -0.5 * as.numeric(determinant(
        nhpp.qp.infomat(par = par,
                        dat = y,
                        u = threshold,
                        np = ny,
                        N = N,
                        q = q)[-1, -1])$modulus)
    }) + log(abs(q2num))
    qmlecontrib <- -as.numeric(det_dphi_mle$modulus) +
      0.5 * as.numeric(determinant(mle_opt$hessian)$modulus)
    logq <- logq + qmlecontrib
    qcor <- sign(q2num) * exp(logq)
    rstar <- ifelse(r == 0,
                    0,
                    r + (logq - log(abs(r)))/r)
    tem.max.opt <- function(psi) {
      para <- c(psi,
                constr.mle.Nquant(psi,
                                  start = pars[which.min(abs(psi - pars[,1]))[1],2:3]))
      pll <- -nhpp.nll(par = para,
                       u = threshold,
                       y = y,
                       ny = ny,
                       N = N,
                       q = q)
      rs <- 2 * pmax(0, maxll - pll)
      logq <- -0.5 * as.numeric(
        determinant(
        nhpp.qp.infomat(par = para,
                        dat = y,
                        u = threshold,
                        np = ny,
                        N = N,
                        q = q)[-1, -1])$modulus) +
        qmlecontrib +
        as.numeric(
          determinant(rbind(c(phi.mle - nhpp_phi(para = para)),
                          nhpp_dphi(para = para)[-1, ]))$modulus)
      rs + (logq - log(sqrt(abs(rs))))
    }
    # Technically, this equation should solves rstar=0
    #  so solving r*rstar=0 gives the second root if rstar != r,
    #   provided the objective function is zero there
    opt.tem <- optim(par = mle[1]+0.1,
                     fn = tem.max.opt,
                     method = "Brent",
                     lower = max(1e-05,
                                 mle[1] - std_error[1]),
                     upper = mle[1] + 1.5* std_error[1],
                     control = list(abstol = 1e-10))
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
    # The TEM suffers from a singularity at the MLE
    # see the thesis of Rongcai Li, UofT, 2001
    if (correction && length(psi) > 10) {
      ans <- spline.corr(ans)
    }
    if(abs(tem.rel) > 1e-5){
      ans$tem.psimax <- predict(smooth.spline(y = ans$psi, x = ans$rstar), 0)$y
      ans$tem.reliability <- 0
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

#' Sufficient directions for the truncated generalized Pareto
#'
#' @param y vector of observations
#' @param par vector of parameter, N-observation \code{q} quantile, scale and shape
#' @param s lower bound of truncation
#' @param N number of years for the maximum
#' @param u threshold
#' @param q probability; quantile level of target
#' @return an n by 3 matrix of sufficient directions
tgp.qp.Vfun <- function(par, y, s, N, u, q = 0.5){
  pos_max <- which.max(y)
  stopifnot(pos_max == length(y))
  y_stop <- y[pos_max]
  y_rest <- y[-pos_max]
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  cst1 <- sigma*((-N/log(q))^xi - 1)/xi
  cst2 <- (-N/log(q))^xi
  Vyl_qp <- -(s*xi - xi*y_stop)/((qp - s)*xi - sigma*cst2)
  Vyl_sigma <- (s*cst2 - y_stop*cst2)/((qp - s)*xi - sigma*cst2)
  Vyl_xi <- ((-(qp - s - cst1)*xi/sigma + 1)^(1/xi)*(1/(-(qp - s - cst1)*xi/sigma + 1)^(1/xi) - 1/(-(qp - y_stop - cst1)*xi/sigma + 1)^(1/xi))*(((sigma*cst2*log(-N/log(q))/xi - sigma*(cst2 - 1)/xi^2)*xi/sigma - (qp - s - cst1)/sigma)/(((qp - s - cst1)*xi/sigma - 1)*xi) + log(-(qp - s - cst1)*xi/sigma + 1)/xi^2) - (-(qp - s - cst1)*xi/sigma + 1)^(1/xi)*((((sigma*cst2*log(-N/log(q))/xi - sigma*(cst2 - 1)/xi^2)*xi/sigma - (qp - s - cst1)/sigma)/(((qp - s - cst1)*xi/sigma - 1)*xi) + log(-(qp - s - cst1)*xi/sigma + 1)/xi^2)/(-(qp - s - cst1)*xi/sigma + 1)^(1/xi) - (((sigma*cst2*log(-N/log(q))/xi - sigma*(cst2 - 1)/xi^2)*xi/sigma - (qp - y_stop - cst1)/sigma)/(((qp - y_stop - cst1)*xi/sigma - 1)*xi) + log(-(qp - y_stop - cst1)*xi/sigma + 1)/xi^2)/(-(qp - y_stop - cst1)*xi/sigma + 1)^(1/xi)))*sigma*(-(qp - y_stop - cst1)*xi/sigma + 1)^(1/xi + 1)/(-(qp - s - cst1)*xi/sigma + 1)^(1/xi)

  Vyr_qp <- ((((qp - s)*xi^2 - sigma*xi*cst2)*(-((qp - s)*xi - sigma*cst2)/sigma)^(1/xi) - ((qp - u)*xi^2 - sigma*xi*cst2)*(-((qp - u)*xi - sigma*cst2)/sigma)^(1/xi))*y_rest - ((qp - s)*u*xi^2 - sigma*u*xi*cst2)*(-((qp - s)*xi - sigma*cst2)/sigma)^(1/xi) - (s*sigma*xi*cst2 - (qp*s - s*u)*xi^2)*(-((qp - u)*xi - sigma*cst2)/sigma)^(1/xi) + ((s - u)*xi^2*y_rest - (qp*s - qp*u)*xi^2 + (s*sigma*cst2 - sigma*u*cst2)*xi)*(-(qp*xi - xi*y_rest - sigma*cst2)/sigma)^(1/xi))/(((qp^2 - qp*s - (qp - s)*u)*xi^2 + sigma^2*(-N/log(q))^(2*xi) + (sigma*u*cst2 - (2*qp*cst2 - s*cst2)*sigma)*xi)*(-((qp - s)*xi - sigma*cst2)/sigma)^(1/xi) - ((qp^2 - qp*s - (qp - s)*u)*xi^2 + sigma^2*(-N/log(q))^(2*xi) + (sigma*u*cst2 - (2*qp*cst2 - s*cst2)*sigma)*xi)*(-((qp - u)*xi - sigma*cst2)/sigma)^(1/xi))

  Vyr_sigma <- -((((qp*cst2 - s*cst2)*xi - sigma*(-N/log(q))^(2*xi))*(-((qp - s)*xi - sigma*cst2)/sigma)^(1/xi) - ((qp*cst2 - u*cst2)*xi - sigma*(-N/log(q))^(2*xi))*(-((qp - u)*xi - sigma*cst2)/sigma)^(1/xi))*y_rest - ((qp*cst2 - s*cst2)*u*xi - sigma*u*(-N/log(q))^(2*xi))*(-((qp - s)*xi - sigma*cst2)/sigma)^(1/xi) - (s*sigma*(-N/log(q))^(2*xi) - (qp*s*cst2 - s*u*cst2)*xi)*(-((qp - u)*xi - sigma*cst2)/sigma)^(1/xi) + ((s*cst2 - u*cst2)*xi*y_rest + s*sigma*(-N/log(q))^(2*xi) - sigma*u*(-N/log(q))^(2*xi) - (qp*s*cst2 - qp*u*cst2)*xi)*(-(qp*xi - xi*y_rest - sigma*cst2)/sigma)^(1/xi))/
    (((qp^2 - qp*s - (qp - s)*u)*xi^2 + sigma^2*(-N/log(q))^(2*xi) + (sigma*u*cst2 - (2*qp*cst2 - s*cst2)*sigma)*xi)*(-((qp - s)*xi - sigma*cst2)/sigma)^(1/xi) - ((qp^2 - qp*s - (qp - s)*u)*xi^2 + sigma^2*(-N/log(q))^(2*xi) + (sigma*u*cst2 - (2*qp*cst2 - s*cst2)*sigma)*xi)*(-((qp - u)*xi - sigma*cst2)/sigma)^(1/xi))

  Vyr_xi <- -sigma*(-(qp - y_rest - cst1)*xi/sigma + 1)^(1/xi + 1)*(((((sigma*cst2*log(-N/log(q))/xi - sigma*(cst2 - 1)/xi^2)*xi/sigma - (qp - u - cst1)/sigma)/(((qp - u - cst1)*xi/sigma - 1)*xi) + log(-(qp - u - cst1)*xi/sigma + 1)/xi^2)/(-(qp - u - cst1)*xi/sigma + 1)^(1/xi) - (((sigma*cst2*log(-N/log(q))/xi - sigma*(cst2 - 1)/xi^2)*xi/sigma - (qp - y_rest - cst1)/sigma)/(((qp - y_rest - cst1)*xi/sigma - 1)*xi) + log(-(qp - y_rest - cst1)*xi/sigma + 1)/xi^2)/(-(qp - y_rest - cst1)*xi/sigma + 1)^(1/xi))/(1/(-(qp - s - cst1)*xi/sigma + 1)^(1/xi) - 1/(-(qp - u - cst1)*xi/sigma + 1)^(1/xi)) -
                                                                      ((((sigma*cst2*log(-N/log(q))/xi - sigma*(cst2 - 1)/xi^2)*xi/sigma - (qp - s - cst1)/sigma)/(((qp - s - cst1)*xi/sigma - 1)*xi) + log(-(qp - s - cst1)*xi/sigma + 1)/xi^2)/(-(qp - s - cst1)*xi/sigma + 1)^(1/xi) - (((sigma*cst2*log(-N/log(q))/xi - sigma*(cst2 - 1)/xi^2)*xi/sigma - (qp - u - cst1)/sigma)/(((qp - u - cst1)*xi/sigma - 1)*xi) + log(-(qp - u - cst1)*xi/sigma + 1)/xi^2)/(-(qp - u - cst1)*xi/sigma + 1)^(1/xi))*(1/(-(qp - u - cst1)*xi/sigma + 1)^(1/xi) - 1/(-(qp - y_rest - cst1)*xi/sigma + 1)^(1/xi))/(1/(-(qp - s - cst1)*xi/sigma + 1)^(1/xi) - 1/(-(qp - u - cst1)*xi/sigma + 1)^(1/xi))^2)*(1/(-(qp - s - cst1)*xi/sigma + 1)^(1/xi) - 1/(-(qp - u - cst1)*xi/sigma + 1)^(1/xi))

  rbind(
    cbind(Vyr_qp, Vyr_sigma, Vyr_xi),
    cbind(Vyl_qp, Vyl_sigma, Vyl_xi)
  )

}

#' Canonical parameter for the truncated generalized Pareto
#' @inheritParams tgp.qp.Vfun
#' @param V matrix of sufficient direction
#' @return a 3 vector of canonical parameters
tgp.qp.phi <- function(par, y, s, N, u, q = 0.5, V){
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  psi_c <- (xi + 1)/(qp*xi - xi*y - sigma*(-N/log(q))^xi)
  crossprod(psi_c, V)
}

#' Derivative of the canonical parameter
#' @inheritParams tgp.qp.Vfun
#' @param V matrix of sufficient direction
#' @return a 3 by 3 matrix
tgp.qp.dphi <- function(par, y, s, N, ny, u, q = 0.5, V){
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  psi_qp <- -(xi + 1)*xi/(qp*xi - xi*y - sigma*(-N/log(q))^xi)^2
  psi_sigma <- (xi + 1)*(-N/log(q))^xi/(qp*xi - xi*y - sigma*(-N/log(q))^xi)^2
  psi_xi <- (sigma*(-N/log(q))^xi*log(-N/log(q)) - qp + y)*(xi + 1)/(qp*xi - xi*y - sigma*(-N/log(q))^xi)^2 + 1/(qp*xi - xi*y - sigma*(-N/log(q))^xi)
  crossprod(cbind(psi_qp, psi_sigma, psi_xi), V)
}

#' Canonical parameter of the inhomogeneous point process (quantile)
#'
#' Monte Carlo approximation for the canonical parameter parametrized
#' in terms of \code{q} quantile of the \code{N} observation maximum
#' This is a discrete approximation that leads to second order accuracy
#' with the tangent exponential model approximation


#' Directional derivative of canonical parameter
#'
#' This term is \code{V} times the derivative
#' of the Monte Carlo approximation to the canonical parameter
#' @inheritParams nhpp.qp.Vfun
#' @param V an \code{n} by \code{p} matrix whose columns span the sufficient directions
nhpp.qp.dphi2 <- function(par, dat, q, N, ny, V){
  y <- dat
  qp <- par[1]
  sigma <- par[2]
  xi <- par[3]
  cst <- sigma*((-N/log(q))^xi - 1)/xi
  dphi_qp <- -(-(qp - y - cst)*xi/sigma + 1)^(1/xi)*(xi + 1)*log(ny*(-(qp - y - cst)*xi/sigma + 1)^(-1/xi - 1)/sigma)/ny + (-(qp - y - cst)*xi/sigma + 1)^(1/xi)*(1+xi)/ny

  dphi_sigma <- ((sigma*(-N/log(q))^xi + qp - y)*(-(qp*xi - xi*y - sigma*(-N/log(q))^xi)/sigma)^(1/xi)*log(-ny/((qp*xi - xi*y - sigma*(-N/log(q))^xi)*(-(qp*xi - xi*y - sigma*(-N/log(q))^xi)/sigma)^(1/xi))) -
(sigma*(-N/log(q))^xi + qp - y)*(-(qp*xi - xi*y - sigma*(-N/log(q))^xi)/sigma)^(1/xi))/(ny*sigma)

  dphi_xi <- -sigma*(-(qp - y - cst)*xi/sigma + 1)^(1/xi + 1)*(((sigma*(-N/log(q))^xi*log(-N/log(q))/xi - cst/xi)*xi/sigma -
 (qp - y - cst)/sigma)*(1/xi + 1)/((qp - y - cst)*xi/sigma - 1) +
  log(-(qp - y - cst)*xi/sigma + 1)/xi^2)*log(ny*(-(qp - y - cst)*xi/sigma + 1)^(-1/xi - 1)/sigma)/ny + sigma*(-(qp - y - cst)*xi/sigma + 1)^(1/xi + 1)*(((sigma*(-N/log(q))^xi*log(-N/log(q))/xi - cst/xi)*xi/sigma -
(qp - y - cst)/sigma)*(1/xi + 1)/((qp - y - cst)*xi/sigma - 1) + log(-(qp - y - cst)*xi/sigma + 1)/xi^2)/ny

  rbind(dphi_qp, dphi_sigma, dphi_xi) %*% V
}

# Nonhomogeneous Poisson process likelihood and HOA approximation
# Parametrized in terms of pth quantile of the T-year maximum

#' Sufficient directions for nonhomogeneous Poisson process
#'
#' Matrix of sufficient directions for the nonhomogeneous
#' Poisson process likelihood parametrized in terms of the
#' p quantile of the T observation maximum
#' @param par vector of parameters estimates, with quantiles of the T\code{ny}-observation, scale and shape parameters
#' @param dat vector of threshold exceedances
#' @param q probability level of the quantile
#' @param N number of observations for the maximum
# nhpp.qp.Vfun <- function(par, dat, q, N, ny){
#   y <- dat
#   stopifnot(length(par) == 3L, par[1:2] > 0)
#   qp <- par[1]
#   sigma <- par[2]
#   xi <- par[3]
#   V_qp <- ny*(-(qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*(1+xi)/sigma^2
#   V_sigma <- -ny*(-(qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 2)*((qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma^2 + ((-N/log(q))^xi - 1)/sigma)*(1/xi + 1)/sigma - ny*(-(qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma^2
#   V_xi <- ny*(-(qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)*(((sigma*(-N/log(q))^xi*log(-N/log(q))/xi - sigma*((-N/log(q))^xi - 1)/xi^2)*xi/sigma - (qp - y - sigma*((-N/log(q))^xi - 1)/xi)/sigma)*(1/xi + 1)/((qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma - 1) + log(-(qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)/xi^2)/sigma
#   cbind(V_qp, V_sigma, V_xi)
# }


#' Directional derivative of the log likelihood
#'
#' This term is \code{V} times the sample
#' space derivative of the log-likelihood
#'
#' @inheritParams nhpp.qp.Vfun
#' @inheritParams nhpp.qp.dphi
# nhpp.qp.phi <- function(par, dat, q, N, ny, V){
#   y <- dat
#   qp <- par[1]
#   sigma <- par[2]
#   xi <- par[3]
#   phi_a <- sigma*log(ny*(-(qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1)/sigma)/(ny*(-(qp - y - sigma*((-N/log(q))^xi - 1)/xi)*xi/sigma + 1)^(-1/xi - 1))
#   rbind(phi_a) %*% V
# }


#' Standard log likelihood for the Maiquetia data
#'
#' We consider a nonhomogeneous Poisson process for threshold exceedances \code{y} above \code{u},
#' with intensity function \deqn{\Lambda(y) = n_y(1+\xi(y-\mu)/\sigma)^(-1/\xi)}
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
#' @param q quantile level
maiq.std.nll <- function(par, u, y, z, ny, N, q){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  mu <- qp + sigma/xi*(1-(-N/log(q))^xi)
  - suppressWarnings(
    ifelse(length(y) == 0, 0,
           sum(mev::pp.ll(par = c(mu, sigma, xi),
               dat = y, u = u, np = ny))) -
    ifelse(length(z) == 0, 0,
           sum(mev::gev.ll(par = c(mu, sigma, xi), dat = z))))
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
maiq.fc.nll <- function(par, u, s, y, z, ny, N, q){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  mu <- qp + sigma/xi*(1-(-N/log(q))^xi)
  Lambda_f <- function(x){
    ny*pmax(1+xi*(x-mu)/sigma,0)^(-1/xi)}
  Lambda_s <- Lambda_f(s)
  Lambda_u <- Lambda_f(u)
  if(any(isTRUE(all.equal(Lambda_s,0, check.attributes = FALSE)),
         isTRUE(all.equal(Lambda_u,0, check.attributes = FALSE)))){
    return(1e10)
  }
  ll <- -length(y)*log(sigma) - (1+1/xi)*sum(log(pmax(1+xi*(y-mu)/sigma, 0))) -
    (length(y)-1)*log(Lambda_u-Lambda_s) - log(Lambda_s) +
    ifelse(length(z) == 0L, 0,
           sum(mev::gev.ll(par = c(mu, sigma, xi), dat = z)))
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
#' @param q probability level for the quantile
#' @param N number of block for the maximum
gevq.Vfun <- function(par, z, N, q){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  Nr <- -N/log(q)
cbind(-1,
      (qp - z)/sigma,
      ((sigma*(Nr)^xi*log(Nr) - qp)*xi + xi*z + (qp*xi - xi*z - sigma*(Nr)^xi)*log(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma))/xi^2
)
}

gevq.phi <- function(par, z, N, q, V){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
 phi_a <- as.numeric(((xi + 1)*(-(qp*xi - xi*z - sigma*(-N/log(q))^xi)/sigma)^(1/xi) - 1)/((qp*xi - xi*z - sigma*(-N/log(q))^xi)*(-(qp*xi - xi*z - sigma*(-N/log(q))^xi)/sigma)^(1/xi)))
 rbind(phi_a) %*% V
}

gevq.dphi <- function(par, z, N, q, V){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  Nr <- -N/log(q)
  dphi_qp <- -((xi^2 + xi)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi) - xi - 1)/((qp^2*xi^2 + xi^2*z^2 - 2*qp*sigma*xi*(Nr)^xi + sigma^2*(Nr)^(2*xi) - 2*(qp*xi^2 - sigma*xi*(Nr)^xi)*z)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi))
  dphi_sigma <- ((sigma*xi*(Nr)^xi + sigma*(Nr)^xi)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi) - sigma*(Nr)^xi - qp + z)/((qp^2*sigma*xi^2 + sigma*xi^2*z^2 - 2*qp*sigma^2*xi*(Nr)^xi + sigma^3*(Nr)^(2*xi) - 2*(qp*sigma*xi^2 - sigma^2*xi*(Nr)^xi)*z)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi))
  dphi_xi <- -((sigma*(Nr)^xi*log(Nr) - qp)*xi^2 + (sigma*(Nr)^xi*log(Nr) - qp)*xi + (xi^2 + xi)*z - (sigma*xi^3*(Nr)^xi*log(Nr) + (sigma*(Nr)^xi*(log(Nr) - 1) - qp)*xi^2 + xi^2*z)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi) + (qp*xi - xi*z - sigma*(Nr)^xi)*log(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma))/((qp^2*xi^4 + xi^4*z^2 - 2*qp*sigma*xi^3*(Nr)^xi + sigma^2*xi^2*(Nr)^(2*xi) - 2*(qp*xi^4 - sigma*xi^3*(Nr)^xi)*z)*(-(qp*xi - xi*z - sigma*(Nr)^xi)/sigma)^(1/xi))
  rbind(dphi_qp, dphi_sigma, dphi_xi) %*% V
}

#' Inequality constraints for optimization
#'
#' This function is used for the Maiquetia
#' full conditional likelihood with stopping rule
#' @inheritParams llmaiq_std
#' @inheritParams llmaiq.fc
#' @keywords internal
constraints <- function(par, u, s, y, z, ny, N, q){
  qp <- par[1]; sigma <- par[2]; xi <- par[3]
  mu <- qp + sigma/xi*(1-(-N/log(q))^xi)
  xv <- ifelse(xi > 0, min(c(u,s,y,z)), max(c(u,s,y,z)))
  as.numeric(c(sigma, 1 - xi, xi + 1, sigma + xi*(xv-mu)))
}


maiq.fc.Vfun <- function(par,
                         u,
                         s,
                         y,
                         z,
                         ny,
                         N,
                         q){
  V1 <- tgp.qp.Vfun(par = par,
                    y = y,
                    u = u,
                    s = s,
                    q = q,
                    N = N)
  if(!missing(z) && length(z) >= 1L){
    V2 <- gevq.Vfun(par = par,
                    z = z,
                    N = N,
                    q = q)
    return(rbind(V1, V2))
  } else{
    return(V1)
  }
}

maiq.fc.phi <- function(par,
                        u,
                        s,
                        y,
                        z,
                        ny,
                        N,
                        q,
                        V){
  tgp.qp.phi(par = par,
             y = y,
             s = s,
             u = u,
             N = N,
             q = q,
             V = V[(1:length(y)),]) +
  gevq.phi(par = par,
           z = z,
           N = N,
           q = q,
           V = V[-(1:length(y)),,drop = FALSE])
}

maiq.fc.dphi <- function(par,
                         u,
                         s,
                         y,
                         z,
                         ny,
                         N,
                         q,
                         V){
  tgp.qp.dphi(par = par,
             y = y,
             s = s,
             u = u,
             N = N,
             q = q,
             ny = ny,
             V = V[(1:length(y)),]) +
    gevq.dphi(par = par,
             z = z,
             N = N,
             q = q,
             V = V[-(1:length(y)),])
}

#' Profile log likelihood of the Maiquetia example
#' with full conditioning
#'
#' The scalar of interest is the \code{q} quantile of the
#' the \code{N}-block maximum, assuming there are \code{ny}
#' blocks.
#'
#' @param s stopping rule
#' @param z vector of block maximum
#' @param y ordered exceedances; the last observation triggers the stopping rules
#' @param psi grid over which to evaluate the profile
#' @param start vector of maximum likelihood or good starting value
#' @param mod string; either \code{profile} for the profile likelihood or \code{tem} for the tangent exponential model approximation.
#' @param N numeric of length one; the number of block maximum for which the quantile is sought
#' @return an object of class \code{eprof}
maiq.fc.profile <- function(psi,
                           start,
                           mod = c("profile","tem"),
                           y,
                           z,
                           s,
                           ny,
                           N = 50,
                           q = 0.5,
                           correction = TRUE,
                           threshold = 0,
                           plot = TRUE, ...){
  names_mle <- c("Nquant","scale","shape")
  param <- "Nquant"
  mod <- match.arg(arg = mod,
                   choices = c("profile", "tem"),
                   several.ok = TRUE)
  ymin <- min(y)
  ymax <- max(y)
  zmin <- min(z)
  zmax <- max(z)
  #' Wrapper for inequality constraints
  #'
  #' Wrapper around \code{constraints} function
  #' @keywords internal
  constraints.fn <- function(par){
    constraints(par = par,
                z = z,
                s = s,
                y = y,
                u = threshold,
                ny = ny,
                N = N,
                q = q)}

  #' Full conditional log likelihood for Maiquetia data
  #'
  #' Wrapper around \code{maiq.fc.nll} for \code{alabama::auglag},
  #' which complains about parameters in \code{hin} or
  #' \code{fn} without the same arguments.
  #'
  #' @keywords internal
  maiq.fc.nll.fn <- function(par){
    maiq.fc.nll(par = par,
                z = z,
                s = s,
                y = y,
                u = threshold,
                ny = ny,
                N = N,
                q = q)
  }
  mle_opt <- Rsolnp::solnp(pars = start,
                           fun = maiq.fc.nll.fn,
                           ineqfun = constraints.fn,
                           ineqLB = rep(0,4),
                           ineqUB = rep(Inf, 4),
                           control = list(trace = 0,
                                          tol = 1e-12))
  mle <- mle_opt$par
  names(mle) <- names_mle
  mle_opt$hessian <-
    numDeriv::hessian(
      fun = maiq.fc.nll.fn,
      x = mle,
      method.args=list(eps=1e-4,
                       d=1e-5,
                       zero.tol=sqrt(.Machine$double.eps/7e-7),
                       r=4,
                       v=2,
                       show.details=FALSE))
  std_error <- sqrt(diag(solve(mle_opt$hessian)))
  # Extract the components, notably V for model `tem`. Keep other components for optimization

  if("tem" %in% mod){
    V <- maiq.fc.Vfun(par = mle,
                       y = y,
                       s = s,
                       N = N,
                       ny = ny,
                       z = z,
                       u = threshold,
                       q = q)
  }
  maxll <- -mle_opt$value[length(mle_opt$value)]
  constr.mle.Nquant <- function(psi, start = mle[-1]) {
      suppressWarnings(Rsolnp::solnp(
        pars = start,
        fun = function(lambda){
          maiq.fc.nll(par = c(psi, lambda),
                      y = y,
                      z = z,
                      s = s,
                      u = threshold,
                      ny = ny,
                      N = N,
                      q = q)},
        ineqfun = function(lambda){
          constraints.fn(c(psi, lambda))},
        ineqUB = rep(Inf, 4),
        ineqLB = rep(0, 4),
        control = list(trace=0))$pars)
    }
  # Missing psi vector
  if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
    #compute profile log-lik on a grid left and right of the MLE
    psi <- unique(seq(from = mle[1] - 2 * std_error,
                      to = mle[1] + 3 * std_error,
                      length.out = 101L))
    psi <- psi[psi > 0]
  }
  pars <- matrix(nrow = length(psi), ncol = 3)
  pars[,1] <- psi
  mid <- which.min(sapply(psi, function(p){abs(p-mle[1])}))
  pars[mid, 2:3] <- constr.mle.Nquant(psi = psi[mid])
  for(i in (mid-1):1){
    pars[i,2:3] <- constr.mle.Nquant(psi = psi[i],
                                     start = pars[i+1,2:3])
  }
  for(i in (mid+1):length(psi)){
    pars[i,2:3] <- constr.mle.Nquant(psi = psi[i],
                                     start = pars[i-1,2:3])
  }
  # Profile log likelihood values for psi
  profll <- apply(pars, 1, function(par) {
    -maiq.fc.nll(par = par,
                 y = y,
                 z = z,
                 s = s,
                 u = threshold,
                 ny = ny,
                 N = N,
                 q = q)
  })
  r <- sign(mle[param] - psi) *
    sqrt(2 * pmax(0, maxll - profll))
  #Sometimes mle difference is -1e-6, which gives NA...
  if ("tem" %in% mod) {
    phi_par <- function(par){
      maiq.fc.phi(par = par,
                  y = y,
                  z = z,
                  N = N,
                  u = threshold,
                  ny = ny,
                  V = V,
                  q = q,
                  s = s)}
    dphi_par <- function(par){
      maiq.fc.dphi(par = par,
                   y = y,
                   z = z,
                   N = N,
                   u = threshold,
                   ny = ny,
                   V = V,
                   q = q,
                   s = s)}
    phi_mle <- phi_par(par = mle)
    dphi_mle <- dphi_par(par = mle)
    det_dphi_mle <- det(dphi_mle)
    q2num <- apply(pars, 1, function(par) {
      det(rbind(c(phi_mle - phi_par(par)),
                dphi_par(par)[-1, ]))
    })
    if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
      warning("Correction factor and likelihood root are of opposite sign - check output")
    }
    # Obtain information matrix through numerical differentiation
    maiq.infomat <- function(par){ #negative of hessian of ll
      numDeriv::hessian(func = function(x){
        #loglik already negated
        maiq.fc.nll(par = x,
                    u = threshold,
                    s = s,
                    y = y,
                    z = z,
                    ny = ny,
                    N = N,
                    q = q)
      },
        x = par,
      method.args = list(eps = 1e-4,
                         d = 1e-5,
                         zero.tol = sqrt(.Machine$double.eps/7e-7),
                         r = 4,
                         v = 2,
                         show.details = FALSE))}
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
      para <- c(psi,
                constr.mle.Nquant(psi, start = pars[which.min(abs(psi - pars[,1]))[1],2:3])) #better starting values?
      pll <- -maiq.fc.nll(par = para, z = z, y = y, ny = ny, u = threshold, s = s, q = q, N = N)
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
      ans <- mev::spline.corr(ans)
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
maiq.std.profile <- function(psi, start, mod = c("profile","tem"),
                                y, z, ny, N = 50, q = 0.5,
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
  constraints.fn <- function(par){
    qp <- par[1]; sigma <- par[2]; xi <- par[3]
    mu <- qp + sigma/xi*(1-(-N/log(q))^xi)
    xv <- ifelse(xi > 0, min(c(y, z, threshold)), max(c(threshold, y, z)))
    as.numeric(c(sigma/1000, 1 - xi, xi + 1, sigma + xi*(xv-mu)))
  }

  #' Full conditional log likelihood for Maiquetia data
  #'
  #' Wrapper around \code{maiq.std.nll} for \code{alabama::auglag},
  #' which complains about parameters in \code{hin} or
  #' \code{fn} without the same arguments.
  #'
  #' @keywords internal
  maiq_std_nll_fn <- function(par){
    maiq.std.nll(par = par, z = z, y = y,
                u = threshold, ny = ny, N = N, q = q)
  }

  mle_opt <- Rsolnp::solnp(pars = start,
                           fun = maiq_std_nll_fn,
                           ineqfun = constraints.fn,
                           ineqLB = rep(0, 4),
                           ineqUB = rep(Inf, 4),
                           control = list(trace = 0))
  mle <- mle_opt$pars
  names(mle) <- names_mle
  mle_opt$hessian <-
    numDeriv::hessian(
      fun = maiq_std_nll_fn,
      x = mle,
      method.args = list(
        eps = 1e-4,
        d = 1e-5,
        zero.tol = sqrt(.Machine$double.eps/7e-7),
        r = 4,
        v = 2,
        show.details = FALSE))
  std_error <- sqrt(diag(solve(mle_opt$hessian)))
  # Extract the components, notably V for model `tem`. Keep other components for optimization
  if("tem" %in% mod){
     V1 <- nhpp.qp.Vfun(par = mle, dat = y, N = N, q = q, ny = ny)
     V2 <- gevq.Vfun(par = mle, z = z, N = N, q = q)
     V <- rbind(V1, V2)
  } else{
    V <- NULL
  }
  maxll <- -mle_opt$value[length(mle_opt$value)]

  constr.mle.Nquant <- function(psi, start = mle[-1]) {
    suppressWarnings(
      Rsolnp::solnp(
        pars = start,
        fun = function(lambda) {
           maiq.std.nll(par = c(psi, lambda),
                        y = y,
                        z = z,
                        u = threshold,
                        ny = ny,
                        N = N,
                        q = q) },
         ineqfun = function(lambda){
           constraints.fn(c(psi, lambda))},
        ineqLB = rep(0,4),
        ineqUB = rep(Inf, 4),
        control = list(trace=0))$pars
    )
  }
  # Missing psi vector
  if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
    #compute profile log-lik on a grid left and right of the MLE
    psi <- unique(seq(from = mle[1]-2.5*std_error,
                      to = mle[1]+3*std_error,
                      length.out = 101L))
    psi <- psi[psi > 0]
  }
  pars <- matrix(nrow = length(psi), ncol = 3)
  pars[,1] <- psi
  mid <- which.min(sapply(psi, function(p){abs(p-mle[1])}))
  pars[mid, 2:3] <- constr.mle.Nquant(psi = psi[mid])
  for(i in (mid-1):1){
    pars[i,2:3] <-
      constr.mle.Nquant(psi = psi[i],
                        start = pars[i+1,2:3])
  }
  for(i in mid:length(psi)){
    pars[i,2:3] <-
      constr.mle.Nquant(psi = psi[i],
                        start = pars[i-1,2:3])
  }
   # Profile log likelihood values for psi
  profll <- apply(pars, 1, function(par) {
    -maiq.std.nll(par = par, y = y, z = z, u = threshold, ny = ny, N = N, q = q)
  })
  r <- sign(mle[param] - psi) * sqrt(2 * pmax(0, maxll - profll))
  #Sometimes mle difference is ~ -1e-6, which gives NA...
  if ("tem" %in% mod) {
    maiq.std.phi <- function(para){
      nhpp.qp.phi(par = para,
                  dat = y,
                  q = q,
                  N = N,
                  ny = ny,
                  V = V1) +
        gevq.phi(par = para,
                 z = z,
                 N = N,
                 q = q,
                 V = V2)
    }
    maiq.std.dphi <- function(para){
      nhpp.qp.dphi2(par = para,
                    dat = y,
                    q = q,
                    N = N,
                    ny = ny,
                    V = V1) +
        gevq.dphi(par = para,
                  z = z,
                  N = N,
                  q = q,
                  V = V2)
    }
    phi.mle <- maiq.std.phi(para = mle)
    q2num <- apply(pars, 1, function(par) {
      det(rbind(c(phi.mle - maiq.std.phi(para = par)),
                maiq.std.dphi(para = par)[-1, ]))
    })
    det_dphi_mle <- determinant(maiq.std.dphi(para = mle))
    if (isTRUE(any(sign(q2num)/det_dphi_mle$sign * sign(r) < 0, na.rm = TRUE))) {
      warning("Correction factor and likelihood root are of opposite sign - check output")
    }
    # Obtain information matrix through numerical differentiation
    maiq.std.infomat <- function(par){ #negative of hessian of ll
      numDeriv::hessian(func = function(x){ #loglik already negated
        maiq.std.nll(par = x,
                     u = threshold,
                     y = y,
                     z = z,
                     ny = ny,
                     N = N,
                     q = q)
      },
      x = par,
      method.args = list(eps = 1e-4,
                         d = 1e-5,
                         zero.tol = sqrt(.Machine$double.eps/7e-7),
                         r = 4,
                         v = 2,
                         show.details = FALSE))}

    logq <- apply(pars, 1, function(par) {
      -0.5 * log(det(maiq.std.infomat(par = par)[-1, -1]))
    }) + log(abs(q2num))
    qmlecontrib <- - as.numeric(det_dphi_mle$modulus) + 0.5 *
      as.numeric(determinant(maiq.std.infomat(par = mle))$modulus)
    logq <- logq + qmlecontrib
    logq[sign(q2num)/det_dphi_mle$sign != sign(r)] <- NA
    qcor <- sign(q2num) * exp(logq)
    rstar <- ifelse(r == 0,
                    0,
                    r + (logq - log(abs(r)))/r)
    tem.max.opt <- function(psi) {
      para <- c(psi,
                constr.mle.Nquant(psi,
                                  start = pars[which.min(abs(psi - pars[,1]))[1],2:3]))
      pll <- -maiq.std.nll(par = para,
                           u = threshold,
                           y = y,
                           z = z,
                           ny = ny,
                           N = N,
                           q = q)
      rs <- 2 * pmax(0, maxll - pll)
      logq <- -0.5 * as.numeric(determinant(
        maiq.std.infomat(par = para)[-1,-1])$modulus) +
        qmlecontrib +
        as.numeric(determinant(rbind(
          c(phi.mle - maiq.std.phi(para = para)),
          maiq.std.dphi(para = para)[-1, ]))$modulus)
      rs + (logq - log(sqrt(abs(rs))))
    }
    # Technically, this equation should solve rstar=0, so solving r*rstar=0 gives the second root if rstar != r
    # #Bespoke implementation: check only to the right of the MLE
    opt.tem <- optim(par = mle[1]+10,
                     fn = tem.max.opt,
                     method = "Brent",
                     lower = max(1e-05,
                                 mle[1] - std_error[1]),
                     upper = mle[1] + 1.5* std_error[1],
                     control = list(abstol = 1e-10))
    tem.max <- opt.tem$par
    tem.rel <- opt.tem$value
  }
  # Return profile likelihood and quantities of interest (modified likelihoods)
  colnames(pars) <- names(mle)
  ans <- list(mle = mle,
              pars = pars,
              psi.max = as.vector(mle[param]),
              param = param,
              std_error = std_error,
              psi = psi,
              pll = profll,
              maxpll = maxll,
              r = r,
              V = V)
  if ("tem" %in% mod) {
    ans$q <- qcor
    ans$rstar <- rstar
    ans$tem.psimax <- as.vector(tem.max)
    ans$normal <- c(ans$psi.max, ans$std_error)
    ans$tem.reliability <- tem.rel
    if (correction && length(psi) > 10) {
      ans <- mev::spline.corr(ans)
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
