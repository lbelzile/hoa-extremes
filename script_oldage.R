#' Likelihood for left-truncated and right-censored generalized Pareto variates
#'
#' Computes the log-likelihood for generalized Pareto or exponential observations.
#'
#' @param par vector of scale and shape
#' @param dat vector of threshold exceedances
#' @param rightcens logical indicating right-censoring (\code{TRUE} for censored)
#' @param slow lower truncation limit
#' @param expo logical; should an exponential model be fitted? Default to \code{FALSE}
#' @return log-likelihood value
gpd_cens <- function(par, dat, rightcens, slow, expo = FALSE){
  if(expo){
    stopifnot(length(par) == 1L)
    shape <- 0
  } else{
    shape <- par[2]
  }
  if(par[1] < 0 || shape < (-1+1e-8)){
    return(1e10)
  }
  g1 <- intersect(which(!rightcens), which(slow > 0))
  g2 <- intersect(which(rightcens), which(slow > 0))
  ll <- 0
  #Contribution from observations in sampling frame (density divided by truncation interval)
  if(sum(!rightcens)>0){
    ll <- sum(evd::dgpd(loc = 0, scale = par[1], shape = shape, x = dat[!rightcens], log = TRUE))
    if(length(g1) > 0){
      ll <- ll - sum(log(1-evd::pgpd(slow[g1], loc = 0, scale = par[1], shape = shape)))  #right censored individuals
    }
  }
  if(sum(rightcens)>0){
    ll <- ll +  sum(log(1-evd::pgpd(dat[rightcens], loc = 0, scale = par[1], shape = shape)))
    if(length(g2) > 0){
      ll <- ll - sum(log(1-evd::pgpd(slow[g2], loc = 0, scale = par[1], shape = shape)))
    }
  }
  if (!is.finite(ll)) {
    return(1e10)
  }  else {
    return(-ll)
  }
}


#' Observed information matrix
#'
#' Observed information matrix for left-truncated right-censored generalized Pareto observations
#' @inheritParams gpd_cens
#' @param theta vector of endpoing and shape parameters for the generalized Pareto
#' @return a two by two matrix
obs.infomat <- function(theta, dat, rightcens, slow){
  endpt = theta[1]; xi = theta[2];
  if(endpt < 0){
    return(matrix(rep(Inf, 4), nrow = 2))
  }
  j11 <- (endpt^2*xi - dat^2 + 2*dat*endpt)/((dat^2*endpt^2 - 2*dat*endpt^3 + endpt^4)*xi)
  j12 <- -dat/((dat*endpt - endpt^2)*xi^2)
  j22 <- (xi - 2*log(-(dat - endpt)/endpt))/xi^3

  jc11 <- -(dat^2 - 2*dat*endpt)/((dat^2*endpt^2 - 2*dat*endpt^3 + endpt^4)*xi)
  jc12 <- -dat/((dat*endpt - endpt^2)*xi^2)
  jc22 <- -2*log(-(dat - endpt)/endpt)/xi^3

  jt11 <- (slow^2 - 2*slow*endpt)/((slow^2*endpt^2 - 2*slow*endpt^3 + endpt^4)*xi)
  jt12 <- slow/((slow*endpt - endpt^2)*xi^2)
  jt22 <- 2*log(-(slow - endpt)/endpt)/xi^3

  info11 <- sum(j11[!rightcens])+sum(jc11[rightcens])+sum(jt11)
  info12 <- sum(j12[!rightcens])+sum(jc12[rightcens])+sum(jt12)
  info22 <- sum(j22[!rightcens])+sum(jc22[rightcens])+sum(jt22)
  -cbind(c(info11, info12), c(info12, info22))
}

#' Profile likelihood for endpoint
#'
#' Profile likelihood for the endpoint of the generalized Pareto distribution
#' for left-truncated and right-censored observations.
#'
#' @inheritsParam gpd_cens
#' @param psi value of the endpoint at which to compute the profile
#' @return a vector of length three containing twice the negative log-likelihood value, the endpoint value and the maximum of the nuisance lambda (i.e., the shape parameter).
prof_gpd_cens <- function(psi, dat, rightcens, slow){
  opt <- optim(fn = function(lambda){
    gpd_cens(par = c(-lambda*psi, lambda),
             dat = dat,
             rightcens = rightcens,
             slow = slow)},
               method = "Brent",
    par = -0.01,
    lower = -0.98,
    upper = -1e-7, control = list(reltol=1e-12))
  res <- -2*opt$value
  attributes(res) <- list("param" = c(psi, opt$par))
  #res
  return(c(res, psi, opt$par))
}
loglik_endpt <- function(theta, dat, rightcens, slow){
  gpd_cens(par = c(-theta[2]*theta[1], theta[2]), dat = dat, rightcens = rightcens, slow = slow)
}

#' Likelihood root function
#'
#' This function returns the likelihood root of the profile log-likelihood for the endpoint of the generalized Pareto distribution with left-truncated and right-censored data.
#' Specifically, \eqn{-r^2/2} is the profile likelihood and the two-sided p-value is\code{qchisq(p, 1)/2}.
#'
#' @param psi parameter of the endpoint at which to compute the p-value
#' @param thetahat maximum likelihood estimates of the endpoint and the shape parameters
#' @inheritParams gpd_cens
#' @return a p-value
rfun_endpt <- function(psi, thetahat, dat, rightcens, slow){
  # Log-likelihood function in terms of parameters (endpt, xi)
  llp <- prof_gpd_cens(psi = psi, dat = dat, rightcens = rightcens, slow = slow)
  sign(thetahat[1]-llp[2])*sqrt(-2*loglik_endpt(thetahat, dat = dat, rightcens = rightcens, slow = slow) - llp[1])
}

modif_prof_endpt_empcov <- function(theta, thetahat, dat, rightcens, slow){

  empcov <- function(par, dat, rightcens, slow){
    zeta = par[1]; xi = par[2]
    ifelse(rightcens,
           log(-dat/zeta + 1)/xi^2,
           -1/xi + log(-dat/zeta + 1)/xi^2)- log(-slow/zeta + 1)/xi^2

  }
  # Log-likelihood function in terms of parameters (endpt, xi)
      - gpd_cens(par = c(-theta[2]*theta[1], theta[2]),
             dat = dat,
             rightcens = rightcens,
             slow = slow) +
      0.5*log(obs.infomat(theta = theta,
                          dat = dat,
                          rightcens = rightcens,
                          slow = slow)[2,2]) -
      log(crossprod(empcov(theta, dat, rightcens, slow),
                    empcov(thetahat, dat, rightcens, slow)))
}


#' Confidence intervals based on (modified) profile likelihoods
#'
#' This code is adapted from the mev package (mev:::confint.eprof)
#' @param object a list containing informations about the profile likelihood in the same format as the \code{hoa} package
#' @param parm string or numerical vector giving the type of interval to consider
#' @param level probability level of the confidence interval
#' @param prob vector of length 2 containing the bounds, by default double-sided
#' @param print logical indicating whether the intervals are printed to the console
#' @param ... additional arguments passed to the function
#' @return a table with confidence intervals.
confint_int <- function (object, parm, level = 0.95, prob = c((1 - level)/2,  1 - (1 - level)/2), print = FALSE, ...)
{
  if (!isTRUE(all.equal(diff(prob), level, check.attributes = FALSE))) {
    warning("Incompatible arguments: `level` does not match `prob`.")
  }
  args <- list(...)
  if ("warn" %in% names(args) && is.logical(args$warn)) {
    warn <- args$warn
  }  else {
    warn <- TRUE
  }
  if (length(prob) != 2) {
    stop("`prob` must be a vector of size 2")
    prob <- sort(prob)
  }
  if (is.numeric(parm)) {
    ind <- parm
    parm <- c("profile", "tem", "modif.tem", "modif.empcov")[ind]
  } else {
    parm <- match.arg(arg = parm, choices = c("profile",
                                              "tem", "modif.tem", "modif.empcov", "r", "rstar"),
                      several.ok = TRUE)
    parm[parm %in% "r"] <- "profile"
    parm[parm %in% "rstar"] <- "tem"
    ind <- which(c("profile", "tem", "modif.tem", "modif.empcov") %in%
                   parm)
  }
  parm <- unique(parm)
  ind <- unique(ind[ind %in% 1:4])
  if (length(ind) == 0) {
    stop("Invalid `parm` argument.")
  }
  qulev <- qnorm(1 - prob)
  conf <- matrix(ncol = 4, nrow = 3)
  i = 1
  if (is.null(object$pll) && is.null(object$r)) {
    break
  }
  if (is.null(object$r)) {
    object$r <- sign(object$psi.max - object$psi) *
      sqrt(2 * (object$maxpll - object$pll))
  } else {
    object$r[is.infinite(object$r)] <- NA
  }
  if (is.null(object$normal)) {
    object$normal <- c(object$psi.max, object$std.error)
  }
  if (requireNamespace("cobs", quietly = TRUE)) {
    fit.r <- cobs::cobs(x = object$r, y = object$psi,
                        constraint = "decrease", lambda = 0, ic = "SIC",
                        pointwise = cbind(0, 0, object$normal[1]),
                        knots.add = TRUE, repeat.delete.add = TRUE,
                        print.mesg = FALSE, print.warn = FALSE)
    pr <- predict(fit.r, c(0, qulev))[, 2]
  } else {
    fit.r <- stats::smooth.spline(x = na.omit(cbind(object$r,
                                                    object$psi)), cv = FALSE)
    pr <- predict(fit.r, c(0, qulev))$y
    pr[1] <- object$normal[1]
  }
  conf[, i] <- pr
  if (warn) {
    if (!any(object$r > qnorm(prob[1]))) {
      warning("Extrapolating the lower confidence interval for the profile likelihood ratio test")
    }
    if (!any(object$r < qnorm(prob[2]))) {
      warning("Extrapolating the upper confidence interval for the profile likelihood ratio test")
    }
  }

  if (!is.null(conf)) {
    colnames(conf) <- c("Profile", "TEM", "Severini (TEM)",
                        "Severini (emp. cov.)")
    rownames(conf) <- c("Estimate", "Lower CI", "Upper CI")
    wrong_below <- which(conf[2, ] > conf[1, ])
    if (length(wrong_below) > 0) {
      conf[2, ][wrong_below] <- NA
    }
    wrong_above <- which(conf[3, ] < conf[1, ])
    if (length(wrong_above) > 0) {
      conf[3, ][wrong_above] <- NA
    }
    if (print) {
      cat("Point estimate for the parameter of interest psi:\n")
      cat("Maximum likelihood          :", round(object$psi.max,
                                                 3), "\n")
      cat("\n")
      cat("Confidence intervals, levels :", prob, "\n")
      cat("Wald intervals               :", round(object$psi.max +
                                                    sort(qulev) * object$std.error, 3), "\n")
      cat("Profile likelihood           :", round(conf[2:3,
                                                       1], 3), "\n")
    }
    return(invisible(conf[, ind]))
  }
}



#' Sample lifetime of semi-supercentenarians
#'
#' Given parameters of a generalized Pareto distribution, sampling window and
#' birth dates with excess lifetimes, sample new observations; excess lifetime
#' at \code{c1} are sampled from an exponential distribution, whereas
#' the birth dates are samples from a jittered histogram-based distribution
#' The new excess lifetime above the threshold are right-censored if they exceed
#' \code{c2}.
#'
#' @param n sample size
#' @param pars vector of length 2 containing the scale and shape parameters of the generalized Pareto distribution
#' @param xcal date at which individual reaches \code{u} years
#' @param c1 date, first day of the sampling frame
#' @param c2 date, last day of the sampling frame
#' @param slow excess lifetime at \code{c1} whenever \code{xcal} precedes the latter.
#' @return list with new birthdates (\code{xcal}), excess lifetime at \code{c1} (\code{slow}),
#' excess lifetime above \code{u} (\code{dat}) and right-censoring indicator (\code{rightcens}).
sample_lifetime <- function(n, pars, xcal, c1, c2, slow){

  sample_dates <- function(n, xcal, c1, c2, slow){
    sample_slow <- function(n, slow){
      sort(round(rexp(n, rate = 1/mean(365.25*slow[slow>0])))/365.25, decreasing = TRUE)
    }
    nday <- as.numeric(xcal[ind]-c1)
    nday <- nday[nday>0]
    nmax <- as.numeric(c2-c1)
    sslow <- sample_slow(round(sum(slow>0)*n/length(slow)), slow = slow)
    xhist <- hist(nday, plot = FALSE)
    bins <- with(xhist, sample(length(mids), n-length(sslow), p=density, replace=TRUE)) # choose a bin
    result <- round(runif(length(bins), xhist$breaks[bins], pmin(xhist$breaks[bins+1], nmax-1)))
    list(xcal = as.Date(round(c(-sslow*365.25, sort(result))), origin = "2009-01-01"),
         slow = c(sslow, rep(0, n-length(sslow))))
  }

  traject <- evd::rgpd(n = n, loc = 0, scale = pars[1], shape = pars[2])
  sdates <- sample_dates(n = n, xcal, c1, c2, slow)
  lifeah <- pmax(1,pmin(round(365.25*(traject + sdates$slow)), as.numeric(c2-sdates$xcal)))
  rcens_new <- lifeah == as.numeric(c2-sdates$xcal)
  list(xcal = sdates$xcal, slow = sdates$slow, dat = lifeah/365.25, rightcens = rcens_new)
}


#Likelihood root function, -r^2/2 is profile likelihood, with confidence intervals given by qchisq(p, 1)/2
rfun_endpt <- function(psi, thetahat, dat, rightcens, slow){
  llp <- prof_gpd_cens(psi = psi, dat = dat, rightcens = rightcens, slow = slow)
  sign(thetahat[1]-llp[2])*sqrt(-2*loglik(thetahat, dat = dat, rightcens = rightcens, slow = slow) - llp[1])
}
