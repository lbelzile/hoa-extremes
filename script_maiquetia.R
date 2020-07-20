# Weighted HOA application

gpdNw.ll <- function(par, dat, w, N) {
  xi = par[2]
  z = par[1]
  euler_gamma = 0.57721566490153231044
  sigma = ifelse(abs(xi) > 1e-8,
                 z * xi/(exp(lgamma(N + 1) + lgamma(-xi + 1) - lgamma(N - xi + 1)) - 1),
                 z / (euler_gamma + psigamma(N + 1)))
  gpdw.ll(par = c(sigma, xi), w = w, dat = dat)
}

gpdw.ll <- function(par, dat, w, tol = 1e-05) {
  sigma = par[1]
  if(sigma < 0){
    return(-Inf)
  }
  xi = par[2]
  n <- length(dat)
  if (abs(xi) > tol) {
    if(xi > -1){
      #Most scenarios
      sum(w * as.vector(-log(sigma) - (1 + 1/xi) * log(pmax(1 + xi/sigma * dat, 0))))
    } else if(isTRUE(all.equal(xi, -1, check.attributes = FALSE))){
      - sum(log(max(dat)) * w)
    } else{
      -Inf
    }
  } else {
    sum(w * (-n * log(sigma) - dat/sigma))
  }
}

gpdNw.score <- function(par, dat, w, N) {
  z = par[1]
  xi = par[2]
  cst <- exp(lgamma(N + 1) + lgamma(1 - xi) - lgamma(N + 1 - xi))
  xizero <- abs(xi) < 1e-6
  if(!xizero){
  as.vector(c(sum(w*(dat * (cst - 1) * (1/xi + 1)/(z^2 * (dat * (cst - 1)/z + 1)) - 1/z)), sum(w*(-(psigamma(N - xi +
  1) * cst - psigamma(-xi + 1) * cst) * dat * (1/xi + 1)/(z * (dat * (cst - 1)/z + 1)) + ((psigamma(N -
  xi + 1) * cst - psigamma(-xi + 1) * cst) * xi * z/(cst - 1)^2 - z/(cst - 1)) * (cst - 1)/(xi *z) +
  log(dat * (cst - 1)/z + 1)/xi^2))))
  } else{
  euler_gamma <- 0.57721566490153231044
  psi1pN <- psigamma(1 + N)
  psip1pN <- psigamma(1 + N, deriv = 1)
  as.vector(c( sum(w*(dat*euler_gamma - z + dat*psi1pN))/z^2,
  sum(w*((6*dat^2*euler_gamma^3 - 12*dat*euler_gamma^2*z - 6*dat*euler_gamma^3*z -
  dat*euler_gamma*pi^2*z + 6*euler_gamma^2*z^2 + pi^2*z^2 +
  6*(3*dat^2*euler_gamma - dat*(2 + 3*euler_gamma)*z + z^2)*
  psigamma(1 + N)^2 +
  6*dat*(dat - z)*psigamma(1 + N)^3 + 6*(dat*euler_gamma - z)*z*
  psigamma(deriv = 1, 1+N) + psigamma(1 + N)*(18*dat^2*euler_gamma^2 -
  dat*(24*euler_gamma + 18*euler_gamma^2 + pi^2)*z + 12*euler_gamma*z^2 +
  6*dat*z*psigamma(deriv = 1, 1+N)))/(12*z^2*(euler_gamma + psigamma(1 + N)))))))
  }
}

gpdw.score <- function(par, w, dat) {
  sigma = par[1]
  xi = as.vector(par[2])
  if (!isTRUE(all.equal(0, xi, check.attributes = FALSE))) {
    c(sum(w*(dat * (xi + 1)/(sigma^2 * (dat * xi/sigma + 1)) - 1/sigma)), sum(w*(-dat * (1/xi + 1)/(sigma *
                                                                                               (dat * xi/sigma + 1)) + log(pmax(dat * xi/sigma + 1, 0))/xi^2)))
  } else {
    c(sum(w*((dat - sigma)/sigma^2)), sum(w*(1/2 * (dat - 2 * sigma) * dat/sigma^2)))
  }
}

gpdNw.infomat <- function(par, dat, N, w, method = c("obs", "exp"), nobs = length(dat)) {
  z = as.vector(par[1])
  xi = as.vector(par[2])
  if (xi < -0.5) {
    return(matrix(NA, 2, 2))
  }
  method <- method[1]  #default corresponds to observed information rather than Fisher information
  xizero <- abs(xi) < 1e-4
  # Fisher information matrix
  cst <- exp(lgamma(N + 1) + lgamma(1 - xi) - lgamma(N + 1 - xi))
  if (method == "exp") {
    if(!xizero){
      sigmaf <- z * xi/(cst - 1)
      Jac <- rbind(c(xi/(cst - 1),
                     -(digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst) * xi * z/(cst - 1)^2 + z/(cst - 1)), c(0, 1))
    } else{
      euler_gamma <- 0.57721566490153231044
      psi1pN <- psigamma(1 + N)
      sigmaf <- z/(euler_gamma + psi1pN)
      Jac <- rbind(c(1/(euler_gamma + psi1pN),
                     -(z*(6*euler_gamma^2 + pi^2 + 12*euler_gamma*psigamma(1 + N) +
                            6*psigamma(1 + N)^2 - 6*psigamma(1 + N, deriv = 1)))/
                       (12*(euler_gamma + psigamma(1 + N))^2)), c(0,1))
    }
    infomat <- t(Jac) %*% gpd.infomat(par = c(sigmaf, xi), dat = dat, method = "exp") %*% Jac
  } else if (method == "obs") {
    # Observed information
    if(!xizero){
k11 <- sum(w*(2 * dat * (cst - 1) * (1/xi + 1)/(z^3 * (dat * (cst - 1)/z + 1)) - dat^2 * (cst -
  1)^2 * (1/xi + 1)/(z^4 * (dat * (cst - 1)/z + 1)^2) - 1/z^2))
k12 <- sum(w*(-(digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst) * dat * (1/xi + 1)/(z^2 * (dat *
  (cst - 1)/z + 1)) + (digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst) * dat^2 * (cst -
  1) * (1/xi + 1)/(z^3 * (dat * (cst - 1)/z + 1)^2) + dat * (cst - 1)/(xi^2 * z^2 * (dat *
  (cst - 1)/z + 1))))
k22 <- sum(w*(-(digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst)^2 * dat^2 * (1/xi + 1)/(z^2 *
  (dat * (cst - 1)/z + 1)^2) + (digamma(N - xi + 1)^2 * cst - 2 * digamma(N - xi + 1) * digamma(-xi +
  1) * cst + digamma(-xi + 1)^2 * cst - psigamma(deriv = 1, N - xi + 1) * cst + psigamma(deriv = 1,
  -xi + 1) * cst) * dat * (1/xi + 1)/(z * (dat * (cst - 1)/z + 1)) - (digamma(N - xi + 1) *
  cst - digamma(-xi + 1) * cst) * ((digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst) * xi *
  z/(cst - 1)^2 - z/(cst - 1))/(xi * z) + (2 * (digamma(N - xi + 1) * cst - digamma(-xi + 1) *
  cst)^2 * xi * z/(cst - 1)^3 - (digamma(N - xi + 1)^2 * cst - 2 * digamma(N - xi + 1) * digamma(-xi +
  1) * cst + digamma(-xi + 1)^2 * cst - psigamma(deriv = 1, N - xi + 1) * cst + psigamma(deriv = 1,
  -xi + 1) * cst) * xi * z/(cst - 1)^2 - 2 * (digamma(N - xi + 1) * cst - digamma(-xi + 1) *
  cst) * z/(cst - 1)^2) * (cst - 1)/(xi * z) + ((digamma(N - xi + 1) * cst - digamma(-xi +
  1) * cst) * xi * z/(cst - 1)^2 - z/(cst - 1)) * (cst - 1)/(xi^2 * z) - 2 * (digamma(N - xi +
  1) * cst - digamma(-xi + 1) * cst) * dat/(xi^2 * z * (dat * (cst - 1)/z + 1)) + 2 * log(dat *(cst - 1)/z + 1)/xi^3))
    } else{
      euler_gamma <- 0.57721566490153231044
      psigamma1pN <- psigamma(1 + N)
      zeta3 <- 1.2020569031595944587
      psigammap1pN <- psigamma(1 + N, deriv = 1)
      psigammapp1pN <-  psigamma(1 + N, deriv = 2)
      k11 <- -sum(w*(-2*dat*euler_gamma + z - 2*dat*psigamma1pN))/z^3
      k12 <- -sum(w*((1/(12*z^3))*(dat*(-12*dat*euler_gamma^2 + (6*euler_gamma*(2 + euler_gamma) + pi^2)*z +
                                       12*(-2*dat*euler_gamma + z + euler_gamma*z)*psigamma1pN +
                                       6*(-2*dat + z)*psigamma1pN^2 -  6*z*psigammap1pN))))
      k22 <- -sum(w*((-96*dat^3*euler_gamma^5 + 24*dat^2*euler_gamma^3*(6*euler_gamma*(1 + euler_gamma) + pi^2)*z -
                     24*dat*euler_gamma^2*(2*euler_gamma^2*(3 + euler_gamma) + (1 + euler_gamma)*pi^2)*
                     z^2 + (12*euler_gamma^4 + 12*euler_gamma^2*pi^2 - pi^4)*z^3 -
                     12*((40*dat^3*euler_gamma + 4*dat*(3 + 5*euler_gamma)*z^2 - z^3 -
                            12*dat^2*(z + 5*euler_gamma*z))*psigamma1pN^4 +
          4*dat*(2*dat^2 - 3*dat*z + z^2)*psigamma1pN^5 +
          2*psigamma1pN^3*(40*dat^3*euler_gamma^2 -  dat^2*(12*euler_gamma*(2 + 5*euler_gamma) + pi^2)*z +
          dat*(4*euler_gamma*(6 + 5*euler_gamma) + pi^2)*z^2 -2*euler_gamma*z^3 + 6*dat*(dat - z)*z*psigammap1pN) +
          z*((12*dat^2*euler_gamma^3 - 12*dat*euler_gamma^2*(1 + euler_gamma)*z + (6*euler_gamma^2 - pi^2)*z^2)*
          psigammap1pN + 3*z^2*psigammap1pN^2 + 4*euler_gamma*(dat*euler_gamma - z)*z*
          (psigammapp1pN + 2*zeta3)) + psigamma1pN^2*(80*dat^3*euler_gamma^3 -
          6*dat^2*euler_gamma*(4*euler_gamma*(3 + 5*euler_gamma) + pi^2)*z + 2*dat*(4*euler_gamma^2*(9 + 5*euler_gamma) +
          (1 + 3*euler_gamma)*pi^2)*z^2 - (6*euler_gamma^2 + pi^2)*z^3 +6*z*(6*dat^2*euler_gamma - 2*dat*(1 + 3*euler_gamma)*z + z^2)*
          psigammap1pN + 4*dat*z^2*(psigammapp1pN + 2*zeta3)) +2*psigamma1pN*(euler_gamma*(20*dat^3*euler_gamma^3 -
          3*dat^2*euler_gamma*(2*euler_gamma*(4 + 5*euler_gamma) + pi^2)*z + dat*(2*euler_gamma^2*(12 + 5*euler_gamma) +
          (2 + 3*euler_gamma)*pi^2)*z^2 - (2*euler_gamma^2 + pi^2)*z^3) +6*euler_gamma*z*(3*dat^2*euler_gamma - dat*(2 + 3*euler_gamma)*z + z^2)*
          psigammap1pN - 2*z^2*(-2*dat*euler_gamma + z)*(psigammapp1pN + 2*zeta3))))/(144*z^3*(psigamma1pN + euler_gamma)^2)))
    }
    infomat <- cbind(c(k11, k12), c(k12, k22))
  }
  colnames(infomat) <- rownames(infomat) <- c("Nmean","shape")
  return(infomat)
}


#' Canonical parameter in the local exponential family approximation
#'
#' @inheritParams gpdN
#' @rdname gpdN.temstat
#' @keywords internal
#' @export
gpdNw.phi <- function(par, dat, w, N, V) {
  xi = par[2]
  z = par[1]
  cst <- exp(lgamma(N + 1) + lgamma(-xi + 1) - lgamma(N - xi + 1))
  matrix(-(cst - 1) * (1/xi + 1)/(z * (dat * (cst - 1)/z + 1)), nrow = 1) %*% (w*V)
}

## Derivative matrix of the canonical parameter in the local exponential family approximation
#' @inheritParams gpdN
#' @rdname gpdN.temstat
#' @keywords internal
#' @export
gpdNw.dphi <- function(par, dat, N, V, w) {
  xi = par[2]
  z = par[1]
  cst <- exp(lgamma(N + 1) + lgamma(-xi + 1) - lgamma(N - xi + 1))
rbind((cst - 1) * (1/xi + 1)/(z^2 * (dat * (cst - 1)/z + 1)) - dat * (cst - 1)^2 * (1/xi + 1)/(z^3 *
(dat * (cst - 1)/z + 1)^2), -(psigamma(N - xi + 1) * cst - psigamma(-xi + 1) * cst) * (1/xi +
1)/(z * (dat * (cst - 1)/z + 1)) + (psigamma(N - xi + 1) * cst - psigamma(-xi + 1) * cst) * dat *
(cst - 1) * (1/xi + 1)/(z^2 * (dat * (cst - 1)/z + 1)^2) + (cst - 1)/(xi^2 * z * (dat * (cst -
1)/z + 1))) %*% (w*V)
}


gpdNw.pll <- function(psi, mod = c("profile","tem","modif"), w,
                      dat, N, correction = TRUE, threshold = NULL, plot = TRUE,
                      ...){
    mle <- optim(par = fit.gpd(xdat = dat, threshold = ifelse(is.null(threshold), 0, threshold))$estimate,
                 fn = function(x, dat, w){-gpdw.ll(dat = dat, par = x, w = w)},
                 gr =function(x, dat, w){-gpdw.score(dat = dat, par = x, w = w)},
                 dat = dat, w = w)$par
    mle <- c((exp(lgamma(N + 1) + lgamma(1 - mle[2]) - lgamma(N + 1 - mle[2])) - 1) * mle[1]/mle[2], mle[2])
    names(mle) <- c("Nmean","shape")
  param <- "Nmean"
  mod <- match.arg(mod, c("profile", "tem", "modif"), several.ok = TRUE)
  #If there is a threshold
  if(!is.null(threshold)){
    stopifnot(is.numeric(threshold), length(threshold) == 1)
    if(min(dat) < threshold){
      dat <- dat[dat>threshold] - threshold
    } else {
      dat <- dat - threshold
    }
  } else {
    threshold <- 0
  }
  # Arguments for parametrization of the log likelihood
  args <- c(param, "shape")
   # Sanity checks to ensure all arguments are provided
  if (missing(N)) {
    stop("Argument `N` missing. Procedure aborted")
  }
  xmin <- min(dat)
  xmax <- max(dat)
  shiftres <- TRUE
  # If maximum likelihood estimates are not provided, find them

  # Extract the components, notably V for model `tem`. Keep other components for optimization
  V <- gpdN.Vfun(par = mle, dat = dat, N = N)

  maxll <- gpdNw.ll(mle, dat = dat, N = N, w=w)
  std.error <- sqrt(solve(gpdNw.infomat(par = mle, dat = dat, method = "exp", N = N))[1])
  constr.mle.Nmean <- function(Nmeant) {
    suppressWarnings(alabama::auglag(par = 0.01, function(lambda, psi, N) {
      -gpdNw.ll(par = c(psi, lambda), dat = dat, N = N, w = w)
    }, gr = function(lambda, psi, N) {
      -gpdNw.score(par = c(psi, lambda), dat = dat, N = N, w = w)[2]
    },
    hin = function(lambda, psi, N){
      sigma = ifelse(abs(lambda > 1e-8),
                     psi * lambda/(exp(lgamma(N + 1) + lgamma(-lambda + 1) - lgamma(N - lambda + 1)) - 1),
                     psi / (0.57721566490153231044 + psigamma(N + 1)))
      if(lambda >= 0){c(1e-8, sigma, 1-lambda, lambda + 1)} else{ c(- sigma/lambda - xmax, sigma, 1 - lambda, lambda + 1)}},
    psi = Nmeant, N = N, control.outer = list(trace = FALSE))$par)
  }
  # Missing psi vector
  if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
    #compute profile log-lik on a grid left and right of the MLE
    psirangelow <- unique(pmax(mean(dat), seq(-3, -0.25, length = 20) * std.error +
                                 mle[1]))
    lowvals <- sapply(psirangelow, function(par) {
      gpdNw.ll(c(par, constr.mle.Nmean(par)), dat = dat, N = N, w = w)
    }) - maxll
    psirangehigh <- seq(2.5, 50, length = 20) * std.error + mle[1]
    highvals <- sapply(psirangehigh, function(par) {
      gpdNw.ll(c(par, constr.mle.Nmean(par)), dat = dat, N = N, w = w)
    }) - maxll
    #Try to do linear interpolation - only works if value inside the interval lowvals or highvals
    lo <- approx(x = lowvals, y = psirangelow, xout = -4)$y
    #Else linear interpolation with linear model fitted to lower values
    lo <- ifelse(is.na(lo), spline(x = lowvals, y = psirangelow, xout = -4)$y, lo)
    psirangelow <- seq(lo, mle[1], length = 20)
    lowvals <- sapply(psirangelow, function(par) {
      gpdNw.ll(c(par, constr.mle.Nmean(par)), dat = dat, N = N, w = w)
    }) - maxll
    #hi <- approx(x = highvals, y = psirangehigh, xout = -4)$y
    #For upper, use spline approx
    hi <- spline(x = highvals, y = psirangehigh, xout = -4)$y
    #Recompute the range
    psirangehigh <- seq(psirangehigh[1], hi, length = 30)
    highvals <- sapply(psirangehigh, function(par) {
      gpdNw.ll(c(par, constr.mle.Nmean(par)), dat = dat, N = N, w = w)
    }) - maxll
    #Remove NAs, inf, etc.
    highvals <- highvals[is.finite(highvals)]
    psirangehigh <- psirangehigh[1:length(highvals)]
    #If could not interpolate, use a simple linear model to predict lower value
    #hi <- ifelse(is.na(hi), spline(x = highvals, y = psirangehigh, xout = -4)$y, hi)
    psi <- as.vector(c(spline(x=c(lowvals, 0), y = c(psirangelow, mle[1]), xout = seq(-4, -0.1, length = 15))$y, mle[1],
                       spline(x=c(0, highvals), y = c(mle[1], psirangehigh), xout = seq(-0.1, highvals[length(highvals)], length = 20))$y))
  } else{
    psi <- psi - threshold
  }

  if (any(as.vector(psi) < 0)) {
    warning("Negative Nmean values provided.")
    psi <- psi[psi > 0]
    if (length(psi) == 0) {
      psi <- mle[1]
    }
  }

  pars <- cbind(psi, sapply(psi, constr.mle.Nmean))
  # Profile log likelihood values for psi
  profll <- apply(pars, 1, function(par) {
    gpdNw.ll(par = par, dat = dat, N = N, w = w)
  })
  r <- sign(mle[param] - psi) * sqrt(2 * (maxll - profll))
  if ("tem" %in% mod) {
    phi.mle <- gpdNw.phi(par = mle, dat = dat, N = N, V = V, w = w)
    q2num <- apply(pars, 1, function(par) {
      det(rbind(c(c(phi.mle) - gpdNw.phi(par = par, dat = dat, V = V, N = N, w = w)),
                gpdNw.dphi(par = par, dat = dat, V = V, N = N, w = w)[-1, ]))
    })
    if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
      warning("Correction factor and likelihood root are of opposite sign - check output")
    }

    logq <- apply(pars, 1, function(par) {
      -0.5 * log(gpdNw.infomat(par = par, dat = dat, method = "obs", N = N, w = w)[-1, -1])
    }) + log(abs(q2num))
    qmlecontrib <- -log(det(gpdNw.dphi(par = mle, dat = dat, V = V, N = N, w = w))) + 0.5 *
      log(det(gpdNw.infomat(par = mle, dat = dat, method = "obs", N = N, w = w)))
    logq <- logq + qmlecontrib
    qcor <- sign(q2num) * exp(logq)
    rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r)))/r)

    tem.max.opt <- function(psi, dat = dat) {
      para <- c(psi, constr.mle.Nmean(psi))
      pll <- gpdNw.ll(par = para, dat = dat, N = N, w = w)
      rs <- 2 * (maxll - pll)
      logq <- -0.5 * log(gpdNw.infomat(par = para, dat = dat, method = "obs", N = N, w = w)[-1,-1]) + qmlecontrib + log(abs(det(rbind(c(c(phi.mle) -
   gpdNw.phi(par = para,dat = dat, V = V, N = N, w=w)), gpdNw.dphi(par = para, dat = dat, V = V, N = N, w=w)[-1, ]))))
      rs + logq - log(sqrt(abs(rs)))
    }
    tem.max <- optim(par = mle[1], fn = tem.max.opt, method = "Brent", dat = dat, lower = max(1e-05, mle[1] - std.error), upper = mle[1] + std.error, control = list(abstol = 1e-10))$par

  }
  if ("modif" %in% mod) {
    # Tangent exponential model approximation of Fraser and Reid to the profile likelihood
    tem.objfunc.Nmean <- function(par) {
      0.5 * log(gpdNw.infomat(par = par, w = w, dat = dat, method = "obs", N = N)[2, 2]) -
        log(abs(gpdNw.dphi(par = par, dat = dat, N = N, w =w, V = V[, 2, drop = FALSE])[2, 1]))
    }
    optim.tem.fn.Nmean <- function(psi) {
      theta.psi.opt <- constr.mle.Nmean(psi)
      param <- c(psi, theta.psi.opt)
      ll <- gpdNw.ll(param, dat = dat, N = N, w = w)
      ll + tem.objfunc.Nmean(param)
    }
    # TEM profile log likelihood values for psi
    proflltem <- profll + suppressWarnings(apply(pars, 1, tem.objfunc.Nmean))
    # Maximum objective function for TEM
    tem.mle.opt <- optim(par = mle[1], fn = optim.tem.fn.Nmean, method = "Brent", lower = max(1e-05, mle[1] - std.error), upper = mle[1] + std.error, control = list(fnscale = -1))
    tem.mle <- c(tem.mle.opt$par, constr.mle.Nmean(tem.mle.opt$par))
    # Severini empirical covariance function adjustment to profile likelihood
    gpdN.score.f <- function(par, dat, N) {
      z = par[1]
      xi = par[2]
      cst <- exp(lgamma(N + 1) + lgamma(1 - xi) - lgamma(N + 1 - xi))
      -(psigamma(N - xi + 1) * cst - psigamma(-xi + 1) * cst) * dat * (1/xi + 1)/(z *(dat * (cst - 1)/z + 1)) +
        ((psigamma(N - xi + 1) * cst - psigamma(-xi +1) * cst) * xi * z/(cst - 1)^2 - z/(cst - 1)) * (cst - 1)/(xi * z) + log(dat *(cst - 1)/z + 1)/xi^2
    }
    # Score at MLE (sums to zero)
    score.Nmean.mle <- w*gpdN.score.f(mle, dat, N)
    empcov.objfunc.Nmean <- function(par) {
      0.5 * log(gpdNw.infomat(par = par, w = w, dat = dat, method = "obs", N = N)[2, 2]) -
        log(abs(sum(score.Nmean.mle * w*gpdN.score.f(par, dat, N = N))))
    }
    profllempcov <- profll + suppressWarnings(apply(pars, 1, empcov.objfunc.Nmean))
    optim.empcov.fn.Nmean <- function(psi) {
      theta.psi.opt <- constr.mle.Nmean(psi)
      param <- c(psi, theta.psi.opt)
      ll <- gpdNw.ll(param, dat = dat, N = N, w = w)
      ll + empcov.objfunc.Nmean(param)
    }
    empcov.mle.opt <- optim(par = mle[1], fn = optim.empcov.fn.Nmean, method = "Brent", 
                            lower = max(1e-05, mle[1] - std.error), upper = mle[1] + std.error, control = list(fnscale = -1))
    empcov.mle <- c(empcov.mle.opt$par, constr.mle.Nmean(empcov.mle.opt$par))
  }

# Return profile likelihood and quantities of interest (modified likelihoods)
colnames(pars) <- names(mle)
ans <- list(mle = mle, pars = pars, psi.max = as.vector(mle[param]), 
            param = param, std.error = std.error,
            psi = psi, pll = profll, maxpll = maxll, r = r)
# Shift by threshold if non-null
if(shiftres){
  ans$psi <- ans$psi + threshold
  ans$mle[1] <- ans$psi.max <- ans$mle[1] + threshold
  ans$pars[,1] <- ans$pars[,1] + threshold
}

if ("tem" %in% mod) {
  ans$q <- qcor
  ans$rstar <- rstar
  ans$tem.psimax <- as.vector(tem.max) + ifelse(shiftres, threshold, 0)
  ans$normal <- c(ans$psi.max, ans$std.error)
   if (correction && length(psi) > 10) {
     ans <- spline.corr(ans)
   }
}
if ("modif" %in% mod) {
  ans$tem.mle <- ifelse(param == "shape", tem.mle[2], tem.mle[1])
  if(shiftres){ ans$tem.mle[1] <- ans$tem.mle[1] + threshold }
  ans$tem.pll <- proflltem
  ans$tem.maxpll <- as.vector(tem.mle.opt$value)
  ans$empcov.mle <- ifelse(param == "shape", empcov.mle[2], empcov.mle[1])
  if(shiftres){ ans$empcov.mle[1] <- ans$empcov.mle[1] + threshold }
  ans$empcov.pll <- as.vector(profllempcov)
  ans$empcov.maxpll <- as.vector(empcov.mle.opt$value)
}
if ("tem" %in% mod) {
  class(ans) <- c("eprof", "fr")
} else {
  class(ans) <- "eprof"
}
ans$family <- "gpd"
ans$threshold <- threshold
if(plot){
  plot(ans)
}
return(invisible(ans))

}

