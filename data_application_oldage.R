setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("script_oldage.R")
setwd("..")
library(mev)
library(evd)
library(tikzDevice)
options(tikzLatexPackages =
          c("\\usepackage{tikz}\n",
            "\\usepackage[active,tightpage,psfixbb]{preview}\n",
            "\\usepackage{amsmath}",
            "\\PreviewEnvironment{pgfpicture}\n",
            "\\setlength\\PreviewBorder{0pt}\n",
            "\\usepackage{fourier}\n"
          )
)
setTikzDefaults(overwrite = FALSE)

data(italcent, package = "mev") 
figures <- FALSE
dwidth <- 6
dheight <- 3.5
# Threshold for 105 years (criterion for inclusion in dataset
u <- 38351L
# Calendar date at which individual reaches 105 years
xcal <- italcent$birth + u
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-2009")
c2 <- lubridate::dmy("01-01-2016")

# Lower truncation level, zero if individual reached 105 between c1 and c2
slow <- as.numeric(pmax(0, c1 - xcal))
# Exceedances
dat <- italcent$numdays - u
rightcens <- italcent$rightcens
expo <- FALSE

#Oldest individual
max(italcent$numdays)/365.25

# Sequence of thresholds, from 105 years (in days), increments of half years
thresh <- seq(from = 0, by = 365.25/2, length = 14L)
# Fit a generalized Pareto model
param_gpd <- matrix(0, nrow = length(thresh), ncol = 9L)
colnames(param_gpd) <- c("loc","scale","shape","deviance", "scale.stderror", "shape.stderror","endpoint", "endpoint.stderror", "nu")
for(i in 1:length(thresh)){
  #For each threshold, compute new threshold exceedances
  datu <- dat - thresh[i]
  # Keep only exceedances
  ind <- which(datu > 0)
  vals <- optim(par = c(5, 0.01), fn = gpd_cens, method = "N", dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25, 
                expo = FALSE, hessian = TRUE)
  endpt_stderror <- sqrt(solve(obs.infomat(c(- vals$par[1]/vals$par[2], vals$par[2]), dat = (datu[ind])/365.25,rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25))[1,1])
  param_gpd[i,] <- c(loc = (u + thresh[i])/365.25, vals$par, -2*vals$value, sqrt(diag(solve(vals$hessian))),
                     ifelse(vals$par[2] >0, Inf, (u + thresh[i])/365.25 - vals$par[1]/vals$par[2]), endpt_stderror,  length(ind))
  #assign(paste0("opt",i), value = vals)
}


# Exponential model
param_exp <- matrix(0, nrow = length(thresh), ncol = 4L)
colnames(param_exp) <- c("loc","scale","deviance", "scale.stderror")
for(i in 1:length(thresh)){
  datu <- dat - thresh[i]
  ind <- which(datu > 0)
  vals <- optim(par = c(500), fn = gpd_cens, method = "Brent", lower = 0, upper = 800, dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = TRUE, hessian = TRUE)
  param_exp[i,] <- c(u + thresh[i], vals$par, -2*vals$value, sqrt(solve(vals$hessian)))
}

# difference in deviance between nested models
round(100*(1-pchisq(param_gpd[,"deviance"] - param_exp[,"deviance"], 1)),1)


# HOA - approximate ancillary stat
Vfun <- function(par, dat){
  endpt <- par[1]; xi <- par[2]; 
  cbind(dat/endpt,
        (dat - endpt)*log(-(dat - endpt)/endpt)/xi)
  }
# Weighted derivative wrt data
phi <- function(theta, dat, rightcens, V){
  endpt <- theta[1]; xi <- theta[2];
  stopifnot(endpt > max(dat))
  phi <- ifelse(rightcens, 
                -1/((dat - endpt)*xi),
                -(xi + 1)/((dat - endpt)*xi))
  colSums(phi*V)
  #phi
}
# Weighted mixed derivative wrt data + parameters
dphi <- function(theta, dat, rightcens, V){
  endpt <- theta[1]; xi <- theta[2];
  stopifnot(endpt > max(dat))
  dphi_endpt <- colSums(V*ifelse(rightcens, 
                                 -1/((dat^2 - 2*dat*endpt + endpt^2)*xi),
                                 -(xi + 1)/((dat^2 - 2*dat*endpt + endpt^2)*xi)))
  dphi_xi  <- colSums(V * 1/((dat - endpt)*xi^2))
  #same derivative wrt xi for censored/uncensored
  cbind(dphi_endpt, dphi_xi)
}
# Log-likelihood function in terms of parameters (endpt, xi)
loglik <- function(theta, dat, rightcens, slow){
  gpd_cens(par = c(-theta[2]*theta[1], theta[2]), dat = dat, 
           rightcens = rightcens, slow = slow)
}
#BN correction factor based on TEM approximation
qf <- function(theta, thetahat, dat, rightcens, loglik, slow){
  V <- Vfun(par = thetahat, dat = dat)
  num <-  det(cbind(
    phi(theta = thetahat, dat = dat, rightcens = rightcens, V = V) - 
      phi(theta = theta, dat = dat, rightcens = rightcens, V = V), 
    dphi(theta, dat = dat, rightcens = rightcens, V = V)[,2]))
  denom <- det(dphi(thetahat, dat = dat, rightcens = rightcens, V = V))
  J <- det(obs.infomat(theta = thetahat, dat = dat, rightcens = rightcens, slow = slow))
  #J <- det(numDeriv::hessian(loglik, x = thetahat, dat = dat, rightcens = rightcens, slow = slow))
  Jll <- #numDeriv::hessian(loglik, x = theta, dat = dat, rightcens = rightcens, slow = slow)[2,2]
    Jll <- obs.infomat(theta = theta, dat = dat, rightcens = rightcens, slow = slow)[2,2]
  return(num/denom*sqrt(J/Jll))
  #list(num, denom, J, Jll)
}
#rstar function; -rstar^2/2 is modified profile likelihood, with confidence intervals given by qchisq(p, 1)/2
rstarfun <- function(psi, thetahat, dat, rightcens, slow){
  llp <- prof_gpd_cens(psi = psi, dat = dat, rightcens = rightcens, slow = slow)
  theta <- llp[2:3]
  rr <- sign(thetahat[1]-llp[2])*sqrt(-2*loglik(thetahat, dat = dat, rightcens = rightcens, slow = slow) - llp[1])
  qp <- qf(theta = theta, thetahat = thetahat, dat = dat, rightcens = rightcens, loglik = loglik, slow = slow)
  (rr + (log(abs(qp))-log(abs(rr)))/rr)^2/2
}

#Likelihood root function, -r^2/2 is profile likelihood, with confidence intervals given by qchisq(p, 1)/2
rfun <- function(psi, thetahat, dat, rightcens, slow){
  llp <- prof_gpd_cens(psi = psi, dat = dat, rightcens = rightcens, slow = slow)
  sign(thetahat[1]-llp[2])*sqrt(-2*loglik(thetahat, dat = dat, rightcens = rightcens, slow = slow) - llp[1])
}

if(figures){
  setwd("figures")
  fig <- "italcent_profile.tex"
  tikz(fig, width = dwidth, height = dheight, standAlone = TRUE)
}
par(mar = c(5, 5, 1, 1), mfrow = c(1, 2))
js <- 1:4
# Consider two thresholds for which xi is negative
t1 <- matrix(0, nrow = length(js), ncol = 3)
for(j in js){
  # Data above threshold (in years) and left-truncation limits (in years)
  datuy <- (dat[dat>thresh[j]] - thresh[j])/365.25
  slowy <- pmax(0, slow - thresh[j])/365.25
  ind <- which(datuy > 0)
  thetahat <- c(-param_gpd[j,2]/param_gpd[j,3], param_gpd[j,3])
  # Find the root psi0 solving rstar(psi0)^2/2=0
  zero_rstar <- optim(fn = rstarfun, par = 25, method = "Brent", 
                      lower = max(datuy), upper = 150, thetahat = thetahat, 
                      dat = datuy[ind], 
                      rightcens = rightcens[ind], slow = slowy[ind])
  psi <- c(seq(max(datuy)+0.4, max(datuy) + 5, length = 30), seq(max(datuy) + 5.5, max(datuy)+500, length = 500))
  # Compute profile curve over range of endpt values
  if(is.finite(zero_rstar$value)){
    psi <- c(zero_rstar$par, psi)
  }
  prof <- t(sapply(psi, function(x){prof_gpd_cens(psi = x, dat = datuy[ind], 
                                             rightcens = rightcens[ind], 
                                             slow = slowy[ind])}))
  # plot(sort(prof[,2]), (prof[,1]-max(prof[,1]))[order(prof[,2])], 
  # ylim = c(-3,0), type = "l", panel.first = abline(h = -qchisq(0.95,1)/2))
  
  # Maximum likelihood estimate of (endpt, xi)
  devmle_ind <- which.max(c(-2*loglik(thetahat, dat = datuy[ind], rightcens = rightcens[ind], slow = slowy[ind]), prof[,1]))
  devmle <- max(c(-2*loglik(thetahat, dat = datuy[ind], rightcens = rightcens[ind], slow = slowy[ind]), prof[,1]))
  
  
  if(devmle_ind > 1){
    thetahat <- prof[devmle_ind-1, -1]
  } 
  # Correction q for TEM approximation
  qvec <- rep(0, nrow(prof))
  for(i in 1:nrow(prof)){
    qvec[i] <- qf(theta = prof[i,2:3], thetahat = thetahat, dat = datuy[ind], 
                  rightcens = rightcens[ind], loglik = loglik, slow = slowy[ind])
  }
  r <- sign(thetahat[1]-prof[,2])*sqrt(2*(devmle - prof[,1]))
  rstar <- r + (log(qvec/r))/r
  

  if(is.finite(zero_rstar$value)){
    rstar <- c(zero_rstar$value,rstar[-1])
  }
  corr <- list(r=r, 
               rstar = rstar, 
               psi = prof[,2], 
               q = qvec,
               psi.max = thetahat[1], 
               normal = as.vector(c(thetahat[1], param_gpd[j,"endpoint.stderror"])))
  corr <- spline.corr(corr, method = "smooth.spline")
  if(is.finite(zero_rstar$value)){
    corr$rstar[1] <- zero_rstar$value
    corr$tem.psimax <- zero_rstar$par
  }
  # Sort entries for plots
  od <- order(corr$psi)
  corr$psi <- corr$psi[od]
  corr$rstar <- corr$rstar[od]
  corr$r <- corr$r[od]
  corr$q <- corr$q[od]
  
  
  class(corr) <- "eprof"
  # plot(corr)
  # Superimpose maximum of rstar
  conftab <- mev:::confint.eprof(corr, method = "smooth.spline")
  if(j %in% 3:4){
    plot(corr$psi + (u+thresh[j])/365.25, -corr$r^2/2, type = 'l',
         xlab = "Limit to lifetime (in years)", 
         ylab = "profile log-likelihood", 
         bty = "l",ylim = c(-3,0), 
         xlim = c(120, conftab[3,2] + (u+thresh[j])/365.25) + 5,
         panel.first = abline(h = -qchisq(0.95,1)/2, col = "gray", lty = 3))
    lines(corr$psi + (u+thresh[j])/365.25, -corr$rstar^2/2, col = "gray60", lwd = 1.5, lty = 2)
    legend(x = "topright", legend = c("$\\ell_{\\mathsf{p}}$", "$\\ell_{\\mathsf{fr}}$"), col = c(1, "gray60"), lty = c(1,2), lwd = c(1,1.5), bty = "n")
  }
  printdate <- function(day){
    dday <- round((day - floor(day))*365.25,0)
    paste0("$", floor(day),"_y\\,", ifelse(dday<100, paste0("\\hphantom{0}",dday),dday), "_d$")}
  printdate <- function(day){ paste0("$", ifelse((day> 1e3)|| (is.na(day)), "\\infty", format(round(day,1), nsmall = 1)),"$")}
  conftab <- conftab + (u + thresh[j]) / 365.25
  r1z <- printdate(thetahat[1] + (u + thresh[j]) / 365.25)
  r1l <- printdate(conftab[2,1])
  r2l <- printdate(conftab[2,2])
  r2z <- printdate(conftab[1,2])
  r1k <- printdate(conftab[3,1])
  r2k <- printdate(conftab[3,2])
  t1[j,2:3] <- c(paste0(r1z, " (", r1l,", ", r1k ,")"),
                 paste0(r2z, " (", r2l,", ", r2k,")"))
  
}

if(figures){
  dev.off()
  system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
  setwd("..")
}
# Create table and caption
t1[,1] <- round(2*(u+thresh[1:length(js)])/ 365.25)/2
# t1 <- cbind(t1, sapply(thresh, function(th){sum(italcent$numdays>th+u)})[1:5])
colnames(t1) <- c("$u$", "$R$", "$R^*$")
library(xtable)
xtab <- xtable(t1, caption = c("Point estimates (95\\% confidence intervals) for the upper limit to lifetime $\\iota$ (in years) based on the profile likelihood ratio statistic $R(\\iota)$ (middle) and the modified likelihood ratio statistic $R^{*}(\\iota)$ for the tangent model approximation (right) using threshold exceedances of $u$ for the Italian semi-super centenarian data set. ", "Estimates and 95\\% confidence intervals for the upper limit to lifetime"),
               label = "tab_italcent_confint")
print(xtab, right = TRUE, booktabs = TRUE, include.rownames=FALSE, sanitize.text.function = identity, 
      sanitize.colnames.function = function(x){paste0("\\multicolumn{1}{c}{",x,"}")}, 
      sanitize.rownames.function = identity, table.placement = "t!",
      file = "tables/Table_italcent_confint.tex")

tab <- cbind(round(param_gpd[,'loc'],1), 
             param_gpd[,'nu'], 
             paste0(format(round(param_gpd[,'scale'],2), nsmall=2)," (", format(round(param_gpd[,'scale.stderror'],2),nsmall=2), ")"),
             paste0("$", format(round(param_gpd[,'shape'],2), nsmall=2),"\\; (", format(round(param_gpd[,'shape.stderror'],2),  nsmall = 2),")","$"),
             #round(param_gpd[,'endpoint'],0), 
             paste0("$",format(round(param_gpd[,'deviance']/2,1),nsmall=1), "$")
             #, paste0(format(round(param_exp[,'scale'],2), nsmall=2)," (", format(round(param_exp[,'scale.stderror'],2),nsmall=2), ")"),       
             # paste0("$",format(round(param_exp[,'deviance']/2,1),nsmall=1), "$")
)
#tab[is.infinite(param_gpd[,"endpoint"]),5] <- "$\\infty$"
colnames(tab) <- c("$u$", "$n_u$", "$\\widehat{\\sigma}$", "$\\widehat{\\xi}$", "$\\ell(\\widehat{\\boldsymbol{\\theta}})$")
library(xtable)
xtab <- xtable(tab[seq(1,by= 2, length.out =7),], caption = c("Maximum likelihood estimates of the generalized Pareto for the Italian super-centenarian data. From left to right, threshold $u$ (in years), number of threshold exceedances $n_u$, estimates (standard errors) of the scale $\\sigma$, shape $\\xi$ parameters, log-likelihood at \\textsc{mle} $\\ell(\\widehat{\\theta})$.", "Parameter estimates for the generalized Pareto distribution for the Italian super-centenarian data."),
               label = "tab_italcent_gpdparam", align = c("l","l",rep("r", ncol(tab)-1)))
print(xtab, right = TRUE, booktabs = TRUE, include.rownames=FALSE, sanitize.text.function = identity, 
      sanitize.colnames.function = function(x){paste0("\\multicolumn{1}{c}{",x,"}")}, 
      sanitize.rownames.function = identity, table.placement = "t!",
      file = "tables/Table_italcent_gpd.tex")





#######################################
##########   Bootstrap   ##############
#######################################
#DiCiccio and Young bootstrap scheme for threshold j=5
if(FALSE){
  set.seed(1234)
  tj  <- 5
  # Compute exceedances and left-truncation
  datuy <- (dat - thresh[tj])/365.25
  slowy <- pmax(0, slow - thresh[tj])/365.25
  ind <- which(datuy > 0)
  dev_mle <- param_gpd[tj, "deviance"]
  mle <- c(-param_gpd[tj, "scale"]/param_gpd[tj,"shape"], param_gpd[tj,"shape"])
  # Compute profile likelihood and graph of latter
  prof <- t(sapply(seq(max(datuy+1), max(datuy+0.1)+100, by = 0.5), 
                   function(x){prof_gpd_cens(psi = x, dat = datuy[ind], 
                                             rightcens = rightcens[ind], 
                                             slow = slowy[ind])}))
  wstat <- dev_mle - prof[,1]
  
  prof_gpd_cens <- function(psi, dat, rightcens, slow){
    opt <- optimize(f = function(lambda){gpd_cens(par = c(-lambda*psi, lambda), 
                                                  dat = dat,  rightcens = rightcens, slow = slow)}, 
                    interval = c(-1e-5, -0.4), tol = 1e-10)
    res <- -2*opt$objective
    attributes(res) <- list("param" = c(psi, opt$minimum))
    #res
    return(c(res, psi, opt$minimum))
  }
  ltrunc <- ifelse(xcal[ind] + thresh[tj] < c1, as.numeric(c1-xcal[ind]-thresh[tj]),0)/365.25
  nobs <- length(xcal[ind])
  N <- nrow(prof)
  wstatboot <- matrix(0, nrow = N, ncol = B + 1)
  mleboot <- array(0, dim = c(N, B + 1L, 2))
  constmleboot <- matrix(0, nrow = N, ncol = B + 1)
  for(i in 1:N){
    wstatboot[i,1] <- wstat[i]
    for(b in 1:B){
      set.seed(b+B*(i-1))
      Fa <- evd::pgpd(ltrunc, scale = -prof[i,2]*prof[i,3], shape = prof[i,3]) 
      simu_dat <-  evd::qgpd(Fa+runif(nobs)*(1-Fa), scale = -prof[i,2]*prof[i,3], shape = prof[i,3])
      simu_rightcens <-  xcal[ind] + thresh[tj] + simu_dat*365.25 > c2
      simu_dat <- ifelse(!simu_rightcens, simu_dat, as.numeric(c2-xcal[ind]-thresh[tj] )/365.25)
      #Compute MLE of the new data set - with no constraint on shape parameter
      opt <- try(optim(par = c(5, 0.01), fn = gpd_cens, method = "N", dat = simu_dat,
                       rightcens = simu_rightcens, slow = slowy, expo = FALSE))
      optc <- try(prof_gpd_cens(psi = prof[i,2], dat = simu_dat, rightcens = simu_rightcens, slow = slowy[ind]))
      if(!is.character(opt)){
        mleboot[i,b+1,] <- opt$par
      } else{
        mleboot[i,b+1,] <- rep(NA, 2)
      }
      # Compute likelihood root
      
      if(!is.character(optc)){
        constmleboot[i, b+1] <- optc[3]
        wb <- try(-2*opt$value - optc[1])
        if(!is.character(wb)){
          wstatboot[i, b+1] <- wb
        } else{
          wstatboot[i,b+1] <- NA
        }
      } else { 
        constmleboot[i, b+1]  <- NA
      }
    }
  }
  bootpval <- (B + 1 - apply(wstatboot, 1, function(x){rank(x, na.last = NA)[1]}))/(B+1)
  save(bootpval, constmleboot, wstatboot, mleboot, file = "Bootstrap_Young.RData") 
  
  if(figures){
    setwd("figures")
    fig <- "italcent_pvalfun.tex"
    tikz(fig, width = dwidth, height = 0.7*dheight, standAlone = TRUE)
    par(mar = c(4,4,1,1))
    plot(prof[,2]+(u+thresh[tj])/365.25, 1-pchisq(wstat, 1),type = "l", xlab = "$\\iota$ (in years)", 
         ylab = "$P$-value function of $R^2(\\iota)$", bty = "l",ylim = c(0,1), yaxs = "i")
    lines(prof[,2]+(u+thresh[tj])/365.25, 1-bootpval/5000, col = "gray60", lwd = 1.5)
    legend(x = "bottomright", legend = c("asymptotic", "bootstrap"), col = c(1, "gray60"), lwd = c(1,1.5), bty = "n")
    
    dev.off()
    system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
    setwd("..")
  }
}
