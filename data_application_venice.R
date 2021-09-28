setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Install packages mev, devtools, ggplot2, eva, qgam
library(mev)
library(ggplot2)
if(!require(eva)){
  install.packages("devtools")
  devtools::install_github("lbelzile/eva")
}
library(tikzDevice)
options(tikzLatexPackages =
          c("\\usepackage{tikz}\n",
            "\\usepackage[active,tightpage,psfixbb]{preview}\n",
            "\\usepackage{amsmath}",
            "\\PreviewEnvironment{pgfpicture}\n",
            "\\setlength\\PreviewBorder{0pt}\n",
            "\\usepackage{fourier}\n",
            "\\DeclareMathAlphabet{\\mathdis}{OT1}{pag}{m}{n}\n"
          )
)
# Switch this boolean to TRUE to generate and save figures
figures <- FALSE
fig_dir <- "figures"
if(!dir.exists("figures")){
  dir.create(fig_dir)
}
#data(venice, package = "mev") # mev version > 1.13.300
load("data/venice.rda")
ny <- nrow(venice)
if(figures){
  setwd("figures")
  fig <- "venice_series.tex"
  tikz(fig, height = 3, width = 6, standAlone = TRUE)
}
  par(mar = c(4,4,1,1), bty = "l")
  matplot(venice$year, cbind(venice$r1, venice$r2),
          pch = c(20,18), col = c("black","grey"),
          ylab = "sea level height (in cm)", xlab = "year")
  venicep <- with(venice, data.frame(sealev = c(r1,r2), year = c(year, year), os = factor(c(rep(1, ny), rep(2, ny)))))
  gamf <- qgam::qgam(sealev ~ s(year, k = 50) + os, qu = 0.5, data = venicep, )
  lines(venice$year, predict(gamf, venicep[1:ny,]), lty = 2)
  # lines(venice$year, predict(lm(r1 ~ year + I((year>1982)*(year-1982)), data = venice)), lty = 1, lwd = 2)
  # lines(venice$year, predict(lm(r1 ~ year + I(year>1982), data = venice)), lty = 1, lwd = 2)
  lines(venice$year, predict(lm(sealev ~ year + os, data = venicep), venicep[1:ny,]), lty = 1, lwd = 2)
  legend(x="topleft", lty = 2:1, lwd = 1:2, legend = c("quantile reg.","linear"), bty = "n")
  if(figures){
  dev.off()
  system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
  setwd("..")
}
# Residual plot
# car::residualPlot(lm(r1 ~ year, data = venice))

source("script_venice.R")

z <- 194 # highest level in 1966
yr <- scale(venice$year)
dat <- as.matrix(venice[-36,2:3])
yrs <- as.numeric(yr)[-36]
#Remove 1926 because there is only one observation and the code doesn't (yet) handle NAs
t0s <- (seq(2020, 2050, by = 5L) - attr(yr,"scaled:center"))/attr(yr,"scaled:scale")
confint_array <- array(NA, dim = c(length(t0s), 3, 2, 2))
for(t0i in 1:length(t0s)){
  t0 <- t0s[t0i]
  # Initial MLE estimates
  mlev <- ismev::rlarg.fit(show = FALSE, r=2, xdat = dat, ydat = data.frame(yrs), mul = 1)$mle
  pexc <- pex(sigma = mlev[3], mu = mlev[1], mut = mlev[2], xi = mlev[4], z = z, t0 = t0)
  # stopifnot(pexc == pgev(q = z, loc = mlev[1]+mlev[2]*t0, scale = mlev[3], shape = mlev[4], lower.tail = FALSE))
  pars <- c(mlev[2], pexc, mlev[3:4])
  mle <- optim(fn = function(x){-rlargp.ll(par = x, dat = dat, z = z, t = yrs, t0 = t0)},
              hessian = TRUE, par = pars, control = list(parscale = pars, reltol = 1e-10), method = "BFGS")
  # refine the maximum likelihood estimates so they solve the score equation
  numDeriv::grad(function(x){rlargp.ll(par = x, dat = dat, z = z, t = yrs, t0 = t0)}, x = mle$par)
  mlev <- mle$par

  # compute the profile log-likelihood
  out <- list()
  psi <- seq(0.0001, 0.035, length.out = 1001)
  npsi <- length(psi)
  cmle <- matrix(NA, nrow = npsi, ncol = 4)
  cmle[,2] <- psi
  clogLik <- q <- Jrest <- numeric(length = npsi)
  V <- rlargp.V(dat = dat, par = mle$par, t = yrs, t0 = t0, z=z)
  detDbot <- det(rlargp.dphi(dat = dat, par = mle$par, z = z, V = V, t = yrs, t0 = t0))
  phiMle <- rlargp.phi(dat = dat, par = mle$par, z = z, V = V, t = yrs, t0 = t0)
  detInfoMle <- det(mle$hessian)
  for(i in 1:npsi){
    opt <- try(optim(fn = function(x){-rlargp.ll(par = c(x[1], psi[i], x[-1]), dat = dat, z = z, t = yrs, t0 = t0)},
                  par = mle$par[-2],
                  gr = function(x){-rlargp.score(par = c(x[1], psi[i], x[-1]), dat = dat, z = z, t = yrs, t0 = t0)[-2]},
                    control = list(parscale = mle$par[-2], reltol = 1e-10), method = "BFGS"))
    if(!is.character(opt) && opt$convergence == 0){
      clogLik[i] <- -opt$value
      cmle[i,-2] <- opt$par
      Jrest[i] <- det(rlargp.infomat(par = c(opt$par[1], psi[i], opt$par[-1]), dat = dat, z = z, t = yrs, t0 = t0)[-2,-2])
      topMat <- rlargp.dphi(dat = dat, par = cmle[i,], z = z, V = V, t = yrs, t0 = t0)
      topMat[2,] <-  phiMle - rlargp.phi(dat = dat, par = cmle[i,], z = z, V = V, t = yrs, t0 = t0)
      q[i] <- det(topMat) / detDbot * sqrt(detInfoMle / Jrest[i])
    }
  }
  out$q <- q
  out$psi.max <- mle$par[2]
  out$r <- sign(out$psi.max - psi)*sqrt(-2*(mle$value+clogLik))
  out$rstar <- out$r + log(out$q/out$r)/out$r
  out$std.error <- sqrt(diag(solve(mle$hessian)))[2]
  out$psi <- psi
  class(out) <- "eprof"
  rstar_ori <- out$rstar
  # Usual blow up and numerical instability...
  out <-  mev::spline.corr(out)
  confint_array[t0i,,,1] <- confint(out)
  confint_array[t0i,,,2] <- confint(out, level = 0.99)
  if(t0i == 1L && figures){
  setwd("figures")
  fig <- "venice_profile.tex"
  tikz(fig, width = 6, height = 4, standAlone = TRUE)
    par(mar = c(4,4,1,1), bty = "l")
    plot(psi, -out$r^2/2, type = "l",
         ylab = "profile log-likelihood",
         xlab = "probability of exceedance of 194cm",
         panel.first = {abline(h=-qchisq(c(0.95,0.99),1)/2, col = scales::alpha("grey",0.8))},
         ylim = c(-4,0), xlim = c(0,0.025), xaxs= 'i')
    lines(psi, -out$rstar^2/2, col = 2, lty = 2)
    legend(x = "topright", col = 1:2, lty = 1:2, legend = c("profile","\\textsc{tem}"), bty = "n")
    dev.off()
    system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
    setwd("..")
  }
}

if(figures){
  setwd("figures")
  fig <- "venice_probexcyear.tex"
  tikz(fig, width = 7, height = 3, standAlone = TRUE)
}
yrv <- seq(2020, 2050, by = 5L)
conflimdf <- data.frame(
  year = rep(yrv, length.out = length(yrv)*2),
  estimator = factor(rep(rep(c("profile","\\textsc{tem}"), each = length(yrv)), length.out = length(yrv)*2)),
  value = c(confint_array[,1,,1]),
  lower95 = c(confint_array[,2,,1]),
  upper95 = c(confint_array[,3,,1]),
  lower99 = c(confint_array[,2,,2]),
  upper99 = c(confint_array[,3,,2])
)
conflimdf$estimator <- relevel(x = conflimdf$estimator, ref = "profile")
ggplot(data = conflimdf, aes(x = year, y = value, col = estimator)) +
  theme_classic() + theme(legend.position = "bottom") +
  scale_colour_grey() +
  geom_point(position = position_dodge2(width = 2, padding = 0.5)) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95),
                width = 2,
                position = position_dodge2(padding = 0.5)) +
  ylab("probability of exceedance")
  # geom_pointrange(aes(ymin = lower99, ymax = upper99),
  #                 fatten = 1,
  #                 position = position_dodge2(width = 2, padding = 0.5))
  ## Base R plot
  # par(mar = c(4,4,1,1), bty = "l")
  # plot(seq(2020, 2050, by = 5L) - 0.5,
  #         y = confint_array[,1,1,1],
  #         xlim= c(2019, 2051),
  #         ylim = c(0, 0.03), yaxs = "i",
  #         pch = 20,
  #         xlab = "year",
  #         ylab = "probability of exceedance of 194cm")
  # points(seq(2020, 2050, by = 5L)-0.5, confint_array[,2,1,1], pch = "-", cex = 2)
  # points(seq(2020, 2050, by = 5L)-0.5, confint_array[,3,1,1], pch = "-", cex = 2)
  # #TEM
  # points(seq(2020, 2050, by = 5L)+0.5, col = 2, confint_array[,1,2,1], pch = 20)
  # points(seq(2020, 2050, by = 5L)+0.5, col = 2, cex = 2, confint_array[,2,2,1], pch = "-")
  # points(seq(2020, 2050, by = 5L)+0.5, col = 2, cex = 2, confint_array[,3,2,1], pch = "-")
  # legend(x = "topleft", col = 1:2, pch = 20, legend = c("profile","\\textsc{tem}"), bty = "n")
  if(figures){
  dev.off()
  system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
  setwd("..")
  }

