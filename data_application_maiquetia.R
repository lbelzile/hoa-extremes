setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('script_maiquetia.R')
setwd('..')
library(evd)
library(revdbayes)
library(mev)
library(xts)
library(lubridate)
library(ggplot2)
library(patchwork)
library(Rsolnp)
library(tikzDevice)
options(tikzLatexPackages =
          c("\\usepackage{tikz}\n",
            "\\usepackage[active,tightpage,psfixbb]{preview}\n",
            "\\usepackage{amsmath}",
            "\\PreviewEnvironment{pgfpicture}\n",
            "\\setlength\\PreviewBorder{0pt}\n",
            "\\usepackage{lmodern}\n",
            "\\usepackage{fourier}\n"
          )
)

setTikzDefaults(overwrite = FALSE)
figures <- FALSE
dwidth <- 6
dheight <- 3.5
# Load data
data(maiquetia, package = "mev")

# Data application: GP distribution for POT - excluding 1999
dates <- seq.Date(from = as.Date("1961-01-01"), to = as.Date("1999-12-31"), by = "day")
subset <- dates < as.Date("1999-01-01") & maiquetia > 0
maiq <- maiquetia[subset]
dates_sub <- dates[subset]
ql <- 0.5 #quantile level for N-year maximum
N <- 50L #N-year time horizon (50 years)
ratio <- length(maiq) / length(maiquetia[lubridate::year(dates) < 1999])

# Threshold selection diagnostics
# Northrop and Coleman
mev::NC.diag(maiq, u = quantile(maiq, seq(0.95, 0.999, by = 0.002)))
abline(h = 0.05)
# Wadsworth's diagnostic plots
useq <- quantile(maiq, seq(0.95, 0.998, by = 0.005))
Wdiag <- mev::W.diag(maiq, u = quantile(maiq, seq(0.95, 0.998, by = 0.005)))

# Extremogram
library(extremogram)
e1 <- extremogram1(maiq, quant = (qul <- ecdf(maiq)(27)), 
             maxlag = 15, start = 1, type = 1, ploting = 0)
B <- 9999
bootquant <- t(replicate(n = B, expr = extremogram:::permboot1(x = sample(maiq), quant = qul,
                          maxlag = 15, type = 1)[-1]))
qu_indep <- apply(bootquant, 2, quantile, probs =  0.95)
# Plot extremogram (not shown)
plot(1:14, e1[-1], ylim = c(0,1),
     yaxs = "i", bty = "l", type = "h",
     xlab = "Lag (in days)", ylab = "Extremogram", lwd = 2)
abline(h = quantile(bootquant, 0.95), col = scales::alpha(1, 0.5))

thresh <- u <- 27

# Fit a generalized Pareto distribution above 27mm via MLE
gpdf <- fit.gpd(maiq, threshold = u, show = TRUE)
# Calculate the profile likelihood for median of 50 year maximum
prof27u <- gpd.pll(param = "Nquant", dat = maiq, threshold = thresh, q = 0.5,
                   mod = c("profile","tem","modif"), N = 50*sum(maiq>=27)/length(maiq)*ratio*365.25, 
                   psi = seq(85, 450, length = 201), plot = FALSE)
# Confidence intervals
confint(prof27u, print = TRUE)

# Plots for Figure 1
if(figures){
  setwd("figures")
  fig <- "maiquetia_init.tex"
  tikz(file = fig, height = 3.5, width = 8, standAlone = TRUE)
}

par(mfrow = c(1,2), mar = c(5.5,4,1,1), bty = "l", pch = 20)
with(prof27u, 
     plot(x = psi, 
          y = pll - maxpll, 
          ylim = c(-4, 0),
          yaxs = "i",
          type = "l",
          ylab = "profile log likelihood", 
          xlab = "median semicentennial \nmaximum (in mm)",     
          yaxt = "n",
          panel.first = {
            abline(v = mle["Nquant"], col = scales::alpha(1, 0.1));
            abline(h = -qchisq(0.95,1)/2, col = scales::alpha(1, 0.3))
            abline(h = -qchisq(0.99,1)/2, col = scales::alpha(1, 0.3))})
)
Axis(x = c(-4,0), labels = paste0("$",-4:0,"$"), -4:0, side = 2)
rug(410.4)
# Plot of profile log likelihood
with(prof27u, lines(psi, -prof27u$rstar^2/2, lty = 2, col = "gray40"))
legend(x = "topright", 
       legend = c("$\\ell_{\\mathsf{p}}$","$\\ell_{\\mathsf{fr}}$"), 
       lty = c(1,2), 
       col = c("black","gray40"),
       bty = "n")

# Wadsworth' diagnostic plot
TradCI <- cbind(Wdiag$MLE[, 3] - qnorm(0.975) * sqrt(diag(Wdiag$Cov)), 
                Wdiag$MLE[, 3] + qnorm(0.975) * sqrt(diag(Wdiag$Cov)))
plot(useq, Wdiag$MLE[, 3], 
     ylim = c(-0.5,0.5), 
     xlab = "quantile",  
     bty = "l", 
     ylab = "shape",     
     panel.first = {abline(v= thresh, col ="gray10")},
     yaxt = "n")
Axis(x = seq(-0.4,0.4, by = 0.2), labels = paste0("$",seq(-0.4,0.4, by = 0.2),"$"), seq(-0.4,0.4, by = 0.2), side = 2)
lines(useq, TradCI[, 1], lty = 2)
lines(useq, TradCI[, 2], lty = 2)

if(figures){
  dev.off()
  system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
  setwd("..")
}

# Probability-probability and quantile-quantile plots for first fit
if(figures){
  setwd("figures")
  fig <- "maiquetia_gpddiag.tex"
  tikz(file = fig, width = dwidth, height = dheight-0.5, standAlone = TRUE)
}
par(mar= c(4,4,1,1), mfrow = c(1,2))
plot(gpdf, which = 1, main = "", add =TRUE, xlab = "theoretical quantiles", ylab = "sample quantiles")
plot(gpdf, which = 2, main = "", add = TRUE, xlab = "theoretical quantiles", ylab = "sample quantiles")
if(figures){
  dev.off()
  system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
  setwd("..")
}
#####################################################
###           Analysis of Section 4.1             ###
#####################################################
# Stopping rule
s <- 282
# All exceedances of the threshold, including 1999
yexc <- maiquetia[maiquetia>u][-sum(maiquetia > u)]
max51to60 <- c(154.0, 49.6, 46.7, 58.3, 70.5, 90.0, 70.1,
               105.7, 37.4, 40.8)
# Extract all yearly maximum for 1951-1999
allmax <- c(max51to60, as.numeric(apply.yearly(xts(maiquetia, order.by = dates), max)))
# Number of years of exceedances
ny <- as.numeric(round(diff(range(dates))/365.25))

# Compute the GEV maximum likelihood estimates 
# and change the parametrization to loc > median of N-year max
tpar <- gev.mle(xdat = allmax, args = c("Nquant","scale","shape"), N = 50, q = 0.5)
# Maximum likelihood estimate of standard likelihood (no stopping rule)
init_std <- optim(par = tpar,
                  fn = maiq_std.nll,
                  z = max51to60,
                  y = yexc,
                  u = u, 
                  ny = ny, 
                  N = 50, 
                  ql = 0.5,
                  method = "Nelder-Mead")
# Maximum likelihood estimates with full conditioning
init_fc <- optim(par = tpar,
                 fn = maiq_fc.nll,
                 z = max51to60,
                 y = yexc,
                 s = 282,
                 u = u, 
                 ny = ny, 
                 N = N, 
                 ql = ql,
                 method = "Nelder-Mead")
# Profile log likelihood for standard likelihood (no stopping rule)
psi_std <- sort(c(seq(130, 900, length.out = 771), 410.4))
prof_std <- profile_qp_std_maiq(psi = psi_std,
                                start = init_std$par,
                                z = max51to60,
                                y = yexc,
                                ny = ny, 
                                N = N, 
                                ql = ql,
                                mod = c("profile","tem"), 
                                threshold = u,
                                plot = FALSE)
# Profile full conditional log likelihood
psi_fc <- sort(c(seq(130, 800, length.out = 301), 410.4))
prof_fc <- profile_qp_fc_maiq(psi = psi_fc,
                              start = init_fc$par,
                              z = max51to60,
                              y = yexc,
                              ny = ny, 
                              N = N, 
                              ql = ql,
                              s = 282, 
                              mod = c("profile","tem"),
                              threshold = u,
                              plot = FALSE
                              )
# Compute confidence intervals (use smoothing spline to interpolate between values)
confint_std <- confint(prof_std, method = "smooth.spline", print = TRUE)
confint_fc <- confint(prof_fc, method = "smooth.spline", print = TRUE)

#' Function to transform GEV parametrization from N-year \code{ql} quantile to location parameter
mle_qp2mu <- function(par, N = 50, ql = 0.5){
  par <- c(par[1]+par[2]/par[3]*(1-(-N/log(ql))^par[3]), par[2], par[3])
  names(par) <- c("loc","scale","shape")
  par
}
# Compute bootstrap confidence bands for P-P plot
set.seed(2134)

# Probability-probability plot
# We have three kinds of data: left-truncated exceedances,
# right-truncated exceedances and yearly maxima
mle <- mle_qp2mu(init_fc$par)
Lambda_f <- function(x, par = mle, ny){
  mu = par[1]; sigma = par[2]; xi = par[3]
  ny*pmax(1+xi*(x-mu)/sigma,0)^(-1/xi)
}
# Transform maxima to uniform scale
punif <- c(evd::pgev(q=max51to60, loc = mle[1], scale = mle[2], shape = mle[3]),
           (Lambda_f(u, mle, ny = ny) - Lambda_f(yexc[-length(yexc)], mle, ny = ny))/(Lambda_f(u, mle, ny = ny)-Lambda_f(s, mle, ny = ny)),
           (Lambda_f(s, mle, ny = ny) - Lambda_f(yexc[length(yexc)], mle, ny = ny))/(Lambda_f(s, mle, ny = ny)))
# Compute envelope
envel <- boot::envelope(mat = t(replicate(sort(runif(length(punif))), n = 1e5)))
if(figures){
  setwd("figures")
  fig <- "profile_hoa_maiquetia.tex"
  tikz(file = fig, height = 0.8*dheight, width = dwidth, standAlone = TRUE)
}
par(mfrow = c(1,2), bty = "l", mar = c(4,4,1,1))
plot(prof_std$psi, -prof_std$rstar^2/2, lwd = 2,
     lty = 2, col = "gray70", type = "l", ylim = c(-3,0.05),
     yaxs = "i", ylab = "profile log likelihood",
     xlab = "median semicentenial maximum",
     yaxt = "n",
     panel.first = {abline(h = -qchisq(0.95,1)/2, col = "grey")})
lines(prof_std$psi, -prof_std$r^2/2, lwd = 2,
      lty = 1, col = "gray70")
lines(prof_fc$psi, -prof_fc$rstar^2/2, lwd = 2,
      lty = 2)
lines(prof_fc$psi, -prof_fc$r^2/2, lwd = 2,
      lty = 1)
Axis(x = c(-3,0), labels = paste0("$",-3:0,"$"), -3:0, side = 2)
rug(410.4)
legend(x = "topright", 
       legend = c("$\\ell_{\\mathrm{fc},\\mathsf{p}}$",
                  "$\\ell_{\\mathrm{fc},\\mathsf{fr}}$",
                  "$\\ell_{\\mathrm{std},\\mathsf{p}}$",
                  "$\\ell_{\\mathrm{std},\\mathsf{fr}}$"), 
       lty = rep(c(1,2), length.out = 4),
       lwd = 2,
       col = rep(c("black","gray70"), each = 2),
       bty = "n")

# Probability-probability plot
plot(y=sort(punif), ppoints(length(punif)),
     xlab = "theorical quantiles", ylab = "empirical quantiles",
     col = scales::alpha("black", 0.8),
     pch = 20, bty = "l", yaxs = "i", xaxs = "i",
     ylim = c(0,1),
     xlim = c(0,1),
     panel.first = {abline(a=0,b=1)})

lines(x=ppoints(length(punif)),envel$overall[1,], lty = 2)
lines(x=ppoints(length(punif)),envel$overall[2,], lty = 2)
if(figures){
  dev.off()
  system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
  setwd("..")
}

# Significance level for test median of 50 years > 410.4
with(prof_fc, pnorm(r[psi == 410.4], lower.tail = TRUE))
with(prof_fc, pnorm(rstar[psi == 410.4], lower.tail = TRUE))

with(prof_std, pnorm(r[psi == 410.4], lower.tail = TRUE))
with(prof_std, pnorm(rstar[psi == 410.4], lower.tail = TRUE))

