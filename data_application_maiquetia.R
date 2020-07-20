setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
library(evd)
library(revdbayes)
library(mev)
library(xts)
library(lubridate)
library(ggplot2)
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

setTikzDefaults(overwrite = FALSE)
figures <- FALSE
dwidth <- 6
dheight <- 3.5
# Load data
data(maiquetia, package = "mev")

# Data application: GP distribution for POT - excluding 1999
dates <- seq.Date(from = as.Date("1961-01-01"), to = as.Date("1999-12-31"), by = "day")
maiq <- maiquetia[lubridate::year(dates) < 1999 & maiquetia > 0]
ratio <- length(maiq) / length(maiquetia[lubridate::year(dates) < 1999])

# Threshold selection diagnostics
mev::NC.diag(maiq, u = quantile(maiq, seq(0.95, 0.999, by = 0.002)))
abline(h = 0.05)
mev::W.diag(maiq, u = quantile(maiq, seq(0.95, 0.998, by = 0.005)))


threshf <- function(dat, p){
  dat <- na.omit(as.vector(dat))
  th <- quantile(dat, p[1])
 dat[dat > th] - th 
}

# Sensitivity of the conclusions to the choice of threshold
sapply(seq(0.96, 0.995, by = 0.0005), function(thp){ 
  c(quantile(maiq, thp), 
    gpd.mle(xdat = threshf(maiq, thp), args = c("shape","Nmean","Nquant"), 
          q = 0.5, N = 100*(1-thp)*ratio*365.25))
})

gpdf <- fit.gpd(maiq, threshold = (thresh <- 27), show = TRUE)

gpd.pll(param = "Nquant", dat = maiq, threshold = 30, q = 0.5,
                mod = c("profile","tem","modif"), N = 50*sum(maiq>=30)/length(maiq)*ratio*365.25, 
                psi = seq(85, 450, length = 201))
gpd.pll(param = "Nquant", dat = maiq, threshold = 35, q = 0.5,
        mod = c("profile","tem","modif"), N = 50*sum(maiq>=35)/length(maiq)*ratio*365.25, 
        psi = seq(90, 500, length = 201))
gpd.pll(param = "Nquant", dat = maiq, threshold = 40, q = 0.5,
        mod = c("profile","tem","modif"), N = 50*sum(maiq>=40)/length(maiq)*ratio*365.25, 
        psi = seq(95, 550, length = 201))
gpd.pll(param = "Nquant", dat = maiq, threshold = 45, q = 0.5,
        mod = c("profile","tem","modif"), N = 50*sum(maiq>=45)/length(maiq)*ratio*365.25, 
        psi = seq(100, 500, length = 201))
gpd.pll(param = "Nquant", dat = maiq, threshold = 50, q = 0.5,
        mod = c("profile","tem","modif"), N = 50*sum(maiq>=50)/length(maiq)*ratio*365.25, 
        psi = seq(100, 1200, length = 201))

if(figures){
  setwd("figures")
  fig <- "maiquetia_gpddiag.tex"
  tikz(file = fig, width = dwidth, height = dheight-0.5, standAlone = TRUE)
}
par(mar= c(4,4,1,1), mfrow = c(1,2))
plot(gpdf, which = 1, main = "", add =TRUE)
plot(gpdf, which = 2, main = "", add = TRUE)
if(figures){
  dev.off()
  system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
  setwd("..")
}

source('scripts/script_maiquetia.R')
# Second application use a cyclical kernel to weight the log-likelihood
ckern <- function(t, tj, nu, nday = 365){
  exp(nu*cos(2*pi*(t-tj)/nday))/exp(nu)
}

we <- ckern(t = yday(dates[which.max(maiquetia)]),
      tj = yday(dates[maiq > thresh]),
      nu = 1)
sum(we)


sub <- (dates < lubridate::dmy("01-01-1999")) & (maiquetia>0)
pthresh <- 0.96
ratio <- sum(sub) / length(maiquetia[(dates < lubridate::dmy("01-01-1999"))])
thresh <- quantile(maiquetia[sub], pthresh)
dat <- maiquetia[intersect(which(sub), which(maiquetia>thresh))] - thresh
time <- dates[intersect(which(sub), which(maiquetia>thresh))]
w1 <- ckern(t = yday(dates[which.max(maiquetia)]),
            tj = yday(time),
            nu = 2, nday = 365)
N <- 50*(1-pthresh)*ratio*365.25
par <- mev::fit.gpd(dat)$estimate
# Weighted likelihood - adapt all functions.
fitw1 <- gpdNw.pll(psi = seq(50, 1000, length.out = 100),
                  mod = c("profile","tem","modif"), N = N, dat = dat, w = w1)
confint(fitw1, print = TRUE)
w2 <- ckern(t = yday(dates[which.max(maiquetia)]),
                tj = yday(time),
                nu = 4, nday = 365)
fitw2 <- gpdNw.pll(psi = NULL,
                   mod = c("profile","tem","modif"), N = N, dat = dat, w = w2)
confint(fitw2, print = TRUE)
w3 <- ckern(t = yday(dates[which.max(maiquetia)]),
            tj = yday(time),
            nu = 7, nday = 365)
fitw3 <- gpdNw.pll(psi = NULL,
                   mod = c("profile","tem","modif"), N = N, dat = dat, w = w3)
confint(fitw3, print = TRUE)


if(figures){
  setwd("figures")
  fig <- "maiquetia_pot.tex"
  tikz(file = fig, width = dwidth, height =dheight, standAlone = TRUE)
}
par(mar = c(4,4,1,1))
plot(fitw2$pars[,1], y = fitw2$pll - fitw2$maxpll, ylab = "profile likelihood", lwd = 1.5,
     xlab = "Mean of semicentennial maximum distribution (in mm)", type = "l", bty = "l", panel.first = abline(h = -qchisq(c(0.95), df = 1)/2, lty = 2, col = "gray"), ylim = c(-4,0))
lines(fitw2$pars[,1], y = -fitw2$rstar^2/2, lty = 3, lwd = 1.5)
lines(fitw2$pars[,1], y = fitw2$tem.pll - fitw2$tem.maxpll, lty = 4, col = "gray70")
lines(fitw2$pars[,1], y = fitw2$empcov.pll - fitw2$empcov.maxpll, lty = 2, col = "gray70")
legend(x = "topright", bty = "n", lty = c(1,3,4,2), lwd = c(1.5,1.5,1,1), col = c(rep(1,2), rep("grey70",2)),
       legend = c("$\\ell_{\\mathsf{p}}$","$\\ell_{\\mathsf{fr}}$",
                  "$\\ell_{\\mathsf{m}}^{\\mathrm{cov}}$","$\\ell_{\\mathsf{m}}^{\\mathrm{tem}}$"))
rug(max(maiquetia))
if(figures){
  dev.off()
  system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
  setwd("..")
}