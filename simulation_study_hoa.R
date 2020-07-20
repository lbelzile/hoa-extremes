setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(simsalapar)
library(parallel)
library(flexsurv)
library(actuar)
library(evd)
library(mev)
library(ggplot2)
library(viridis)
library(xtable)
library(tikzDevice)
options(tikzLatexPackages =
          c("\\usepackage{tikz}\n",
            "\\usepackage[active,tightpage,psfixbb]{preview}\n",
            "\\usepackage{amsmath}",
            "\\PreviewEnvironment{pgfpicture}\n",
            "\\setlength\\PreviewBorder{0pt}\n",
            "\\usepackage{fourier}\n",
            "\\DeclareMathAlphabet{\\mathsf}{OT1}{pag}{m}{n}\n"
          )
)
setTikzDefaults(overwrite = FALSE)
figures <- FALSE
save <- FALSE
# Function to obtain the penultimate approximations for target
obtain_true_vals <- function(Nn = 9000, p = 1/Nn){
  true_vals <- array(dim = c(12,3), 
                     dimnames = list(dist = c("burr", "weibull", "gengamma.orig", "norm", "lnorm", "stud",
                                              "gevp", "gevz", "gevn", "gpdp", "exp", "gpdn"),
                                     parameter = c("quant", "Nquant", "Nmean")))
  for(dist_ind  in 1:12){
    dist <- c("burr", "weibull", "gengamma.orig", "norm", "lnorm", "t", rep("gev", 3), "gpd", "exp", "gpd")[dist_ind]
    dist_param <- switch(dist_ind,
                         list(shape1 = 2, shape2 = 5), #burr
                         list(scale = 1, shape = 2 / 3), #weibull
                         list(shape = 0.54, scale = 1.83, k = 1.11 / 0.54), # gengamma
                         list(mean = 0, sd = 1), #norm
                         list(meanlog = 0, sdlog = 1), #lnorm
                         list(df = 10, ncp = 0), #student
                         list(loc = 0, scale = 1, shape = 0.1), #gev
                         list(loc = 0, scale = 1, shape = 0), #gev
                         list(loc = 0, scale = 1, shape = -0.1), #gev
                         list(loc = 0, scale = 1, shape = 0.1), #gpd, positive shape
                         list(rate = 1), #exponential
                         list(loc = 0, scale = 1, shape = -0.1) #gpd, neg shape
    )
    penult <- do.call(what = "smith.penult", args = c(family = dist, model = "bm", m = Nn, dist_param))
    #print(penult$shape)
    #Create grid of values of psi, includes the quantity of interest for the null,
    #as well as nearby values to check power of test
    true_vals[dist_ind, ] <- c(do.call(paste0("q", dist), args = c(dist_param, p = 1-p)),
                               Nquant = qgev(0.5, loc = penult$loc, scale = penult$scale, shape = penult$shape),
                               Nmean = gevN.mean(c(penult$loc, penult$scale, penult$shape), N = 1))
  }
  true_vals
}
true_values <- obtain_true_vals()

# Simulation settings
npy <- 90L
ny <- 20L
Ntot <- 100L

# Function that computes the results for each simulation setting
doOne <-  function(n, dist_ind, level, N, q, p, parameter, method, type, vari, clvind) {
  res <- array(dim = c(3, 2, 3, 6, 5, 3),
               dimnames =  list(
                 "parameter" = c("quant", "Nquant", "Nmean"),
                 "method" = c("bm", "pot"),
                 "level" = c("l1","l2","l3"),
                 "vari" = c("PE", "LoCI", "UpCI", "RB", "CP", "CIW"),
                 "type" = c("wald","profile", "tem", "modif.tem", "modif.empcov"),
                 "covlev" = c(0.1,0.05,0.01)
               )
  )
  dist <- c("burr", "weibull", "gengamma.orig", "norm", "lnorm", "t", rep("gev", 3), "gpd", "exp", "gpd")[dist_ind]
  #Create sample
  dist_param <- switch(dist_ind,
                       list(shape1 = 2, shape2 = 5), #burr
                       list(scale = 1, shape = 2 / 3), #weibull
                       list(shape = 0.54, scale = 1.83, k = 1.11 / 0.54), # gengamma
                       list(mean = 0, sd = 1), #norm
                       list(meanlog = 0, sdlog = 1), #lnorm
                       list(df = 10, ncp = 0), #student
                       list(loc = 0, scale = 1, shape = 0.1), #gev
                       list(loc = 0, scale = 1, shape = 0), #gev
                       list(loc = 0, scale = 1, shape = -0.1), #gev
                       list(loc = 0, scale = 1, shape = 0.1), #gpd, positive shape
                       list(rate = 1), #exponential
                       list(loc = 0, scale = 1, shape = -0.1) #gpd, neg shape
  )
  dat <- do.call(paste0("r", dist), args = c(n = n, dist_param))
  #Create grid of values of psi, includes the quantity of interest for the null,
  #as well as nearby values to check power of test
  psi <- c(do.call(paste0("q", dist), args = c(p = 0.99, dist_param)),
           seq(do.call(paste0("q", dist), args = c(p = 0.999, dist_param)), 
               do.call(paste0("q", dist), args = c(p = 0.9999, dist_param)), length = 15),
           seq(do.call(paste0("q", dist), args = c(p = 0.99991, dist_param)),
               do.call(paste0("q", dist), args = c(p = 0.9999999999, dist_param)), length = 30))
  for (para in parameter) {
    #For threshold exceedances, compute the excesses over threshold at level pu
    for (lev in level) {
      u <- as.vector(quantile(dat, 1-switch(lev, 60/n, 40/n, 20/n))) #sample size 20, 40, 60
      exceedances <- dat[dat > u] - u
      #For block maxima with sample of size n at least 1000
      m <- switch(lev, 30, 45, 90) #Sample size n/30, n/45, n/90
      blm <- apply(matrix(dat, ncol = m), 1, max)
      #Run the profile log likelihood routine
      for (meto in method) {
        #we target the distribution of the maxima of N * n observations
        profile <- switch(meto,
                          bm = try(gev.pll(psi = psi, param = para, dat = blm, N = N * n / m, q = q, p = p*m, mod = c("tem","modif"))),
                          #divide by number of blocks to get number of block maxima samples
                          pot = try(gpd.pll(psi = psi[(psi-u)>0], param = para, dat = dat, mod = c("tem","modif"), N = switch(lev, 60, 40, 20) * N, 
                                            q = q, p = p/switch(lev, 60/n, 40/n, 20/n), threshold = u))
                          #for GPD, this is the average maximum of N exceedances
        )
        if (is.character(profile)) { next; }
        ind1 <- switch(para,
                       quant = 1,
                       Nquant = 2,
                       Nmean = 3)
        ind2 <- switch(meto, bm = 1, pot = 2)
        conf_inter <- array(NA, dim = c(3,5,3))
        for(clv in clvind){
          conf_inter_temp <- try(confint(profile, level = c(0.9,0.95,0.99)[clv]))
          if (!is.character(conf_inter)) {
            conf_inter[,2:5,clv] <- conf_inter_temp
          }
          wald <- try( c(profile$psi.max, 
                         profile$psi.max - qnorm(c(0.95,0.975,0.995)[clv])*profile$std.error, 
                         profile$psi.max + qnorm(c(0.95,0.975,0.995)[clv])*profile$std.error))
          if(!is.character(wald)){
            conf_inter[,1,clv] <- wald 
          }
          if(meto == "pot"){ #FIX Must add back the threshold to relate to the original data!
            conf_inter[,,clv] <- conf_inter[,,clv] + u
          }
          res[ind1, ind2, lev, , , clv] <-
            rbind(conf_inter[,,clv], 
                  (conf_inter[1, , 1] - true_values[dist_ind, ind1]) / true_values[dist_ind, ind1],  # relative bias
                  I(conf_inter[2,, clv] <= true_values[dist_ind, ind1]) * I(conf_inter[3,, clv] >= true_values[dist_ind, ind1]),  
                  # coverage proba (binary) FIXED
                  conf_inter[3,, clv]-conf_inter[2,, clv]  # confidence interval width
            )
        }
      }
    }
  }
  res
}

# Number of simulation replications
Nsim <- 13000L
#Set list of variables for the simulation study
varListTEM <- varlist(
  n.sim = list(type = "N", expr = quote(N), value = Nsim),
  # replications
  n = list(type = "frozen", value = c(npy*ny)),
  # sample size
  dist_ind = list(type = "grid",  value = 1:12),
  N = list(type = "frozen", value = c(Ntot/ny)),
  q = list(type = "frozen", value = c(0.5)),
  p = list(type = "frozen", value = c(1/(npy*Ntot))),
  parameter = list(type = "inner", value = c("quant", "Nquant", "Nmean")),
  method = list(type = "inner", value = c("bm", "pot")),
  level = list(type = "inner", value =  c("l1","l2","l3")),
  type = list(type = "inner", value = c("wald", "profile", "tem", "modif.tem", "modif.empcov")),
  vari = list(type = "inner", value = c("PE", "LowCI", "UpCI", "RB", "CP", "CIW")),
  clvind = list(type = "inner", value = c(0.1,0.05,0.01))
)
if(save){
file <- simsalapar::doLapply(vList = varListTEM, doOne = doOne, repFirst = TRUE)
test <- mkAL(file, vList = varListTEM, repFirst = TRUE)
dimname <- dimnames(test[[1]]$value)
for(i in 1:length(test)){
   dimnames(test[[i]]$value) <- dimname
}
valTEM <- getArray(test, comp = "value")
rm(test)

# Transform the Wald intervals
for(i1 in 1:3){#parameter
  for(i2 in 1:2){ #method (BM vs POT)
    for(i3 in 1:3){ #level of threshold
      for(i4 in 1:3){ 
        #for(i4 in 1:5){#estimator
        for(i5 in 1:12){ #distind      
          stderrWald <- (valTEM[i1,i2,i3,"UpCI","wald",i4,i5,] - valTEM[i1,i2,i3,"PE","wald",i4,i5,])/qnorm(1-switch(i4, 0.1, 0.05,0.01)/2)
          valTEM[i1,i2,i3,"UpCI","wald",i4,i5,] <- exp(log(valTEM[i1,i2,i3,"PE","wald",i4,i5,]) + qnorm(1-switch(i4, 0.1, 0.05,0.01)/2)*stderrWald/valTEM[i1,i2,i3,"PE","wald",i4,i5,])
          valTEM[i1,i2,i3,"LoCI","wald",i4,i5,] <- exp(log(valTEM[i1,i2,i3,"PE","wald",i4,i5,]) - qnorm(1-switch(i4, 0.1, 0.05,0.01)/2)*stderrWald/valTEM[i1,i2,i3,"PE","wald",i4,i5,])
          valTEM[i1,i2,i3,"CIW","wald",i4,i5,] <- valTEM[i1,i2,i3,"UpCI","wald",i4,i5,] - valTEM[i1,i2,i3,"LoCI","wald",i4,i5,]
          valTEM[i1,i2,i3,"CP","wald",i4,i5,] <- I(valTEM[i1,i2,i3,"UpCI","wald",i4,i5,]> true_values[i5,i1]) *  I(valTEM[i1,i2,i3,"LoCI","wald",i4,i5,] <  true_values[i5,i1]) 
        }
      }
    }
  }
}

#Transform confidence interval width to relative width (relative to profile)
for(j0 in c(1,3:5)){
  valTEM[,,,"CIW",j0,,,] <- valTEM[,,,"CIW",j0,,,]/valTEM[,,,"CIW",2,,,]
}




# names(dimnames(valTEM))
# str(valTEM)
#dim 1: parameter: quant, Nquant, Nmean
#dim 2: method: bm, pot
#dim 3: level: (1, 2, 3)
#dim 4: type: PE, LoCI, UpCI, RB, CP, CIW
#dim 5: vari: wald, profile, tem, modif.tem, modif.empcov
#dim 6: clvind: 0.1, 0.05, 0.01
#dim 7: dist_ind  1:12
#dim 8: Nsim

#Remove $ from the list of escape characters
MyescapeLatex <- function(x) {
   x <- gsub("| ", "\\textbar\\ ", x, fixed = TRUE)
   x <- gsub("|", "\\textbar", x, fixed = TRUE)
  x <- gsub("#", "\\#", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x <- gsub("~ ", "\\textasciitilde\\ ", x, fixed = TRUE)
  x <- gsub("~", "\\textasciitilde ", x, fixed = TRUE)
  x <- gsub("^", "\\verb|^|", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  #x <- gsub("{", "\\{", x, fixed = TRUE)
  #x <- gsub("}", "\\}", x, fixed = TRUE)
  x <- gsub(" '", " `", x, fixed = TRUE)
  x <- gsub(" \"", " ``", x, fixed = TRUE)
  x <- gsub("... ", "\\dots\\ ", x, fixed = TRUE)
  x <- gsub("...",  "\\dots", x, fixed = TRUE)
  x <- gsub(" - ", " -- ", x, fixed = TRUE)
  x
}
#Bias and standard error
# bias <- apply(valTEM[,,,"RB",-1,,,], 1:6, function(x){ mean(x, na.rm=TRUE)})[,,,,"0.1",]
# sd.bias <- apply(valTEM[,,,"RB",-1,,,], 1:6, function(x){ sd(x, na.rm=TRUE)/sqrt(length(x))})[,,,,"0.1",]
# for(dist_ind in 1:12){
#   for(qty_ind in 1:3){
#     sd.bias[qty_ind,,,,dist_ind] <- sd.bias[qty_ind,,,, dist_ind]/true_values[dist_ind, qty_ind]
#   }
# }


nsim <- dim(valTEM)['n.sim']

nomerror <- array(0, dim = c(6,3,2,3,5,12))
for(i1 in 1:3){#parameter
 for(i2 in 1:2){ #method (BM vs POT)
   for(i3 in 1:3){ #level of threshold
     for(i4 in 1:5){#estimator
       for(i5 in 1:12){ #distind      
       nomerror[1,i1,i2,i3,i4,i5] <- 100*mean(I(valTEM[i1,i2, i3,"LoCI",i4,"0.01",i5,] > true_values[i5,i1]), na.rm = TRUE)
       nomerror[2,i1,i2,i3,i4,i5] <- 100*mean(I(valTEM[i1,i2, i3,"LoCI",i4,"0.05",i5,] > true_values[i5,i1]), na.rm = TRUE)
       nomerror[3,i1,i2,i3,i4,i5] <- 100*mean(I(valTEM[i1,i2, i3,"LoCI",i4,"0.1", i5,] > true_values[i5,i1]), na.rm = TRUE)
       nomerror[4,i1,i2,i3,i4,i5] <- 100*mean(I(valTEM[i1,i2, i3,"UpCI",i4,"0.1", i5,] < true_values[i5,i1]), na.rm = TRUE)
       nomerror[5,i1,i2,i3,i4,i5] <- 100*mean(I(valTEM[i1,i2, i3,"UpCI",i4,"0.05",i5,] < true_values[i5,i1]), na.rm = TRUE)
       nomerror[6,i1,i2,i3,i4,i5] <- 100*mean(I(valTEM[i1,i2, i3,"UpCI",i4,"0.01",i5,] < true_values[i5,i1]), na.rm = TRUE)
      }   
   }
  }
 }
}
dimnames(nomerror) <- c( list(rate = c(0.5,2.5,5,95,97.5,99.5)), dimnames(valTEM)[c(1,2,3,5,7)])
save(nomerror, file = "TEM_nominalerror.RData")
} else{
  load("TEM_nominalerror.RData") 
}

# Produce tables

setwd("..")
setwd("tables")
#Table Penultimate approximation parameters
Nn <- 9000L
dist <- c("burr", "weibull", "gengamma.orig", "norm", "lnorm", "t", rep("gev", 3), "gpd", "exp", "gpd")
penult <- t(sapply(1:6, function(dist_ind){
  distF <- dist[dist_ind]
  dist_param <- switch(dist_ind,
                       list(shape1 = 2, shape2 = 5), #burr
                       list(scale = 1, shape = 2 / 3), #weibull
                       list(shape = 0.54, scale = 1.83, k = 1.11 / 0.54), # gengamma
                       list(mean = 0, sd = 1), #norm
                       list(meanlog = 0, sdlog = 1), #lnorm
                       list(df = 10, ncp = 0), #student
                       list(loc = 0, scale = 1, shape = 0.1), #gev
                       list(loc = 0, scale = 1, shape = 0), #gev
                       list(loc = 0, scale = 1, shape = -0.1), #gev
                       list(loc = 0, scale = 1, shape = 0.1), #gpd, positive shape
                       list(rate = 1), #exponential
                       list(loc = 0, scale = 1, shape = -0.1) #gpd, neg shape
  )
  c(sapply(c(1-c(60,40,20)/(20*90)), function(qu){do.call(what = "smith.penult", args = c(family = distF, model = "pot", 
                                                                                          qu = qu, dist_param))$shape}),
    sapply(c(30, 45, 90, Nn), function(N){
      do.call(what = "smith.penult", args = c(family = distF, model = "bm", m = N, dist_param))$shape
      
    }))
}))

penult <- cbind(penult, c(1/10,0,0,0,0,1/10))
dimnames(penult) <- list(Distribution = c("Burr($c=2, k=5$)", "Weibull($\\lambda=1,\\alpha=2/3$)", "Gen. gamma($\\beta=1.83, \\gamma_1=1.11, \\gamma_2=0.54$)", "Normal", "Lognormal", "Student($\\nu=10$)"),
                         "Approx." = c(paste0("POT, $q=", round(1-c(60,40,20)/(20*90),3),"$"), 
                                       paste0("$\\xi_{",c(30,45,90,Nn),"}$"), "$\\xi_{\\infty}$"))

dimnames(penult) <- list(Distribution = c("Burr", "Weibull", "Gen. gamma", "Gaussian", "Lognormal", "Student"),
                         " " = c(paste0("\\small{$q=", round(1-c(60,40,20)/(20*90),3),"$}"), 
                                       paste0("\\small{$m = ",c(30,45,90,Nn),"$}"), "$\\xi_{\\infty}$"))
ftabu <- ftable(apply(penult, 1:2, function(x){paste0("$",formatC(x, format='f', digits=2),"$")}))
ftabu[,8] <- c("$0.1$", rep("$0$", 4), "$0.1$")
tabu <- toLatex(ftabu, 
                escapeFUN =   identity, booktabs = TRUE, align = "*{1}{l}*{7}{r}*{1}{l}",
                fontsize = "small", center = TRUE, caption = "Penultimate shape parameters for six distributions, based on threshold exceedances with threshold at $q$ percentile (first three columns), block maxima with maximum of $m$ observations (fourth to sixth column).  The penultimate shape parameter for the maximum of $9000$ observations is the reference, still far from the tail index $\\xi_{\\infty}$ in the last column.", label = "table:penult_shape", placement  = "t!")  
cat(tabu, file = "Table_penult_shape_TEMsimu.tex", append = FALSE)


# Produce tables of coverage probability and confidence interval width
for(meto in c("CP","CIW")){
  obj <- switch(meto,
                CP = apply(valTEM[,,, "CP",,,,], 1:6, function(x){ mean(x, na.rm = TRUE)})*100,
                RB = apply(valTEM[,,,"RB",-1,,,], 1:6, function(x){ mean(x, na.rm=TRUE, trim = 0.1)})*100,
                CIW = apply(valTEM[,,,"CIW",-2,,,], 1:6, function(x){ mean(x, na.rm=TRUE, trim = 0.1)})*100)
 
  for(met in c("bm","pot")){
    for(levl in 1:3){
      tabu <- obj[,met,paste0("l",levl),,,-c(switch(met, bm = 10:12, pot = 7:9))]
      if(met == "pot"){
        dimnames(tabu) <- list( Parameter = c("Quantile", "$N$-obs. median", "$N$-obs. mean"),
                                Method = switch(meto,
                                                CP = c("Wald", "profile", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                                                CIW = c("Wald", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                                                RB =  c("\\textsc{mle}", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)")),
                                "Conf. level (%)" = c(90,95,99), 
                                "$F$" = paste0("$F_",c(1:9),"$"))
                                # Distribution = c("Burr", "Weibull", "Gen. gamma", 
                                #                  "Normal", "Log-normal", "Student ", "$\\mathsf{GP}(\\xi =0.1)$",
                                #                  "Exponential","$\\mathsf{GP}(\\xi=-0.1)$"))
      } else {
        dimnames(tabu) <- list( Parameter = c("Quantile", "$N$-obs. median", "$N$-obs. mean"),
                                Method = switch(meto,
                                                CP = c("Wald", "profile", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                                                CIW = c("Wald", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                                                RB =  c("\\textsc{mle}", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)")),
                                "Conf. level (%)" = c(90,95,99), 
                                "$F$" = paste0("$F_",c(1:9),"$"))
                                # Distribution = c("Burr", "Weibull", "Gen. gamma", 
                                #                  "Normal", "Log-normal", "Student", "$\\mathsf{GEV}(\\xi=0.1)$",
                                #                  "$\\mathsf{GEV}(\\xi = 0)$","$\\mathsf{GEV}(\\xi=-0.1)$"))
        
      }
      
      ndigits <- switch(meto, CP = 1, CIW = 0, RB = 0)
      ftab1 <- ftable(apply(round(tabu, ndigits), 1:4, function(x){  paste0("$",formatC(x, format=switch(meto, CP = 'f', CIW = 'f', RB = "d"), digits=ndigits),"$")}),
                      col.vars = c(1, 3), row.vars = c(4, 2))
      #ftab1 <- ftable(round(tabu,1), col.vars = c(1, 3), row.vars = c(4, 2))
      str <- switch(met, bm = paste0("with $m=", switch(levl, 30, 45, 90), ", k=", 20*90/switch(levl, 30, 45, 90), "$."),
                    pot = paste0("with $k=", switch(levl, 60,40,20),"$."))
      #% not escaped in caption by escapeFUN
      fftab1 <- toLatex(ftab1, booktabs = TRUE, fontsize = "small", escapeFUN = MyescapeLatex, 
                        caption = paste0(switch(meto,
                                                CP ="Coverage probability (\\%), ", 
                                                CIW = "Truncated mean ($\\alpha=0.1$) of the ratio of the confidence interval width relative to the width of profile confidence interval (in \\%), ", 
                                                RB = "Truncated mean ($\\alpha=0.1$) of relative bias (in \\%), "), 
                                         switch(met, pot = "peaks-over-threshold", bm = "block maximum"), " method ", str, 
                          " The largest standard error",
switch(meto, CP = "", CIW = ", obtained using a nonparametric bootstrap, ", RB = ", obtained using a nonparametric bootstrap, "), "is ", round(max(switch(meto, CP = bootstdCP, RB = bootstdRB, CIW = bootstdCIW)[,switch(met,pot = 2, bm = 1),levl,,,]),2), "\\%.",
" The distributions (from top to bottom) are Burr ($F_1$), Weibull ($F_2$), generalized gamma ($F_3$), normal ($F_4$), lognormal ($F_5$), Student $t$ ($F_6$), ", 
                                     switch(met, 
                                               bm = "$\\mathsf{GEV}(\\xi =0.1)$ ($F_7$), Gumbel ($F_8$) and $\\mathsf{GEV}(\\xi=-0.1)$ ($F_9$)",
                                               pot = "$\\mathsf{GP}(\\xi =0.1)$ ($F_7$), exponential ($F_8$) and $\\mathsf{GP}(\\xi=-0.1)$ ($F_9$)"),"."),
                        label = paste0("tab_",meto,"_l",levl,"_",met), align = "*{2}{l}*{9}{r}", center = TRUE, placement  = "t!")
      cat(fftab1, file = paste0("Table_tem_", met,"_",
                                switch(met, 
                                       pot = switch(levl, 60,40,20), 
                                       bm = paste0("m",switch(levl, 30, 45, 90)))
                                , "_", meto,".tex"))
    }
  }
}

# Produce table of relative bias
for(meto in c("RB")){
  obj <- apply(valTEM[,,,"RB",-1,,,], 1:6, function(x){ mean(x, na.rm=TRUE, trim = 0.1)})*100
  
  for(met in c("bm","pot")){
    for(levl in 1:3){
      tabu <- obj[,met,paste0("l",levl),,,-c(switch(met, bm = 10:12, pot = 7:9))]
      if(met == "pot"){
        dimnames(tabu) <- list( Parameter = c("Quantile", "$N$-obs. median", "$N$-obs. mean"),
                                Method = switch(meto,
                                                CP = c("Wald", "profile", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                                                CIW = c("Wald", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                                                RB =  c("\\textsc{mle}", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)")),
                                "Conf. level (%)" = c(90,95,99), 
                                "$F$" = paste0("$F_",c(1:9),"$"))
        # Distribution = c("Burr", "Weibull", "Gen. gamma", 
        #                  "Normal", "Log-normal", "Student ", "$\\mathsf{GP}(\\xi =0.1)$",
        #                  "Exponential","$\\mathsf{GP}(\\xi=-0.1)$"))
      } else {
        dimnames(tabu) <- list( Parameter = c("Quantile", "$N$-obs. median", "$N$-obs. mean"),
                                Method = switch(meto,
                                                CP = c("Wald", "profile", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                                                CIW = c("Wald", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                                                RB =  c("\\textsc{mle}", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)")),
                                "Conf. level (%)" = c(90,95,99), 
                                "$F$" = paste0("$F_",c(1:9),"$"))
        # Distribution = c("Burr", "Weibull", "Gen. gamma", 
        #                  "Normal", "Log-normal", "Student", "$\\mathsf{GEV}(\\xi=0.1)$",
        #                  "$\\mathsf{GEV}(\\xi = 0)$","$\\mathsf{GEV}(\\xi=-0.1)$"))
        
      }
      
      ndigits <- switch(meto, CP = 1, CIW = 0, RB = 0)
      ftab1 <- ftable(apply(round(tabu[,,1,], ndigits), 1:3, function(x){  paste0("$",formatC(x, format=switch(meto, CP = 'f', CIW = 'f', RB = "d"), digits=ndigits),"$")}),
                      col.vars = 3, row.vars = c(1,2))
      #ftab1 <- ftable(round(tabu,1), col.vars = c(1, 3), row.vars = c(4, 2))
      str <- switch(met, bm = paste0("with $m=", switch(levl, 30, 45, 90), ", k=", 20*90/switch(levl, 30, 45, 90), "$."),
                    pot = paste0("with $k=", switch(levl, 60,40,20),"$."))
      #% not escaped in caption by escapeFUN
      fftab1 <- toLatex(ftab1, booktabs = TRUE, fontsize = "small", escapeFUN = MyescapeLatex, 
                        caption = paste0(switch(meto,
                                                CP ="Coverage probability (\\%), ", 
                                                CIW = "Truncated mean ($\\alpha=0.1$) of the ratio of the confidence interval width relative to the width of profile confidence interval (in \\%), ", 
                                                RB = "Truncated mean ($\\alpha=0.1$) of relative bias (in \\%), "), 
                                         switch(met, pot = "peaks-over-threshold", bm = "block maximum"), " method ", str, 
                                         " The largest standard error",
                                         switch(meto, CP = "", CIW = ", obtained using a nonparametric bootstrap, ", RB = ", obtained using a nonparametric bootstrap, "), "is ", round(max(switch(meto, CP = bootstdCP, RB = bootstdRB, CIW = bootstdCIW)[,switch(met,pot = 2, bm = 1),levl,,,]),2), "\\%.",
                                         " The distributions (from top to bottom) are Burr ($F_1$), Weibull ($F_2$), generalized gamma ($F_3$), normal ($F_4$), lognormal ($F_5$), Student $t$ ($F_6$), ", 
                                         switch(met, 
                                                bm = "$\\mathsf{GEV}(\\xi =0.1)$ ($F_7$), Gumbel ($F_8$) and $\\mathsf{GEV}(\\xi=-0.1)$ ($F_9$)",
                                                pot = "$\\mathsf{GP}(\\xi =0.1)$ ($F_7$), exponential ($F_8$) and $\\mathsf{GP}(\\xi=-0.1)$ ($F_9$)"),"."),
                        label = paste0("tab_",meto,"_l",levl,"_",met), align = "*{2}{l}*{9}{r}", center = TRUE, placement  = "t!")
      cat(fftab1, file = paste0("Table_tem_", met,"_",
                                switch(met, 
                                       pot = switch(levl, 60,40,20), 
                                       bm = paste0("m",switch(levl, 30, 45, 90)))
                                , "_", meto,".tex"))
    }
  }
}

# Produce tables of empirical error rates
for(met in c("pot", "bm")){
  for(levl in 1:3){
    if(met == "pot"){
      tabu <- nomerror[,c(1,3),met,levl,,c(1:6,10:12)]
      dimnames(tabu) <- list("Error rate" = c("0.5","2.5","5", "5","2.5","0.5"),
                             Parameter = c("Quantile", "$N$-obs. mean"),
                             Method = c("Wald", "profile", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                             "$F$" = paste0("$F_",c(1:9),"$") 
                             #"Burr", "Weibull", "Gen. gamma", "Normal", "Log-normal", "Student ", "GP(0.1)", "Exponential","GP(-0.1)")
      )
    } else if(met == "bm"){
      tabu <- nomerror[,c(1,3),met,levl,,c(1:9)]
      dimnames(tabu) <- list("Error rate" = c("0.5","2.5","5", "5","2.5","0.5"),
                             Parameter = c("Quantile", "$N$-obs. mean"),
                             Method = c("Wald", "profile", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                             "$F$" = paste0("$F_",c(1:9),"$"))
      
    }
    # Tables
    str <- switch(met, bm = paste0("with $m=", switch(levl, 30, 45, 90), ", k=", 20*90/switch(levl, 30, 45, 90), "$."),
                  pot = paste0("with $k=", switch(levl, 60,40,20),"$."))
    ftabu <- ftable(apply(round(2*tabu,0)/2, 1:4, function(x){paste0("$",formatC(x, format='f', digits=1),"$")}),row.vars = 4:3, col.vars = 2:1)
    fftab1 <- toLatex(ftabu, booktabs = TRUE, fontsize = "footnotesize", escapeFUN = MyescapeLatex,  
                      caption = paste0("One-sided empirical error rate (\\%) for lower (first to third columns) and upper (fourth to sixth columns) confidence limits, ", 
                                       switch(met, pot = "peaks-over-threshold", bm = "block maximum"), " method ", str,"  
The distributions (from top to bottom) are Burr ($F_1$), Weibull ($F_2$), generalized gamma ($F_3$), normal ($F_4$), lognormal ($F_5$), Student $t$ ($F_6$), ", 
                                       switch(met, 
                                              bm = "$\\mathsf{GEV}(\\xi =0.1)$ ($F_7$), Gumbel ($F_8$) and $\\mathsf{GEV}(\\xi=-0.1)$ ($F_9$)",
                                              pot = "$\\mathsf{GP}(\\xi =0.1)$ ($F_7$), exponential ($F_8$) and $\\mathsf{GP}(\\xi=-0.1)$ ($F_9$)"),"."),
                      label = paste0("tab_nomerr_l",levl,"_",met), align = "*{2}{l}*{12}{r}", center = TRUE, placement  = "t!")
    cat(fftab1, file = paste0("Table_tem_", met,"_",
                              switch(met, 
                                     pot = switch(levl, 60,40,20), 
                                     bm = paste0("m",switch(levl, 30, 45, 90)))
                              , "_nomerr.tex"), append = FALSE)
    
  }
}


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")

# Produce figures of empirical error rate

#for(quant in c("quant","Nquant","Nmean")){
for(met in c("pot", "bm")){
    if(met == "pot"){
      tabu <- nomerror[,c(1,3),met,,,c(1:6,10:12)]
  dimnames(tabu) <- list("Error rate" = c("0.5","2.5","5", "5","2.5","0.5"),
                         Parameter = c("Quantile", "$N$-obs. mean"),
                          "Sample size" = c(60,40,20),
                          method = c("Wald", "profile", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                          "$F$" = paste0("$F_",c(1:9),"$") 
                         #"Burr", "Weibull", "Gen. gamma", "Normal", "Log-normal", "Student ", "GP(0.1)", "Exponential","GP(-0.1)")
        )
    } else if(met == "bm"){
      tabu <- nomerror[,c(1,3),met,,,c(1:9)]
      dimnames(tabu) <- list("Error rate" = c("0.5","2.5","5", "5","2.5","0.5"),
                             Parameter = c("Quantile", "$N$-obs. mean"),
                             "Sample size" = c(60,40,20),
                             method = c("Wald", "profile", "\\textsc{tem}","Severini (\\textsc{tem})", "Severini (cov.)"),
                             "$F$" = paste0("$F_",c(1:9),"$"))
   
    }
  # Figures 

  # transform matrix to data frame (long)
  tabudf <- reshape2::melt(tabu, id.vars="row.names")
  tabudf$error <- (tabudf$value - tabudf$`Error rate`)/tabudf$`Error rate`
  tabudf$`Relative error rate` <- rep(c(0.5,2.5,5,95,97.5,99.5), length.out = nrow(tabudf))
  tabudf$`Relative error rate` <- as.factor(tabudf$`Relative error rate`)
  tabudf$tail <- factor(rep(c(rep("lower",3),rep("upper",3)), length.out = nrow(tabudf)))
  
  tabudf$Parameter <- as.factor(tabudf$Parameter)
  levels(tabudf$Parameter) <- c("return level","mean")
  tabudf$method <- as.factor(tabudf$method)
  levels(tabudf$method) <- c("Wald","profile","tem","modif (tem)","modif (cov)")
  tabudf$dist <- as.factor(tabudf$`$F$`) 
  levels(tabudf$dist) <- switch(met,
                                pot = c("Burr", "Weibull", "ggamma", "normal", "log-normal", "Student ", "GP(0.1)", "exp","GP(-0.1)"),
                                bm = c("Burr", "Weibull", "ggamma", "normal", "log-normal", "Student ", "GEV(0.1)", "GEV(0)","GEV(-0.1)"))
  tabudf$`$F$` <- NULL
  tabudf$`Sample size` <- as.factor(tabudf$`Sample size`)
  tabudf$scenario <- interaction(tabudf$`Sample size`,tabudf$method)
  dist.labs <- c("$\\xi=0.1$", "$\\xi=0$", "$\\xi=-0.1$")
  levels(tabudf$dist)[7:9] <- dist.labs
  if(figures){
    setwd("figures")
    fig <- paste0("error_rate_true",met,".tex")
    tikz(fig, width = 8, height = 5, standAlone = TRUE)
  }
  p1 <- ggplot(data = tabudf[(tabudf$Parameter == "mean")&(tabudf$`Sample size` %in% c("20","60"))&(tabudf$dist %in% levels(tabudf$dist)[7:9])&(tabudf$method != "Wald"),],
                 mapping = aes(x = `Relative error rate`, y = error, 
                               color = method, group=scenario)) +
    geom_abline(slope = 0, intercept = 0) + 
    geom_point() +
    stat_summary(fun=identity, geom="line") +
    facet_grid(dist ~ `Sample size`*tail, 
               scales = "free_x")   +
    theme_minimal() +
    ylab("relative error rate of one-sided confidence limit")+ 
    xlab("quantile") + 
    scale_color_viridis(discrete = TRUE, option = "D")+
    theme(legend.position="bottom")
  print(p1)
  if(figures){
  dev.off()
  system(paste0("lualatex '", getwd(),"/", fig, "'; rm *.aux; rm *.log"))
  setwd("..")
  }
}
