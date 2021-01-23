###########
##TEM derivations, N-exceedances mean parametrization
dat, z, xi, N, u, k, sigma = var('dat','z','xi','N','u','k','sigma')
gpdll = -log(sigma)-(1+1/xi)*log(1+xi/sigma*dat)
gpdNmll = gpdll.substitute(sigma=z*xi/(gamma(N+1)*gamma(1-xi)/gamma(N+1-xi)-1))
gpdF = 1-(1+xi/sigma*dat)^(-1/xi)
gpdNmF = gpdF.substitute(sigma=z*xi/(gamma(N+1)*gamma(1-xi)/gamma(N+1-xi)-1))
gpdNmF2 = 1-(1+dat/z*(gamma(N+1)*gamma(1-xi)/gamma(N+1-xi)-1))^(-1/xi)
gpdNmf = gpdNmF.diff(dat)
Vxi = -gpdNmF.diff(xi)/gpdNmf

#Tangent derivative with respect to nuisance
gpdNmll.diff(dat, xi)
gpdNmll.diff(dat, z)
gpdNmll.diff(dat)


## TEM derivations for the median of the

## TEM derivations for the maximum of k Pareto variates
gpdFm = gpdF^k
gpdF.substitute(



## TEM derivations for probability of exceedance for GEV using r-largest observations
dat, datr, pex, xi, mu, mut, sigma, r, t, z, t0 = var('dat','datr','pex','xi','mu','mut','sigma','r','t','z', 't0')
#for the case r=2, we can write the log-likelihood in full
rlargll = -r*log(sigma) - (1+1/xi)*log(1+xi/sigma*(dat-mu-mut*t)) - (1+1/xi)*log(1+xi/sigma*(datr-mu-mut*t)) - (1+xi/sigma*(datr-mu-mut*t))^(-1/xi)
intensF = (1+xi/sigma*(dat-mu-mut*t))^(-1/xi)
intensrF = (1+xi/sigma*(datr-mu-mut*t))^(-1/xi)
# z is the determined level for which we seek probability of exceedance pex
iF = intensF.substitute(mu = -mut*t0+z-sigma/xi*((-log(1-pex))^(-xi)-1))
irF = intensrF.substitute(mu = -mut*t0+z-sigma/xi*((-log(1-pex))^(-xi)-1))

dydmut = -iF.diff(mut)/iF.diff(dat)
dydpex = -iF.diff(pex)/iF.diff(dat)
dydsigma = -iF.diff(sigma)/iF.diff(dat)
dydxi = -iF.diff(xi)/iF.diff(dat)

dyrdmut = -(irF.diff(mut)-iF.diff(mut))/irF.diff(datr)
dyrdpex = -(irF.diff(pex)-iF.diff(pex))/irF.diff(datr)
dyrdsigma = -(irF.diff(sigma)-iF.diff(sigma))/irF.diff(datr)
dyrdxi = -(irF.diff(xi)-iF.diff(xi))/irF.diff(datr)


rlargpexll = rlargll.substitute(mu = -mut*t0+z-sigma/xi*((-log(1-pex))^(-xi)-1))
dlldy = rlargpexll.diff(dat);
dlldyr = rlargpexll.diff(datr);

#Mixed derivatives
dlldydmut = dlldy.diff(mut);
dlldydpex = dlldy.diff(pex);
dlldydsigma = dlldy.diff(sigma);
dlldydxi = dlldy.diff(xi);
dlldyrdmut = dlldyr.diff(mut);
dlldyrdpex = dlldyr.diff(pex);
dlldyrdsigma = dlldyr.diff(sigma);
dlldyrdxi = dlldyr.diff(xi);

# Score vector

rlargpexlldmut = rlargpexll.diff(mut)
rlargpexlldpex = rlargpexll.diff(pex)
rlargpexlldsigma = rlargpexll.diff(sigma)
rlargpexlldxi = rlargpexll.diff(xi)

# Observed information matrix

rlargpexlldmutdmut = rlargpexll.diff(mut,2)
rlargpexlldmutdpex = rlargpexll.diff(mut,pex)
rlargpexlldmutdsigma = rlargpexll.diff(mut,sigma)
rlargpexlldmutdxi = rlargpexll.diff(mut,xi)
rlargpexlldpexdpex = rlargpexll.diff(pex, 2)
rlargpexlldpexdsigma = rlargpexll.diff(pex, sigma)
rlargpexlldpexdxi = rlargpexll.diff(pex, xi)
rlargpexlldsigmadsigma = rlargpexll.diff(sigma, 2)
rlargpexlldsigmadxi = rlargpexll.diff(sigma, xi)
rlargpexlldxidxi = rlargpexll.diff(xi,2)

# TEM derivations for scale of GEV
mu,sigma,xi,dat,y,z = var('mu','sigma','xi','dat','y','z')
gevll = -log(sigma)-(1/xi+1)*log(1+xi*(x-mu)/sigma)-(1+xi*(x-mu)/sigma)^(-1/xi)
#CDF and PDF of the GEV distribution
F = exp(-(1+xi*(y-mu)/sigma)^(-1/xi))
f = exp(-(1+xi*(y-mu)/sigma)^(-1/xi))/sigma*(1+xi*(y-mu)/sigma)^(-1/xi-1)
#V function
(-F.diff(mu)/f).simplify_full()
(-F.diff(sigma)/f).simplify_full()
(-F.diff(xi)/f).simplify_full()
#Phi: data space derivative
d0 = gevll.diff(x)
d01 = gevll.diff(x,mu)
d02 = gevll.diff(x,sigma)
d03 = gevll.diff(x,xi)


#### January 9th, 2021 Derivation of TEM for Maiquetia

y, s, z, Nu, ny, u, mu, sigma, xi, qp, N, ql = var('y','s','z','Nu','ny','u','mu','sigma','xi','qp','N','ql')
Lambd(x) = ny*(1+xi*(x-mu)/sigma)^(-1/xi)
lambd(x) = ny*(1+xi*(x-mu)/sigma)^(-1/xi-1)/sigma
intens = ny*(1+xi*(y-mu)/sigma)^(-1/xi-1)/sigma
# For NHPP, we can get V through derivative of -dot(Lambda) wrt parameters
intensQ = intens.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
Vy_qp = intensQ.diff(qp)
Vy_sigma = intensQ.diff(sigma)
Vy_xi = intensQ.diff(xi)
# Limits for case xi=0 (technically not needed, because only evaluated at MLE)
Vy_qp0 = lim(Vy_qp, xi=0)
Vy_sigma0 = lim(Vy_sigma, xi=0)
Vy_xi0 = lim(Vy_xi, xi=0)

# HOA components (psi and dpsi) for the NHPP psi = log{-dot(Lambda)}/{-dot(Lambda)}
psif = log(intensQ)/intensQ
# Derivatives of psi wrt parameters
psi_qp = psif.diff(qp)
psi_sigma = psif.diff(sigma).simplify_full()
psi_xi = psif.diff(xi)

# HOA components for the (truncated) GP
dens_fc_nu =  (lambd(y)/Lambd(s)).substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
cdf_fc_nu =  ((Lambd(s) - Lambd(y))/Lambd(s)).substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
(cdf_fc_nu.diff(y) - dens_fc_nu).simplify_full()
V_fc_nu_qp = (- cdf_fc_nu.diff(qp)/dens_fc_nu).simplify_full()
V_fc_nu_sigma = (- cdf_fc_nu.diff(sigma)/dens_fc_nu).simplify_full()
V_fc_nu_xi = (- cdf_fc_nu.diff(xi)/dens_fc_nu).simplify_full()

dens_fc =  (lambd(y)/(Lambd(u) - Lambd(s))).substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
cdf_fc =  ((Lambd(u) - Lambd(y))/(Lambd(u) - Lambd(s))).substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
(cdf_fc.diff(y) - dens_fc).simplify_full()
V_fc_qp = (- cdf_fc.diff(qp)/dens_fc).simplify()
V_fc_sigma = (- cdf_fc.diff(sigma)/dens_fc).simplify()
V_fc_xi = (- cdf_fc.diff(xi)/dens_fc).simplify()

# The sample space derivative kills the conditioning terms, so contribution for y_nu is the same

log(dens_fc).diff(y).simplify_full()
log(dens_fc).diff(y, qp).simplify_full()
log(dens_fc).diff(y, sigma).simplify_full()
log(dens_fc).diff(y, xi).simplify_full()

# GEV contribution (the mev package has a function @gevN.Vfun, etc., but it profiles over the scale parameter
gevll = -log(sigma)-(1/xi+1)*log(1+xi*(z-mu)/sigma)-(1+xi*(z-mu)/sigma)^(-1/xi)
#CDF and PDF of the GEV distribution
F = 1-exp(-(1+xi*(z-mu)/sigma)^(-1/xi))
f = exp(-(1+xi*(z-mu)/sigma)^(-1/xi))/sigma*(1+xi*(z-mu)/sigma)^(-1/xi-1)
Fqp = F.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
fqp = f.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
#V function
(-Fqp.diff(qp)/fqp).simplify_full()
(-Fqp.diff(sigma)/fqp).simplify_full()
(-Fqp.diff(xi)/fqp).simplify()
#Phi: data space derivative
phi = gevll.substitute(mu = qp +sigma/xi*(1-(-N/log(ql))^xi)).diff(z)
dphi_qp = gevll.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi)).diff(z,qp)
dphi_sigma = (gevll.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi)).diff(z,sigma)).simplify_full()
dphi_xi = gevll.substitute(mu =  qp+sigma/xi*(1-(-N/log(ql))^xi)).diff(z,xi)


