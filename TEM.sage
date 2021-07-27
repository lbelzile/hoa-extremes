
################################################################################################
##TEM derivations, N-exceedances median parametrization
## Generalized Pareto distribution
dat, xi, N, u, k, tau, zetau, qp, q = var('dat','xi','N','u','k','tau','zetau','qp','q')
gpdll = -log(tau)-(1+1/xi)*log(1+xi/tau*dat)
gpdNqll = gpdll.substitute(tau=(qp-u)*xi/((1-q^(1/(zetau*N)))^(-xi)-1))
gpdF = 1-(1+xi/tau*dat)^(-1/xi)
gpdNqF = gpdF.substitute(tau=(qp-u)*xi/((1-q^(1/(zetau*N)))^(-xi)-1))
gpdNqf = gpdNqF.diff(dat)
Vqp = (-gpdNqF.diff(qp)/gpdNqf).simplify()
Vxi = (-gpdNqF.diff(xi)/gpdNqf).simplify()

#Tangent derivative with respect to nuisance
gpdNqll.diff(dat, xi)
gpdNqll.diff(dat, qp)
gpdNqll.diff(dat)

################################################################################################

#### Derivation of TEM for Maiquetia 

y, s, z, Nu, ny, u, mu, sigma, xi, qp, N, ql = var('y','s','z','Nu','ny','u','mu','sigma','xi','qp','N','ql')
intens = ny*(1+xi*(y-mu)/sigma)^(-1/xi-1)/sigma
# Non-homogeneous Poisson process (NHPP)
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
psi_qp = (Vy_qp/intensQ^2).simplify()
psi_sigma = (Vy_sigma/intensQ^2).simplify_full()
psi_xi = Vy_xi/intensQ^2


# GEV contribution (the mev package has a function @gevN.Vfun, etc., 
# but it's parametrized in terms of (mu, qp, xi) rather than (qp, sigma, xi)
gevll = -log(sigma)-(1/xi+1)*log(1+xi*(z-mu)/sigma)-(1+xi*(z-mu)/sigma)^(-1/xi)
#CDF and PDF of the GEV distribution
F = 1-exp(-(1+xi*(z-mu)/sigma)^(-1/xi))
f = exp(-(1+xi*(z-mu)/sigma)^(-1/xi))/sigma*(1+xi*(z-mu)/sigma)^(-1/xi-1)
Fqp = F.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
fqp = f.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
#V function
Vgev_qp = (-Fqp.diff(qp)/fqp).simplify_full()
Vgev_sigma = (-Fqp.diff(sigma)/fqp).simplify_full()
Vgev_xi = (-Fqp.diff(xi)/fqp).simplify()
#Phi: data space derivative
phi = gevll.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi)).diff(z)
dphi_qp = gevll.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi)).diff(z, qp)
dphi_sigma = (gevll.substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi)).diff(z, sigma)).simplify()
dphi_xi = gevll.substitute(mu =  qp+sigma/xi*(1-(-N/log(ql))^xi)).diff(z, xi)


# Left and right truncated exceedances 
# the conditional distributions are truncated generalized Pareto
Lambd(x) = (1+xi*(x-mu)/sigma)^(-1/xi)
lambd(x) = (1+xi*(x-mu)/sigma)^(-1/xi-1)/sigma
dens_fc_lt =  (lambd(y)/Lambd(s)).substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
cdf_fc_lt =  ((Lambd(s) - Lambd(y))/Lambd(s)).substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
dens_fc_rt =  (lambd(y)/(Lambd(u) - Lambd(s))).substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))
cdf_fc_rt =  ((Lambd(u) - Lambd(y))/(Lambd(u) - Lambd(s))).substitute(mu = qp+sigma/xi*(1-(-N/log(ql))^xi))

# Sufficient directions for left truncation

Vyl_qp = (-cdf_fc_lt.diff(qp)/dens_fc_lt).simplify_full()
Vyl_sigma = (-cdf_fc_lt.diff(sigma)/dens_fc_lt).simplify_full()
Vyl_xi = (-cdf_fc_lt.diff(xi)/dens_fc_lt).simplify()

# Sufficient directions for right truncation

Vyr_qp = (-cdf_fc_rt.diff(qp)/dens_fc_rt).simplify_full()
Vyr_sigma = (-cdf_fc_rt.diff(sigma)/dens_fc_rt).simplify_full()
Vyr_xi = (-cdf_fc_rt.diff(xi)/dens_fc_rt).simplify()

psi = log(dens_fc_lt).diff(y).simplify_full()
psi_qp =  psi.diff(qp)
psi_sigma = psi.diff(sigma)
psi_xi = psi.diff(xi)

################################################################################################


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


################################################################################################

##HOA derivations, endpoint of generalized Pareto
dat, zeta, xi, u, sigma, x, l = var('dat','zeta','xi','u','sigma','x', 'l')
gpdll = -log(sigma)-(1+1/xi)*log(1+xi/sigma*dat)
gpdF = 1-(1+xi/sigma*dat)^(-1/xi)
surv(x) = (1-x/zeta)^(-1/xi)
gpdendptll = gpdll.substitute(sigma = -xi*zeta)
gpdendptF = gpdF.substitute(sigma = -xi*zeta)

indll1 = (gpdendptll - log(surv(l))).diff(xi)
indll2 = (log(surv(dat)) - log(surv(l))).diff(xi)


