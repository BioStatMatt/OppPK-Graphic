## Bayes prediction
library('mvtnorm')  ## dmvnorm
library('cubature') ## adaptIntegrate
library('ellipse')  ## ellipse
source('model.R')

## Here we assume that a one-compartment model is sufficient, and that 
## prior distribution for the log of the v_1 and k_10 parameters have 
## a normal distribution with mean and covariance as follows. This is a
## weakly informative prior that has arisen from a prior study (in submission).
lpk_mean_d <- c(lv_1=3.223, lk_10=-1.650)
lpk_vcov_d <- 300 * matrix(c( 0.00167, -0.00128,
                             -0.00128,  0.00154), 2,2)

## Prior distribution for error variance 
err_shape_d <- 2
err_scale_d <- sqrt(754)/2
err_mean_d  <- err_shape_d * err_scale_d

lpr_mean_d  <- c(lpk_mean_d, log(err_mean_d)) 

## Given a series of concentration measurements for a particular patient, 
## we can use the corresponding posterior distribution to summarize pharmacodynamic
## target attainment, and make predictions about how the patient will respond to
## dose/frequency changes. The function below evaluates the log posterior density
## as a function of the two PK parameters (lk_10, lv_1).

## lpr - log parameter vector: c(lk_10=??, lv_1=??)
## err - error standard deviation
## mu  - prior mean c(lk_10=??, lv_1=??)
## sig - prior variange-covariance (2x2 PD matrix)
log_prior <- function(lpr, mu=lpk_mean_d, sig=lpk_vcov_d,
                      err_shape=err_shape_d, err_scale=err_scale_d)
  dmvnorm(lpr[1:2], mu, sig, log = TRUE) + 
  dgamma(lpr[3], shape=err_shape, scale=err_scale, log=TRUE)

## lpr - log parameter vector: c(lk_10=??, lv_1=??)
## ivt - list describing sequence of doses
## dat - concentration data: data.frame(time_h=??, conc_mg_dl=??)
## err - error standard deviation [mg/dl]
## ini - initial concentrations
log_likelihood <- function(lpr, ivt, dat, ini=c(0,0)) {
  epr <- exp(lpr)
  sol <- pk_solution(v_1=epr[1], k_10=epr[2], ivt=ivt)
  dat$pconc_g_l   <- sol(dat$time_h)[1,]
  dat$pconc_mg_dl <- 1000*dat$pconc_g_l
  with(dat, sum(dnorm(conc_mg_dl, pconc_mg_dl, epr[3], log=TRUE)))
}

## lpr - log parameter vector: c(lk_10=??, lv_1=??, lerr=??)
## ivt  - list describing sequence of doses
## dat  - concentration data: data.frame(time_h=??, conc_mg_dl=??)
log_posterior <- function(lpr, ivt, dat) {
  dat <- na.omit(dat)
  if(nrow(dat) < 1) {
    log_prior(lpr) 
  } else {
    log_prior(lpr) + log_likelihood(lpr, ivt, dat)
  }
}


## plot the contours of a univariate function of two variables ('x' and 'y')
## fun - function with prototype 'function(x, ...) {...}' that outputs a single 
##       numeric value, and where x is a numeric vector of length 2, i.e., 'c(x,y)'
## xlm - limits of 'x' variable
## ylm - limits of 'y' variable
## pts - number of grid points to evaluate in the 'x' and 'y' directions
## nrm - normalize the values returned from 'fun'?
## ... - arguments passed to 'contour'
contourf <- function(fun, xlm, ylm, pts=c(100,100), nrm=TRUE, ...) {
  xsq <- seq(xlm[1], xlm[2], length.out=pts[1])
  ysq <- seq(ylm[1], ylm[2], length.out=pts[2])
  grd <- expand.grid(x=xsq, y=ysq)
  val <- apply(grd, 1, fun)
  if(nrm) { ## normalize
    ## approximate normalization constant
    ## using cubature::adaptIntegrate
    nrm <- adaptIntegrate(fun, 
      lowerLimit=c(xlm[1],ylm[1]),
      upperLimit=c(xlm[2],ylm[2]))$integral
    val <- val/nrm
  }
  val <- matrix(val, nrow=pts[1], ncol=pts[2])
  contour(x=xsq,y=ysq,z=val, ...)
}

## suppose that a patient received the default dosing pattern
## and has the following concentration measurements at time 1h
#dat <- data.frame(time_h = c(1,4,8), conc_mg_dl = c(82.7,50.4,30.6))
dat <- data.frame(time_h = c(1,4), conc_mg_dl = c(82.7,80.4))
#dat <- data.frame(time_h = c(1), conc_mg_dl = c(82.7))
#dat <- data.frame(time_h = c(8), conc_mg_dl = c(30.6))

## then the Bayes estiamte for v_1 and k_10 is the posterior mode, which
## we can find as follows:
est <- optim(lpr_mean_d, log_posterior, ivt=ivt_d,
             dat=dat, control = list(fnscale=-1), hessian=TRUE)

## plot prior predictions for the concentration-time curve
par(mfrow=c(1,2))
tms <- seq(1e-3, 8, 0.1)
sol_prior <- pk_solution()
con_prior <- sol_prior(tms)*1000
plot(tms, con_prior[1,], xlab="Time (h)",
     ylab="Central Concentration (mg/dL)",
     ylim=range(con_prior), type='l', lty=2,
     main="Concentration vs. Time")
sol_posterior <- pk_solution(k_10=exp(est$par['lk_10']),
                             v_1=exp(est$par['lv_1']))
con_posterior <- sol_posterior(tms)*1000
lines(tms, con_posterior[1,])
points(dat$time_h, dat$conc_mg_dl, pch=16)
legend('topright', c('Prior', 'Posterior', 'Measured'),
       lty=c(2,1,NA), pch=c(NA, NA, 16), bty='n')


## plot contours of the prior distribution for k_10 and v_1
ell_prior <- ellipse(log_par_vcov_d, centre=lpk_mean_d)
## ell_posterior is computed using posterior Laplace approximation
ell_posterior <- ellipse(solve(-est$hessian), centre=est$par)
plot(range(c(ell_prior[,1],ell_posterior[,1])),
     range(c(ell_prior[,2],ell_posterior[,2])), type="n", 
     xlab=expression(log~V[1]),
     ylab=expression(log~k[10]),
     main="95% Credible Region")
lines(ell_prior, lty=2)
lines(ell_posterior, lty=1)
legend('topright', c('Prior', 'Posterior'),
       lty=c(2,1), bty='n')

## plot contours of the prior/posterior distribution for k_10 and v_1
# contourf(function(x) exp(log_prior(x)),
#          xlm, ylm, pts=c(60,60), lty=2,
#          levels=c(0.05, 0.20, 0.80),
#          xlab=expression(log~V[1]),
#          ylab=expression(log~k[10]),
#          main="Prior Contours")


# contourf(function(x) exp(log_posterior(x, ivt, dat)),
#          xlm, ylm, pts, lty=1,
#          levels=c(0.05, 0.20, 0.80),
#          xlab=expression(log~V[1]),
#          ylab=expression(log~k[10]),
#          main="Posterior Contours")


## compute finite-difference gradient (c.f., nlme::fdHess)
fdGrad <- function (pars, fun, ...,
                    .relStep = (.Machine$double.eps)^(1/2), 
                    minAbsPar = 0) {
  ##pars <- as.numeric(pars)
  npar <- length(pars)
  incr <- ifelse(abs(pars) <= minAbsPar, .relStep, 
                 (abs(pars)-minAbsPar) * .relStep)
  ival <- do.call(fun, list(pars, ...))
  diff <- rep(0,npar)
  sapply(1:npar, function(i) {
    del <- rep(0,npar)
    del[i] <- incr[i]
    (do.call(fun, list(pars+del, ...))-ival)/incr[i]
  })
}

grd <- fdGrad(est$par, function(pars) {
  tms <- seq(1e-3, 8, 0.1)
  sol <- pk_solution(v_1=exp(pars[1]), k_10=exp(pars[2]))
  log(sol(tms)[1,]*1000) ## g/l -> ug/ml
})

sde <- sqrt(diag(grd %*% solve(-est$hessian) %*% t(grd)))
sde <- ifelse(is.nan(sde), 0, sde)

tms <- seq(1e-3, 8, 0.1)
sol_posterior <- pk_solution(k_10=exp(est$par['lk_10']),
                             v_1=exp(est$par['lv_1']))
con_posterior <- apply(sol_posterior(tms)*1000, 2, function(x) pmax(0,x))
par(mfrow=c(1,1))
plot(tms, con_posterior[1,], xlab="Time (h)",
     ylab="Central Concentration (mg/dL)",
     ylim=c(0, max(exp(log(con_posterior[1,])+1.96*sde), na.rm=TRUE)),
     type='l', lwd=2, col="#882255",
     main="Concentration vs. Time")
points(dat$time_h, dat$conc_mg_dl, pch=16)
polygon(c(tms,rev(tms)), 
        c(exp(log(con_posterior[1,]) + 1.96*sde),
          rev(exp(log(con_posterior[1,]) - 1.96*sde))),
        col="#CC667755", border=NA)
legend('topright', c("Predicted", "95% Credible Band", "Measured"),
       lwd=c(2,4,NA), pch=c(NA,NA,16), col=c('#882255','#CC667755','black'),
       border=NA, bty='n')

