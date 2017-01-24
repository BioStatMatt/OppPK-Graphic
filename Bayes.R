## Bayes prediction
library('mvtnorm')  ## dmvnorm
library('cubature') ## adaptIntegrate
library('ellipse')  ## ellipse (only needed for example)
source('model.R')

## Here we assume that a one-compartment model is sufficient, and that 
## prior distribution for the log of the v_1 and k_10 parameters have 
## a normal distribution with mean and covariance as follows. This is a
## weakly informative prior that has arisen from a prior study (in submission).
lpk_mean_d <- c(lv_1=3.223, lk_10=-1.650)
lpk_vcov_d <- 300 * matrix(c( 0.00167, -0.00128,
                              -0.00128,  0.00154), 2,2)

## Prior distribution for error standard deviation (normal)
ler_mean_d <- 2.33
ler_sdev_d <- 0.32

## Vectof of all three parameters
lpr_mean_d  <- c(lpk_mean_d, ler_mean_d) 

## Given a series of concentration measurements for a particular patient, 
## we can use the corresponding posterior distribution to summarize pharmacodynamic
## target attainment, and make predictions about how the patient will respond to
## dose/frequency changes. The function below evaluates the log posterior density
## as a function of the two PK parameters (lk_10, lv_1).

## Log prior density
## lpr - log parameter vector
## err - error standard deviation
## mu  - prior mean for PK parameters
## sig - prior variange-covariance for PK parameters (2x2 PD matrix)
## ler_mean - prior mean for error variance
## ler_sdev - prior standard deviation for error variance
log_prior <- function(lpr, mu=lpk_mean_d, sig=lpk_vcov_d,
                      ler_mean=ler_mean_d, ler_sdev=ler_sdev_d)
  dmvnorm(lpr[1:2], mu, sig, log = TRUE) + 
  dnorm(lpr[3], mean=ler_mean, sd=ler_sdev, log=TRUE)

## Log likelihood function
## lpr - log parameter vector
## ivt - list describing sequence of doses
## dat - concentration data: data.frame(time_h, conc_mg_dl)
## ini - initial concentrations
log_likelihood <- function(lpr, ivt, dat, ini=c(0,0)) {
  epr <- exp(lpr)
  sol <- pk_solution(v_1=epr[1], k_10=epr[2], ivt=ivt)
  dat$pconc_g_l   <- sol(dat$time_h)[1,]
  dat$pconc_mg_dl <- 1000*dat$pconc_g_l
  with(dat, sum(dnorm(conc_mg_dl, pconc_mg_dl, epr[3], log=TRUE)))
}

## Log posterior density
## lpr - log parameter vector
## ivt  - list describing sequence of doses
## dat  - concentration data: data.frame(time_h, conc_mg_dl)
log_posterior <- function(lpr, ivt, dat) {
  dat <- na.omit(dat)
  if(nrow(dat) < 1) {
    log_prior(lpr) 
  } else {
    log_prior(lpr) + log_likelihood(lpr, ivt, dat)
  }
}

## Function to compute finite-difference gradient (c.f., nlme::fdHess)
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

# Function to calculate the fraction of time spent above 'th'
mic_stat <- function(pk_pars, ivt, tms, con, th){
  #Values of the plotted posterior concentrations
  conc <- con[1,]
  
  # Get PK solution equation evaluated at parameters
  soln <- pk_solution(v_1=exp(pk_pars[1]), k_10=exp(pk_pars[2]), ivt=ivt)
  # Use PK solution to define function that computes 
  #   the concentrations centered by threshold
  f_mic <- function(times = tms){ 
    val <- apply(soln(times)*1000, 2, function(x) pmax(0,x)) - th
    return(val[1,]) 
  }
  
  # Convert start/end dose times to numeric vectors
  ibe <- sapply(ivt, `[[`, 'begin')
  ied <- sapply(ivt, `[[`, 'end')
  
  # Initialize time spent above 4*mic
  t_above <- 0
  
  # Time between start of dose and end of dose
  for(i in 1:length(ibe)){
    # Concentration values at interval endpoints, centered by threshold
    cb <- ifelse(i > 1, conc[tms == ibe[i]] - th, 0-64)
    ce <- unique(conc[tms == ied[i]] - th) #Printing two copies for some reason? Inserted unique to fix for now
    if(cb < 0 && ce > 0){
      # Crosses to above threshold during interval
      root <- uniroot(f_mic, lower = ibe[i], upper = ied[i])$root #Must move from below to above
      t_above <- t_above + (ied[i] - root)
    }else if(cb >= 0 && ce >= 0){ 
      # Above during whole interval
      t_above <- t_above + (ied[i] - ibe[i])
    }  
  }
  
  # Time between end of one dose and start of the next
  for(j in 1:length(ied)){
    ce <- unique(conc[tms == ied[j]] - th)
    c_next <- ifelse(j < length(ied), conc[tms == ibe[j+1]] - th, conc[tms == max(tms)] - th)
    if(ce > 0 && c_next < 0){
      # Crosses to below threshold during interval
      ulim <- ifelse(j < length(ied), ibe[j+1], max(tms))
      root <- uniroot(f_mic, lower = ied[j], upper = ulim)$root #Must move from above to below
      t_above <- t_above + (root - ied[j])
    }else if(ce >= 0 && c_next >= 0){ 
      # Above during whole interval
      t_add <- ifelse(j < length(ied), ibe[j+1] - ied[j], max(tms) - ied[j])
      t_above <- t_above + t_add
    }  
  }
  
  # Use time spent above 4*mic to compute proportion
  frac_time <- t_above/max(tms)
  
  return(frac_time)
}


## Plot posterior estimated concentration-time curve with approximate
## Wald-type posterior (1-alp)*100% credible bands
## est - object returned from 'optim' with 'hessian=TRUE'
## ivt - list describing sequence of doses
## dat - concentration data: data.frame(time_h, conc_mg_dl)
## alp - credible level (1-alp)%
## cod - additional time following last dose (h)
plot_post_conc <- function(est, ivt, dat, alp=0.05, cod=12, thres=64) {

  ## Compute plotting times
  ## - ensure peak and trough times
  ## - avoid time zero
  tms <- sapply(ivt, function(x) c(x$begin, x$end))
  tms <- c(tms, max(tms)+cod)
  tms <- unlist(sapply(1:(length(tms)-1), function(i) { 
    s1 <- seq(tms[i], tms[i+1], 1/10)
    if(tms[i+1] %% 1/10)
      s1 <- c(s1, tms[i+1])
    return(s1)
    }))
  tms <- pmax(1e-3, tms)
  
  ## Approximate standard deviation of log concentration-time curve
  grd <- fdGrad(est$par, function(pars) {
    sol <- pk_solution(v_1=exp(pars[1]), k_10=exp(pars[2]), ivt=ivt) 
    log(sol(tms)[1,]*1000) ## mulitply by 1000: g/l -> ug/ml
  })
  sde <- sqrt(diag(grd %*% solve(-est$hessian) %*% t(grd)))
  sde <- ifelse(is.nan(sde), 0, sde)
  
  ## Plot posterior estiamte
  col_bg <- "#AAAAB5"
  col_fg <- "#3E465A"
  sol <- pk_solution(v_1=exp(est$par[1]), k_10=exp(est$par[2]), ivt=ivt) 
  con <- apply(sol(tms)*1000, 2, function(x) pmax(0,x))
  par(mfrow=c(1,1))
  plot(tms, con[1,], xlab="Time (h)", ylab="Central Concentration (mg/dL)",
       ylim=c(0, max(exp(log(con[1,])+qnorm(1-alp/2)*sde), na.rm=TRUE)),
       type='n', main="Concentration vs. Time")
  
  ## Plot 95% credible bands
  polygon(c(tms,rev(tms)), 
          c(exp(log(con[1,]) + qnorm(1-alp/2)*sde),
            rev(exp(log(con[1,]) - qnorm(1-alp/2)*sde))),
          col=col_bg, border=NA)
  
  ## Plot posterior estimate
  lines(tms, con[1,], lwd=2, col=col_fg)
  
  ## Plot measured points
  points(dat$time_h, dat$conc_mg_dl, pch=16)
  
  ## Create legend
  legend('topleft', c("Predicted", "95% Credible Band", "Measured"),
         lwd=c(2,4,NA), pch=c(NA,NA,16), col=c(col_fg,col_bg,'black'),
         border=NA, bty='n')
  
  
  ## Add MIC statistic information to plot
  # Obtain statistic
  frac_mic <- mic_stat(pk_pars = est$par, ivt, tms, con, th = thres)
  
  # SE of logit(statistic)
  grd_mic <- fdGrad(est$par, function(pars) {
    mic <- mic_stat(pk_pars = pars, ivt, tms, con, th = thres) 
    log(mic/(1-mic)) ## constrain between 0 and 1
  })
  sde_mic <- sqrt(diag(t(grd_mic) %*% solve(-est$hessian) %*% grd_mic))
  
  # Get CI for logit transformed statistic then backtransform to original scale
  ci_logit_mic <- log(frac_mic/(1-frac_mic)) + c(-1,1)*qnorm(1-alp/2)*sde_mic
  ci_mic <- exp(ci_logit_mic)/(1 + exp(ci_logit_mic))
  
  #Plotting elements for mic statistic
  abline(h = thres, lty = 2)
  
  legend("topright", bty = 'n',
         legend = c(paste("FracTime>4*MIC:", round(frac_mic, 3)), 
                    paste("95% CI: (", round(ci_mic[1], 3), ",", round(ci_mic[2], 3), ")")))
}


## Example
## suppose that a patient received the default dosing pattern
## and has the following concentration measurements at time 1h
#dat <- data.frame(time_h = c(1,4,40), conc_mg_dl = c(82.7,80.4,60))
## dat <- data.frame(time_h = c(1,4,8), conc_mg_dl = c(82.7,50.4,30.6))
## dat <- data.frame(time_h = c(1), conc_mg_dl = c(82.7))
## dat <- data.frame(time_h = c(8), conc_mg_dl = c(30.6))
#system.time({
#est <- optim(lpr_mean_d, log_posterior, ivt=ivt_d,
#             dat=dat, control = list(fnscale=-1), hessian=TRUE)
#plot_post_conc(est, ivt_d, dat)
#})

## plot prior and posterior predictions for the concentration-time curve
# par(mfrow=c(1,2))
# tms <- seq(1e-3, 8, 0.1)
# sol_prior <- pk_solution()
# con_prior <- sol_prior(tms)*1000
# plot(tms, con_prior[1,], xlab="Time (h)",
#      ylab="Central Concentration (mg/dL)",
#      ylim=range(con_prior), type='l', lty=2,
#      main="Concentration vs. Time")
# sol_posterior <- pk_solution(k_10=exp(est$par['lk_10']),
#                              v_1=exp(est$par['lv_1']))
# con_posterior <- sol_posterior(tms)*1000
# lines(tms, con_posterior[1,])
# points(dat$time_h, dat$conc_mg_dl, pch=16)
# legend('topright', c('Prior', 'Posterior', 'Measured'),
#        lty=c(2,1,NA), pch=c(NA, NA, 16), bty='n')

## plot approximate 95% credible region for prior and posterior for k_10 and v_1
# ell_prior <- ellipse(lpk_vcov_d, centre=lpk_mean_d)
# ell_posterior <- ellipse(solve(-est$hessian), centre=est$par) ## Laplace approximation
# plot(range(c(ell_prior[,1],ell_posterior[,1])),
#      range(c(ell_prior[,2],ell_posterior[,2])), type="n",
#      xlab=expression(log~V[1]),
#      ylab=expression(log~k[10]),
#      main="95% Credible Region")
# lines(ell_prior, lty=2)
# lines(ell_posterior, lty=1)
# legend('topright', c('Prior', 'Posterior'),
#        lty=c(2,1), bty='n')
