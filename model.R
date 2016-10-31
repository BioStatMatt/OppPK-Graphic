## This code implements a two-compartment model with PK parameters
## k_10 - [1/h] elimination rate from central compartment
## k_12 - [1/h] distribution rate from central to peripheral compartment
## k_21 - [1/h] distribution rate from peripheral to central compartment
## v_1  - [l]   volume of central compartment
## v_2  - [l]   volume of peripheral compartment
## The additional parameter k_R [g/h] represents the rate of intravenous
## infusion; typically a fixed, known quantity.

## The model is expressed as a system of ODEs as follows, where m1 and m2
## are the masses of drug [mg] in the central and peripheral compartments,
## respectively. This system is based on the concept of conservation of
## mass.
## dm1/dt = -k_10*m1 - k_12*m1 + k_21*m2 + k_R
## dm2/dt =          + k_12*m2 - k_21*m2

## Substituting m1=c1*v_1 and m2=c2*v_2 allows the model to be expressed
## in terms of drug concentration, rather than mass:
## dc1/dt = -k_10*c1 - k_12*c1         + k_21*c2*v_2/v_1 + k_R/v_1
## dc2/dt =          + k_12*c1*v_1/v_2 - k_21*c2

## However, we generally only measure the central compartment
## concentration. Thus, there is usually no information about
## the v_2 parameter, and we can simplify the model by assuming that 
## v_1=v_2, which gives the following, with four free parameters:
## dc1/dt = -k_10*c1 - k_12*c1 + k_21*c2 + k_R/v_1
## dc2/dt =          + k_12*c1 - k_21*c2

## default PK parameters
## setting k_12_d and k_21_d close to zero (e.g., 0.001) mimics
## a one compartment model
k_10_d <- 0.35
k_12_d <- 3.50
k_21_d <- 1.50
v_1_d  <- 10.0

## A sequence of intravenous infusions is represented by a list of 'begin'
## and 'end' times [h] for the infusion, as well as the rate of infusion [g/h]
## The following list encodes a sequence of five half-hour infusions totaling
## 3g of drug each (6 g/h), separated by 8 hours:
ivt_d <- list(list(begin=0.0, end=0.5, k_R=6),
              list(begin=8.0, end=8.5, k_R=6),
              list(begin=16.0, end=16.5, k_R=6),
              list(begin=24.0, end=24.5, k_R=6),
              list(begin=32.0, end=32.5, k_R=6))

## Specify the initial drug concentrations in each compartment
init_d <- c(0,0)

## The 'pk_basic_solution' function implements the 
## two-compartment model solution for fixed parameters k_10,
## k_12, k_21, v_1, k_R, and c_0 (initial concentrations). 
## This function is based on the solution to a nonhomogeneous
## linear system of ODEs. The basic solution could also be 
## computed using numerical integration, which is commonly used
## for more complicated models:
pk_basic_solution <-
function(k_10, k_12, k_21, v_1, k_R, c_0=c(0,0)) {
  K <- matrix(c(-(k_10+k_12), k_21, k_12, -k_21),
              2,2, byrow=TRUE)
  T_K <- -(k_10 + k_12 + k_21)
  D_K <- (k_10+k_12)*k_21 - k_12*k_21
  lambda <- T_K/2 + c(-1,1)*(T_K^2/4 - D_K)^(1/2)
  N_21 <- (k_21+lambda[2])/k_21
  N_22 <- (k_21+lambda[1])/k_21
  N <- matrix(c(-1,N_21,-1,N_22),2,2)
  gamma <- c(k_R/k_10/v_1, k_12*k_R/k_10/k_21/v_1) 
  r <- c(( N_22*(c_0[1]-(k_R/v_1/k_10)) +
             (c_0[2]-(k_R/v_1/k_10)*k_12/k_21)),
         (-N_21*(c_0[1]-(k_R/v_1/k_10)) -
            (c_0[2]-(k_R/v_1/k_10)*k_12/k_21)))
  r <- r * k_21/(diff(lambda))
  c_1 <- function(t)
    gamma[1] + r[1]*N[1,1]*exp(lambda[1]*t) +
               r[2]*N[1,2]*exp(lambda[2]*t)
  c_2 <- function(t)
    gamma[2] + r[1]*N[2,1]*exp(lambda[1]*t) +
               r[2]*N[2,2]*exp(lambda[2]*t)
  return(list(c_1=c_1, c_2=c_2))
}

## The 'pk_solution' function impelments a piecewise solution
## to the two-compartment model, taking into account the sequence
## of dosing events:
pk_solution <-
  function(k_10=k_10_d, k_12=k_12_d, k_21=k_21_d, v_1=v_1_d,
           ivt=ivt_d, init=init_d) {
  ## create a list of event times
  ibe <- sapply(ivt, `[[`, 'begin')
  ied <- sapply(ivt, `[[`, 'end')
  prd <- sort(unique(c(0, c(ibe,ied), Inf)))
  rits <- list()
  ## compute basic solution in each interval
  for(i in 1:(length(prd)-1)) {
    civt <- sapply(ivt, function(iv) {
      if(prd[i] >= iv$begin && prd[i] < iv$end) {
        iv$k_R
      } else { 0 }
    })

    rit <- list(begin=prd[i], end=prd[i+1],
                idose=sum(civt))
    if(i == 1) {
      rit$init <- init
    } else {
      rit$init <- c(rits[[i-1]]$c_1(rits[[i-1]]$end-rits[[i-1]]$begin),
                    rits[[i-1]]$c_2(rits[[i-1]]$end-rits[[i-1]]$begin))
    }
    
    sol <- pk_basic_solution(k_10, k_12, k_21, v_1,
                             k_R=rit$idose, c_0=rit$init)
    
    rits[[i]] <- c(rit, sol)
  }
  
  # both c_1 and c_2
  ret <- function(tms) {
    sapply(tms, function(t) {
      val <- NA
      for(rit in rits) {
        if(t >= rit$begin && t <= rit$end) {
          val <- c(rit$c_1(t-rit$begin),
                   rit$c_2(t-rit$begin))
          break
        }
      }
      val
    })
  }
  return(ret)
}

## The following code chunk plots the model solution for the default
## PK parameters, dosing schedule, and initial concentrations.
# sol <- pk_solution()
# tms <- seq(0, 40, 0.01)
# con <- sol(tms)
# plot(tms, con[1,], xlab="Time (h)", ylab="Concentration (g/L)",
#      ylim=range(con), type='l')
# lines(tms, con[2,], type='l', lty=2)
# legend('topleft', c('Central', 'Peripheral'), lty=c(1,2), bty='n')

## The code below creates a similar figure using model parameters that
## are more typical for hospitalized patients who have received 
## intravenous piperacillin. Here we assume that a one-compartment model
## is sufficient, and that the log of the v_1 and k_10 parameters have 
## a normal distribution in the population, with mean and covariance as
## follows:
log_par_mean_d <- c(lv_1=3.223, lk_10=-1.650)
log_par_vcov_d <- 754 * matrix(c( 0.00167, -0.00128,
                                 -0.001285, 0.00154), 2,2)

## set the corresponding default parameters
k_10_d <- exp(log_par_mean_d['lk_10'])
v_1_d  <- exp(log_par_mean_d['lv_1'])
k_12_d <- 0.001
k_21_d <- 0.001


## sample (N=500) the population distribution of lv_1 and lk_10
# log_par_samp <- t(log_par_mean_d + t(chol(log_par_vcov_d)) %*% matrix(rnorm(1000),2,500))
# sol_samp <- apply(log_par_samp, 1, function(x) pk_solution(v_1=exp(x[1]), k_10=exp(x[2])))
# tms <- seq(0, 40, 0.01)
# con <- pk_solution()(tms)
# plot(tms, con[1,], xlab="Time (h)", ylab="Central Concentration (g/L)",
#      ylim=range(con), type='l')
# for(sol_ in sol_samp)
#   lines(tms, sol_(tms)[1,], col=rgb(0,0,0,0.))