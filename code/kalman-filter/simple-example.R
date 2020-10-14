#!/usr/bin/env R

library(tidyverse)
source("covidhub-common.R")

forecast_dates <- "2020-10-12"
hopdir <- file.path("hopkins", forecast_dates)
tdat <- load_hopkins(hopdir) 

nyc <- tdat %>% filter(location == 36) %>% filter(target_type == "wk ahead inc case")




nyc2 <- nyc %>% mutate(time = lubridate::decimal_date(target_end_date))
case_data <- nyc2 %>% select(time, value) %>% rename(reports = value)

cov_data0 <- tibble(time = c(case_data$time, max(case_data$time) * 1.05))

bspline_basis <- pomp::periodic.bspline.basis(
  cov_data0$time,
  nbasis = 6,
  degree = 3,
  period = 1,
  names = "xi%d"
) %>%
  as_tibble()

cov_data <- bind_cols(cov_data0, bspline_basis)

PsystemSEIR <- function(pvec, covf,
                        init.vars,
                        time.steps = c(1995, 1995 + 1 / 52)){
  
  PModel <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
      xi1t <- covf$xi1(t)
      xi2t <- covf$xi2(t)
      xi3t <- covf$xi3(t)
      xi4t <- covf$xi4(t)
      xi5t <- covf$xi5(t)
      xi6t <- covf$xi6(t)
      beta_t <- beta_mu * (1 + exp(xi1t * b1 + xi2t * b2 + xi3t * b3 + xi4t * b4 + xi5t * b5 + xi6t * b6))
      eta <- 365 / 8  
      gamma <- 365 / 5
      S <- exp(lS)
      E <- exp(lE)
      I <- exp(lI)
      F <- rbind(c(-beta_t * I / N,    0, -beta_t * S / N, 0),
                 c( beta_t * I / N, -eta,  beta_t * S / N, 0),
                 c(                0,  eta,            -gamma, 0),
                 c(                0,    0,             gamma, 0))
      
      
      f <- c(0, beta_t * (S / N) * (I / N), eta * E / N, gamma * I / N)
      Q <- rbind(c(f[1] + f[2],       -f[2],           0,     0),
                 c(      -f[2], f[2] + f[3],       -f[3],     0),
                 c(          0,       -f[3], f[3] + f[4], -f[4]),
                 c(          0,           0,       -f[4],  f[4]))
      
      P <- rbind(c(Pss, Pse, Psi, Psc),
                 c(Pse, Pee, Pei, Pec),
                 c(Psi, Pei, Pii, Pic),
                 c(Psc, Pec, Pic, Pcc))
      
      dlS <- (- beta_t * S * I / N) / S
      dlE <- (beta_t * S * I / N - eta * E) / E
      dlI <- (iota + eta * E -  gamma * I) / I
      dC <- (gamma * I)
      dP <-  F %*% P + P %*% t(F) + Q
      
      list(c(dlS = dlS, dlE = dlE, dlI = dlI, dC = dC, dPss = dP[1,1], dPse = dP[1,2], 
             dPsi = dP[1,3], dPsc = dP[1,4], dPee = dP[2,2], dPei = dP[2,3], 
             dPec = dP[2,4], dPii = dP[3,3], dPic = dP[3,4], dPcc = dP[4,4]))
    })
  }
  deSolve::lsoda(init.vars, time.steps, PModel, pvec)
}

genfun <- function(y) {
  approxfun(cov_data$time, y)
}
covf <- apply(cov_data, 2, genfun)

pvec <-   params <- c(
  N = 20e6,
  beta_mu = 500,
  beta_sd = 0.001,
  b1 = 3,
  b2 = 3,
  b3 = 6,
  b4 = 5,
  b5 = 2,
  b6 = 0,
  iota = 2,
  rho = 0.5,
  S_0 = 20e6, 
  E_0 = 10,
  I_0 = 10,
  tau = 0.001
)


init.vars <- c(lS=log(3e-2 * 20e6), lE=log(1.6e-4 * 20e6), lI=log(1.6e-4 * 20e6), C = 0,
               Pss = 1, Pse = 0, Psi = 0, Psc = 0, Pee = 1, Pei = 0, Pec = 0, Pii = 1, Pic = 0, Pcc = 1)

out <- PsystemSEIR(pvec = pvec, covf = covf, init.vars = init.vars, time.steps = case_data$time)

iterate_f_and_P <- function(xhat, PN, pvec, covf, time.steps){
  P <- PN / pvec["N"]
  xhat_trans <- c(log(xhat[c("S", "E", "I")]), xhat["C"])
  if(!all(is.finite(xhat_trans))) {
    browser()
  }
  names(xhat_trans)[1:3] <- c("lS", "lE", "lI")
  init.vars <- c(xhat_trans, Pss = P[1,1], Pse = P[1,2], 
                 Psi = P[1,3], Psc = P[1,4], Pee = P[2,2], Pei = P[2,3], Pec = P[2,4], 
                 Pii = P[3,3], Pic = P[3,4], Pcc = P[4,4])
  ret <- PsystemSEIR(pvec = pvec, init.vars = init.vars, covf, time.steps)[2, ]
  xhat_new <- c(exp(ret[c("lS", "lE", "lI")]), ret["C"])
  names(xhat_new)[1:3] <- c("S", "E", "I")
  P_new  <- with(as.list(ret),        
                 rbind(c(Pss, Pse, Psi, Psc),
                       c(Pse, Pee, Pei, Pec),
                       c(Psi, Pei, Pii, Pic),
                       c(Psc, Pec, Pic, Pcc)))
  PN_new <- P_new * pvec["N"]
  list(xhat = xhat_new, PN = PN_new)
}

xhat <-  c(S=(3e-2 * 20e6), E=(1.6e-4 * 20e6), I=(1.6e-4 * 20e6), C = 0)
P <- with(as.list(init.vars[-c(1:4)]),        
          rbind(c(Pss, Pse, Psi, Psc),
                c(Pse, Pee, Pei, Pec),
                c(Psi, Pei, Pii, Pic),
                c(Psc, Pec, Pic, Pcc)))

iterate_f_and_P(xhat = xhat, PN = P, pvec = pvec, covf = covf, 
                time.steps = c(2020.1, 2020.1 + 1 / 52))

# Initialize
xhat0 <- matrix(xhat, ncol = 1)
rownames(xhat0) <- names(xhat)
Phat0 <- P
z_1 <-  case_data$reports[1]
H <- matrix(c(0, 0, 0, pvec["rho"]), ncol = 4)
R <- z_1 * (1 - pvec["rho"])

# Predict
XP_1_0 <- iterate_f_and_P(xhat0[, 1], PN = Phat0, pvec = pvec, covf = covf,
                          time.steps = case_data$time[c(1, 2)])
xhat_1_0 <- XP_1_0$xhat
P_1_0 <- XP_1_0$PN
# Update

K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R[1])
ytilde_1 <- z_1 - H %*% xhat_1_0
xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
P_1_1 <- (diag(4) - K_1 %*% H) %*% P_1_0

## Now calculate for each step in simulation


T <- nrow(case_data)
z <- case_data$reports

ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(4, T))
rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat)
P_kk <- P_kkmo <- array(NA_real_, dim = c(4, 4, T))

K[, 1] <- K_1
xhat_kkmo[, 1] <- xhat_1_0
xhat_kk[, 1] <- xhat_1_1
P_kk[, , 1] <- P_1_1
P_kkmo[, , 1] <- P_1_0
S[, 1] <- H %*% P_kkmo[, , 1] %*% t(H) + R
ytilde_kk[, 1] <- z[1] - H %*% xhat_kk[, 1]
ytilde_k[, 1] <- ytilde_1

for (i in seq(2, T)){
  xhat_init <- xhat_kk[, i - 1]
  xhat_init["C"] <- 0
  PNinit <- P_kk[,,i - 1]
  PNinit[, 4] <- PNinit[4, ] <- 0
  XP <- iterate_f_and_P(xhat_init, PN = PNinit, pvec = pvec, covf = covf,
                        time.steps = case_data$time[c(i - 1, i)])
  xhat_kkmo[, i] <- XP$xhat
  P_kkmo[, , i] <- XP$PN
  #R <- xhat_kkmo["C", i] * pvec["rho"] * (1 - pvec["rho"])
  R <- max(5, z[i - 1] * (1 - pvec["rho"]))
  S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
  K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
  ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
  xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
  xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
  P_kk[, , i] <- (diag(4) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
  ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
}

log_lik <- function(Sigma, resids){
  -0.5 * sum(resids ^ 2 / Sigma + log(Sigma) + log(2 * pi))
}

log_lik(S, ytilde_k)





