#!/usr/bin/env R

options(repos = "https://cran.microsoft.com/snapshot/2020-02-29/")
library(pomp)

source("../model-setup/makepompmodel.R")
source("../model-setup/setparsvars.R")

pomp_data <- data.frame(time = c(1, 100), cases = NA, hosps = NA, deaths = NA)

n_knots <- 2
est_these_pars = c("log_sigma_dw", "min_frac_dead", "max_frac_dead", "log_half_dead",
                   "log_theta_cases", "log_theta_deaths")
est_these_inivals = c("E1_0", "Ia1_0", "Isu1_0", "Isd1_0")
# est_these_inivals = ""  # to not estimate any initial values
knot_coefs <-  paste0("b", 1:n_knots)
est_these_pars <- c(est_these_pars, knot_coefs)

# Set the parameter values and initial conditions
par_var_list <- setparsvars(est_these_pars = est_these_pars, 
                            est_these_inivals = est_these_inivals,
                            population = 9e6,
                            rnaught = 4)  # set R0 at beginning of epidemic

pomp_covar <- covariate_table(
  t = pomp_data$time,
  seas = bspline.basis(
    x=t,
    nbasis=n_knots,
    degree=0
  ),
  rel_beta_change = 2,
  trend_sim = 1,  # this is a placeholder only needed for simulation
  fit = 0,  # 1 = fitting; 0 = simulating
  times="t",
  order = "constant"
)

pomp_model <- makepompmodel(par_var_list = par_var_list, 
                            pomp_data = pomp_data, 
                            pomp_covar = pomp_covar,
                            n_knots = n_knots)

coef(pomp_model) <- par_var_list$allparvals
out <- pomp::simulate(pomp_model, times = 1:100)

sd <- as.data.frame(out)

reactions <- list(infectS=list(paste0("rate = S * exp(log_beta_s) * ((Isd1 + Isd2 + Isd3 + Isd4)",
                                     " + (Isu1 + Isu2 + Isu3 + Isu4)",
                                     " + (Isd1 + Isd2 + Isd3 + Isd4)",
                                     " + 1/(1 + exp(trans_e)) * (E1 + E2 + E3 + E4)",
                                     " + 1/(1 + exp(trans_a)) * (Ia1 + Ia2 + Ia3 + Ia4)",
                                     " + 1/(1 + exp(trans_c)) * (C1 + C2 + C3 + C4)",
                                     " + 1/(1 + exp(trans_h)) * (H1 + H2 + H3 + H4));"), c(S=-1, E1=1)),
                 exitE1=list("rate = E1 * exp(log_g_e);", c(E1=-1, E2 = 1)),
                 exitE2=list("rate = E2 * exp(log_g_e);", c(E2=-1, E3 = 1)),
                 exitE3=list("rate = E3 * exp(log_g_e);", c(E3=-1, E4 = 1)),
                 enterIa1=list("rate = E4 * exp(log_g_e) * 1/(1+exp(frac_asym));", c(E4=-1, Ia1 = 1)),
                 enterIsu1=list("rate = E4 * exp(log_g_e) * (1 - 1/(1+exp(frac_asym))) * (1 - detect_frac);", c(E4=-1, Isu1 = 1)),
                 enterIsd1=list("rate = E4 * exp(log_g_e) * (1 - 1/(1+exp(frac_asym))) * detect_frac;", c(E4=-1, Isd1 = 1)),
                 exitIa1=list("rate = Ia1 * exp(log_g_a);", c(Ia1=-1, Ia2=1)),
                 exitIa2=list("rate = Ia2 * exp(log_g_a);", c(Ia2=-1, Ia3=1)),
                 exitIa3=list("rate = Ia3 * exp(log_g_a);", c(Ia3=-1, Ia4=1)),
                 exitIa4=list("rate = Ia4 * exp(log_g_a);", c(Ia4=-1, R=1)),
                 exitIsu1=list("rate = Isu1 * exp(log_g_su);", c(Isu1=-1, Isu2=1)),
                 exitIsu2=list("rate = Isu2 * exp(log_g_su);", c(Isu2=-1, Isu3=1)),
                 exitIsu3=list("rate = Isu3 * exp(log_g_su);", c(Isu3=-1, Isu4=1)),
                 exitIsu4=list("rate = Isu4 * exp(log_g_su);", c(Isu4=-1, R=1)),
                 exitIsd1=list("rate = Isd1 * g_sd;", c(Isd1=-1, Isd2=1)),
                 exitIsd2=list("rate = Isd2 * g_sd;", c(Isd2=-1, Isd3=1)),
                 exitIsd3=list("rate = Isd3 * g_sd;", c(Isd3=-1, Isd4=1)),
                 exitIsd4=list("rate = Isd4 * g_sd;", c(Isd4=-1, C1=1, C_new = 1)),
                 exitC1=list("rate = C1 * g_c;", c(C1=-1, C2=1)),
                 exitC2=list("rate = C2 * g_c;", c(C2=-1, C3=1)),
                 exitC3=list("rate = C3 * g_c;", c(C3=-1, C4=1)),
                 exitC4enterR=list("rate = C4 * g_c * (1 - 1/(1+exp(frac_hosp)));", c(C4=-1, R=1)),
                 exitC4enterH=list("rate = C4 * g_c * 1/(1+exp(frac_hosp));", c(C4=-1, H1=1, H_new = 1)),
                 exitH1=list("rate = H1 * exp(log_g_h);", c(H1 = -1, H2 = 1)),
                 exitH2=list("rate = H2 * exp(log_g_h);", c(H2 = -1, H3 = 1)),
                 exitH3=list("rate = H3 * exp(log_g_h);", c(H3 = -1, H4 = 1)),
                 exitH4enterR=list("rate = H4 * exp(log_g_h) * (1 - frac_dead);", c(H4 = -1, R=1)),
                 exitH4enterD=list("rate = H4 * exp(log_g_h) * frac_dead;", c(H4 = -1, D=1, D_new = 1)))

sir_rproc <- do.call(pomp::gillespie_hl, reactions)

stoich <- sir_rproc@v
mode(stoich) <- "character"

pull_rate <- function(x){
  ret <- strsplit(x[[1]], split = "=")[[1]][2]
  ret <- strsplit(ret, split = ";")[[1]][1]
  paste0("(", ret, ")")
}
rates <- lapply(reactions, pull_rate)[colnames(stoich)]

macrorate <- function(x) {
  test <- x != 0
  paste(rates[test], x[test], sep = "*", collapse = "+")
}
rhs <- apply(stoich, 1, macrorate)
lhs <- paste0("D", names(rhs))

frates <- sapply(rates, function(x) parse(text = x))



R <- length(rates)
N <- nrow(stoich)
statenms <- rownames(stoich)
jvars <- outer(statenms, statenms, function(x,y) paste0("xi", x, "_", y))
A <- matrix("", N, N)

for (i in seq(1, N)){
  for (k in seq(1, N)){
    for (j in seq(1, R)){
      dfj_dthetak <- D(frates[j], statenms[k])
      dval <- deparse(dfj_dthetak, width.cutoff = 500L)
      sval <- stoich[i, j]
      if (dval != "0" && sval != "0"){
        A[i, k] <- paste(A[i, k], "+",  sval, "*", dval)
      }
    }
    if (nchar(A[i, k]) > 0){
      A[i, k] <- paste0(" + (", A[i, k], ")")
    }
  }
}
Aformulas <- sapply(as.character(A), function(x) parse(text = x))
vf_formulas <- sapply(rhs, function(x) parse(text = x))


eval_model <- function(xhat, params, time, N, vf_expressions, Aformulas, stoich, frates){
  aparams <- as.list(c(xhat, t = time, params))
  
  diag_speedup <- with(aparams, 1 + exp(log_max_diag)  *  (t ^ exp(log_diag_inc_rate)) / (exp(log_half_diag) ^ exp(log_diag_inc_rate))  + 
                         (t ^ exp(log_diag_inc_rate)))
  g_sd <- with(aparams, diag_speedup * exp(log_g_sd))  ## shortened time in symptomatic stage prior to diagnosis
  g_c  <- with(aparams, exp(log_g_c) / diag_speedup) ## increased time in symptomatic stage post diagnosis
  detect_frac <- with(aparams, 1 / (1+exp(max_detect_par)) * (t ^ exp(log_detect_inc_rate))  / ( (exp(log_half_detect) ^ exp(log_detect_inc_rate)) + (t ^ exp(log_detect_inc_rate))) + 
                        exp(base_detect_frac))
  frac_dead <- with(aparams,  1/(1+exp(min_frac_dead)) + (1/(1+exp(max_frac_dead)) - 1/(1+exp(min_frac_dead)) ) * (1 - t/(t+exp(log_half_dead)) ))
  
  fvals <- sapply(frates, function(x) with(aparams, eval(x)))
  B <- stoich %*% diag(fvals) %*% t(stoich)
  
  aparams2 <- c(list(g_sd = g_sd, g_c = g_c, detect_frac = detect_frac), aparams)
  tmpf <- function(x){
    ret <- with(aparams2, eval(x))
    if(is.null(ret)){
      ret <- 0
    }
    ret
  }
  Avals <- sapply(Aformulas, function(x) tmpf(x))
  vf <- sapply(vf_expressions, function(x) tmpf(x))
  list(vectorfield = vf,
       Jacobian = matrix(Avals, nrow = N, ncol = N),
       B = B)
}

xhat0 <- numeric(length(statenms))
xhat0 <- as.numeric(sd[1, statenms])
names(xhat0) <- statenms

stoichn <- stoich
mode(stoichn) <- "numeric"
mvals <- eval_model(xhat = xhat0, params = par_var_list$allparvals, N = N, vf_expressions = vf_formulas, Aformulas = Aformulas, time = 2,
                    frates = frates, stoich = stoichn)

iterate_f_and_P <- function(xhat, P, params, N, vf, jac, time,  stoich, frates, dt){
  inds <- seq_len(N)
  PModel <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
      xhat <- x[inds]
      xhat[xhat < 0 ] <- 0
      P <- matrix(x[-inds], nrow = N, ncol = N)
      zinds <- which(diag(P) < 0)
      P[zinds, zinds] <- 0
      mvals <- eval_model(xhat = xhat, params = parms, N = N, vf_expressions = vf, Aformulas = jac, time = t,
                          frates = frates, stoich = stoich)
      dxhat <-  mvals$vectorfield
      dP <-  (mvals$Jacobian %*% P + P %*% t(mvals$Jacobian) + mvals$B)
      list(c(dxhat, as.numeric(dP)))
    })
  }
  init.vars <- c(xhat, as.numeric(P))
  out <- deSolve::lsoda(init.vars, c(0, 0 + dt), PModel, params)[2, -1]
  xhat_next <- out[inds]
  P_next <- matrix(out[-inds], nrow = N, ncol = N)
  list(xhat = xhat_next, P = P_next)
}

P0 <- diag(N)
P0[20,20] <- 0
xP <- iterate_f_and_P(xhat = xhat0, P = P0, params = par_var_list$allparvals, N = N, vf = vf_formulas, 
                         jac = Aformulas, time = 2, frates = frates, stoich = stoichn, dt = 1)

kfnll <-
  function(data,
           params, 
           dt = 1,
           N, 
           vf,
           jac,
           frates,
           stoich,
           xhat0,
           P0,
           obsvars = 20,
           R = 1, 
           just_nll = FALSE) {
    
    H <- matrix(0, nrow = 1, ncol = N)
    H[1, obsvars] <- 1
    
    # Predict
    xP <- iterate_f_and_P(xhat = xhat0, P = P0, params = params, N = N, vf = vf, 
                          jac = jac, frates = frates, stoich = stoich, 
                          time = data$time[1] - dt, dt = dt)
    xhat_1_0 <- xP$xhat
    P_1_0 <- xP$P
    # Update
    
    K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R)
    ytilde_1 <- data$cases[1] - H %*% xhat_1_0
    xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
    P_1_1 <- (diag(N) - K_1 %*% H) %*% P_1_0
    
    T <- nrow(data)
    nobsv <- length(obsvars)
    ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(nobsv, T))
    K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(N, T))
    rownames(xhat_kk) <- names(xhat0)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(N, N, T))
    
    K[, 1] <- K_1
    xhat_kkmo[, 1] <- xhat_1_0
    xhat_kk[, 1] <- xhat_1_1
    P_kk[, , 1] <- P_1_1
    P_kkmo[, , 1] <- P_1_0
    S[, 1] <- H %*% P_kkmo[, , 1] %*% t(H) + R
    ytilde_kk[, 1] <- data$cases[1] - H %*% xhat_kk[, 1]
    ytilde_k[, 1] <- ytilde_1
    
    for (i in seq(2, T)){
      xhat <- xhat_kk[, i - 1]
      xhat[obsvars] <- 0
      P <- P_kk[,, i - 1]
      P[obsvars, ] <- 1e-4
      P[, obsvars] <- 1e-4
      xP <- iterate_f_and_P(xhat, P, params = params, N = N, vf = vf, 
                            jac = jac, frates = frates, stoich = stoich, 
                            time = data$time[i - 1], dt = data$time[i]  - data$time[i - 1])
      xhat_kkmo[, i] <- xP$xhat
      P_kkmo[, , i] <- xP$P
      S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
      K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
      ytilde_k[, i] <- data$cases[i] - H %*% xhat_kkmo[, i, drop = FALSE]
      xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
      P_kk[, , i] <- (diag(N) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
      ytilde_kk[i] <- data$cases[i] - H %*% xhat_kk[, i, drop = FALSE]
    }
    
    nll <- 0.5 * sum(ytilde_k ^ 2 / S + log(S) + log(2 * pi))
    if (!just_nll){
      list(nll = nll, xhat_kk = xhat_kk, P_kk = P_kk, ytilde_k = ytilde_k)
    } else {
      nll
    }
    
  }

kfo <- kfnll(data = sd[-1, ], xhat = xhat0, P0 = P0, params = par_var_list$allparvals, N = N, vf = vf_formulas, 
                      jac = Aformulas, frates = frates, stoich = stoichn, dt = 1)

plot(sd$time[-1], kfo$xhat_kk["E1",], xlab = "Time", ylab = "First latent compartment", type = 'l')
points(sd$time[-1], sd$E1[-1], col = "blue")
legend("topright", pch = c(1, NA), lty = c(NA, 1), col = c("blue", 1), legend = c("Simulation (Truth)", "Kalman filter based on cases"))

plot(sd$time[-1], kfo$xhat_kk["C_new",], xlab = "Time", ylab = "cases", type = 'l')
points(sd$time[-1], sd$cases[-1], col = "blue")
legend("topright", pch = c(1, NA), lty = c(NA, 1), col = c("blue", 1), legend = c("Simulation (Truth)", "Kalman filter based on cases"))


nll_wrap <- function(logbeta){
  params <- par_var_list$allparvals
  params["log_beta_s"] <- logbeta
  res <- kfnll(data = sd[-1, ], xhat = xhat0, P0 = P0, params = params, N = N, vf = vf_formulas, 
        jac = Aformulas, frates = frates, stoich = stoichn, dt = 1, just_nll = TRUE)
  print(paste("logbeta = ", logbeta, ",  nll = ", res, ", time = ", Sys.time()))
  res
}

library(bbmle)

m0 <- mle2(minuslogl = nll_wrap, start = list(logbeta = -12), method = "L-BFGS-B", lower = -20, upper = -10)

