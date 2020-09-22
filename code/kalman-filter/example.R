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
                            population = 1e6,
                            rnaught = 6)  # set R0 at beginning of epidemic

pomp_covar <- covariate_table(
  t = pomp_data$time,
  seas = bspline.basis(
    x=t,
    nbasis=n_knots,
    degree=0
  ),
  rel_beta_change = 1,
  trend_sim = 10,  # this is a placeholder only needed for simulation
  fit = 1,  # 1 = fitting; 0 = simulating
  times="t",
  order = "constant"
)

pomp_model <- makepompmodel(par_var_list = par_var_list, 
                            pomp_data = pomp_data, 
                            pomp_covar = pomp_covar,
                            n_knots = n_knots)

coef(pomp_model) <- par_var_list$allparvals
out <- pomp::simulate(pomp_model, times = 1:30)

sd <- as.data.frame(out)

reactions <- list(infectS=list(paste0("rate = psi * phi * 1 / Nthis * beta * ((Isd1 + Isd2 + Isd3 + Isd4)",
                                     " + (Isu1 + Isu2 + Isu3 + Isu4)",
                                     " + (Isd1 + Isd2 + Isd3 + Isd4)",
                                     " + bL * (E1 + E2 + E3 + E4)",
                                     " + bIa * (Ia1 + Ia2 + Ia3 + Ia4)",
                                     " + bC * (C1 + C2 + C3 + C4)",
                                     " + bH * (H1 + H2 + H3 + H4));"), c(S=-1, E1=1)),
                 exitE1=list("rate = E1 * gammaE;", c(E1=-1, E2 = 1)),
                 exitE2=list("rate = E2 * gammaE;", c(E2=-1, E3 = 1)),
                 exitE3=list("rate = E3 * gammaE;", c(E3=-1, E4 = 1)),
                 enterIa1=list("rate = E4 * gammaE * a;", c(E4=-1, Ia1 = 1)),
                 enterIsu1=list("rate = E4 * gammaE * (1 - a) * (1 - q);", c(E4=-1, Isu1 = 1)),
                 enterIsd1=list("rate = E4 * gammaE * (1 - a) * q;", c(E4=-1, Isd1 = 1)),
                 exitIa1=list("rate = Ia1 * gammaIa;", c(Ia1=-1, Ia2=1)),
                 exitIa2=list("rate = Ia2 * gammaIa;", c(Ia2=-1, Ia3=1)),
                 exitIa3=list("rate = Ia3 * gammaIa;", c(Ia3=-1, Ia4=1)),
                 exitIa4=list("rate = Ia4 * gammaIa;", c(Ia4=-1, R=1)),
                 exitIsu1=list("rate = Isu1 * gammaIsu;", c(Isu1=-1, Isu2=1)),
                 exitIsu2=list("rate = Isu2 * gammaIsu;", c(Isu2=-1, Isu3=1)),
                 exitIsu3=list("rate = Isu3 * gammaIsu;", c(Isu3=-1, Isu4=1)),
                 exitIsu4=list("rate = Isu4 * gammaIsu;", c(Isu4=-1, R=1)),
                 exitIsd1=list("rate = Isd1 * gammaIsd * s;", c(Isd1=-1, Isd2=1)),
                 exitIsd2=list("rate = Isd2 * gammaIsd * s;", c(Isd2=-1, Isd3=1)),
                 exitIsd3=list("rate = Isd3 * gammaIsd * s;", c(Isd3=-1, Isd4=1)),
                 exitIsd4=list("rate = Isd4 * gammaIsd * s;", c(Isd4=-1, C1=1, C_new = 1)),
                 exitC1=list("rate = C1 * gammaC / s;", c(C1=-1, C2=1)),
                 exitC2=list("rate = C2 * gammaC / s;", c(C2=-1, C3=1)),
                 exitC3=list("rate = C3 * gammaC / s;", c(C3=-1, C4=1)),
                 exitC4enterR=list("rate = C4 * gammaC / s * (1 - h);", c(C4=-1, R=1)),
                 exitC4enterH=list("rate = C4 * gammaC / s * h;", c(C4=-1, H1=1, H_new = 1)),
                 exitH1=list("rate = H1 * gammaH;", c(H1 = -1, H2 = 1)),
                 exitH2=list("rate = H2 * gammaH;", c(H2 = -1, H3 = 1)),
                 exitH3=list("rate = H3 * gammaH;", c(H3 = -1, H4 = 1)),
                 exitH4enterR=list("rate = H4 * gammaH * (1 - m);", c(H4 = -1, R=1)),
                 exitH4enterD=list("rate = H4 * gammaH * m;", c(H4 = -1, D=1, D_new = 1)))

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

eval_vf_and_jac <- function(xhat, params, N, frates, Aformulas){
  aparams <- as.list(c(xhat, params))
  tmpf <- function(x){
    ret <- with(aparams, eval(x))
    if(is.null(ret)){
      ret <- 0
    }
    ret
  }
  Avals <- sapply(Aformulas, function(x) tmpf(x))
  vf <- sapply(frates, function(x) tmpf(x))
  list(vectorfield = vf,
       Jacobian = matrix(Avals, nrow = N, ncol = N))
}

