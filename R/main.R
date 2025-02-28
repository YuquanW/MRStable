#' Title
#'
#' @param beta_exp Exposure
#' @param beta_out Outcome
#' @param se_exp Std error exposure
#' @param se_out Std error outcome
#' @param pi_thr Threshold
#' @param rd Randomized lasso
#' @param over.dispersion Over dispersion
#'
#' @return A list
#' @export
#'
#' @examples
#' attach(bmi.bmi)
#' MRStable(beta.exposure, beta.outcome, se.exposure, se.outcome)
MRStable <- function(beta_exp,
                     beta_out,
                     se_exp,
                     se_out,
                     pi_thr = 0.6,
                     rd = 50,
                     over.dispersion = T) {
  m <- length(beta_exp)
  iv.sig <- .stablility_selection_stg1(beta_exp,
                                       se_exp,
                                       rd,
                                       pi_thr)
  if (length(iv.sig)==0) {
    stop("No significant IVs were selected.")
  } else if (length(iv.sig) <= 10) {
    warning("The number of significant IVs is less than 10.")
  }
  t <- min(abs(beta_exp)[iv.sig]/se_exp[iv.sig])
  #t <- quantile(abs(beta_exp)/se_exp, 1-sqrt(0.5/m))
  beta_exp_cor <- .fix_point(beta_exp, se_exp, iv.sig, t)
  list(iv.sig, .divw(beta_exp_cor, beta_out, se_exp, se_out, over.dispersion))
}

#' @rdname MRStable
#' @export

MRStable_V <- function(beta_exp,
                       beta_out,
                       se_exp,
                       se_out,
                       pi_thr = 0.6,
                       rd = 50,
                       majority = F,
                       over.dispersion.stg1 = T,
                       over.dispersion.stg2 = T) {
  m <- length(beta_exp)
  iv.sig <- .stablility_selection_stg1(beta_exp,
                                        se_exp,
                                        rd,
                                        pi_thr)

  t <- min(abs(beta_exp)[iv.sig]/se_exp[iv.sig])
  # t <- quantile(abs(beta_exp)/se_exp, 1-sqrt(0.5/m))
  # print(t)
  beta_exp_cor <- .fix_point(beta_exp, se_exp, iv.sig, t)
  if (majority == T) {
    iv.valid <- .stablility_selection_stg2_M(beta_exp_cor,
                                             beta_out,
                                             se_exp,
                                             se_out,
                                             iv.sig)
  } else {
    iv.valid <- .stablility_selection_stg2_P(beta_exp_cor,
                                             beta_out,
                                             se_exp,
                                             se_out,
                                             iv.sig)
  }

  if (length(iv.valid)==0) {

    warning("No valid IVs were selected. Estimating causal effect with significant IVs.")
    list(iv.sig, .divw(beta_exp_cor, beta_out, se_exp, se_out, over.dispersion.stg1))

  } else {
    list(iv.valid, .divw(beta_exp_cor, beta_out, se_exp, se_out, over.dispersion.stg2))
  }
}

#' @rdname MRStable
#' @export

ldsc_divw <- function(beta_exp, beta_out, se_exp, se_out, scale_exp, scale_out, over.dispersion = F) {
  m <- length(beta_exp)
  se_exp <- sqrt(scale_exp)*se_exp
  se_out <- sqrt(scale_out*se_out^2)
  .divw(beta_exp, beta_out, se_exp, se_out, over.dispersion)
}

#' @rdname MRStable
#' @export

ada_ldsc_divw <- function (beta_exp, beta_out, se_exp, se_out, scale_exp, scale_out,
          n, maxit = 10000)
{
  m <- length(beta_exp)
  se_exp <- sqrt(scale_exp) * se_exp
  se_out <- sqrt(scale_out) * se_out
  divw.init <- ldsc_divw(beta_exp, beta_out, se_exp, se_out,
                         1, 1)
  beta.init <- divw.init$beta.hat
  alpha.init <- beta_out - beta.init * beta_exp
  w <- 1/(abs(alpha.init)^2 + 1e-16)
  lambda <- 1
  gamma <- rep(0, m)
  for (i in 1:maxit) {
    beta.old <- beta.hat
    alpha.hat <- (abs(beta_out - beta.old * gamma) - lambda * w * se_out^2) *
      (abs(beta_out - beta.old * gamma) > lambda * w * se_out^2) *
      sign(beta_out - beta.old * gamma)
    gamma <- (beta.old * (beta_out - alpha.hat)/se_out^2 +
                beta_exp/se_exp^2)/(beta.old^2/se_out^2 + 1/se_exp^2)
    beta.hat <- sum((beta_out - alpha.hat) * gamma/se_out^2)/sum(gamma^2/se_out^2)
    if (abs(beta.hat - beta.old)/abs(beta.old) < 1e-07) {
      break
    }
  }
  iv.valid <- which(alpha.hat == 0)
  divw.res <- ldsc_divw(beta_exp[iv.valid], beta_out[iv.valid],
                        se_exp[iv.valid], se_out[iv.valid], 1, 1)
  list(beta.hat = divw.res$beta.hat, beta.se = divw.res$beta.se,
       beta.p.value = divw.res$beta.p.value,
       iv.invalid = which(alpha.hat != 0))
}

#' @rdname MRStable
#' @export

ldsc_mcp_divw <- function(beta_exp, beta_out, se_exp, se_out, scale_exp, scale_out, n, a = 3, maxit = 10000) {
  m <- length(beta_exp)
  nlam <- 50
  se_exp <- sqrt(scale_exp)*se_exp
  se_out <- sqrt(scale_out)*se_out
  divw.init <- ldsc_divw(beta_exp, beta_out, se_exp, se_out, 1, 1)
  beta.init <- divw.init$beta.hat
  alpha.init <- beta_out - beta.init*beta_exp
  ny <- min(1/se_out^2)
  lambda <- exp(seq(-0.5*log(n), 0, length.out = nlam))
  bic <- NULL
  alpha.all <- matrix(0, m, nlam)
  for (i in 1:nlam) {
    beta.hat <- beta.init
    alpha.hat <- alpha.init
    for (j in 1:maxit) {
      beta.old <- beta.hat
      gamma <- (beta.old*(beta_out-alpha.hat)/se_out^2+beta_exp/se_exp^2)/(beta.old^2/se_out^2+1/se_exp^2)
      alpha.hat <- (abs(beta_out-beta.old*gamma)-ny*lambda[i]*se_out^2)*
        (abs(beta_out-beta.old*gamma)>ny*lambda[i]*se_out^2)*
        (abs(beta_out-beta.old*gamma)<=a*lambda[i])*
        sign(beta_out-beta.old*gamma)/
        (1-ny/a*se_out^2)+
        (beta_out-beta.old*gamma)*(abs(beta_out-beta.old*gamma)>a*lambda[i])
      beta.hat <- sum((beta_out-alpha.hat)*gamma/se_out^2)/sum(gamma^2/se_out^2)
      if (abs(beta.hat - beta.old)/abs(beta.old) < 1e-7) {
        break
      }
    }
    valid.iv <- which(alpha.hat==0)
    fn_beta <- function(beta) {
      sum((beta_out[valid.iv] - beta*beta_exp[valid.iv])^2/
            (beta^2*se_exp[valid.iv]^2+se_out[valid.iv]^2))
    }
    sol <- optim(beta.init, fn_beta, method = "BFGS")
    ll <- sol$value
    bic <- c(bic, ll+log(n)*sum(alpha.hat!=0))
    alpha.all[, i] <- alpha.hat
  }
  lambda.final <- lambda[which.min(bic)]
  alpha.final <- alpha.all[, which.min(bic)]
  iv.valid <- which(alpha.final==0)
  divw.res <- ldsc_divw(beta_exp[iv.valid],
                        beta_out[iv.valid],
                        se_exp[iv.valid],
                        se_out[iv.valid], 1, 1)
  list(beta.hat = divw.res$beta.hat,
       beta.se = divw.res$beta.se,
       beta.p.value = divw.res$beta.p.value,
       iv.invalid = which(alpha.final!=0))
}

