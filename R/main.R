#' Title
#'
#' @param beta_exp Exposure
#' @param beta_out Outcome
#' @param se_exp Std error exposure
#' @param se_out Std error outcome
#' @param pi_thr1 Threshold 1
#' @param pi_thr2 Threshold 2
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
                     pi_thr1 = 0.6,
                     rd = 20,
                     over.dispersion = T) {
  m <- length(beta_exp)
  iv.sig <- .stablility_selection_stg1(beta_exp,
                                       se_exp,
                                       rd,
                                       pi_thr1)
  t <- quantile(abs(beta_exp)/se_exp, 1-sqrt(0.5/m))
  beta_exp_cor <- .fix_point(beta_exp, se_exp, iv.sig, t)
  .divw(beta_exp_cor, beta_out, se_exp, se_out, iv.sig, over.dispersion)
}

#' @rdname MRStable
#' @export

MRStable_V <- function(beta_exp,
                       beta_out,
                       se_exp,
                       se_out,
                       pi_thr1 = 0.6,
                       pi_thr2 = 0.4,
                       rd = 20,
                       dp = 20,
                       over.dispersion = T) {
  m <- length(beta_exp)
  iv.sig <- .stablility_selection_stg1(beta_exp,
                                       se_exp,
                                       rd,
                                       pi_thr1)
  t <- quantile(abs(beta_exp)/se_exp, 1-sqrt(0.5/m))
  beta_exp_cor <- .fix_point(beta_exp, se_exp, iv.sig, t)
  iv.valid <- .stablility_selection_stg2(beta_exp_cor,
                                         beta_out,
                                         se_exp,
                                         se_out,
                                         iv.sig,
                                         dp,
                                         pi_thr2)
  .divw(beta_exp_cor, beta_out, se_exp, se_out, iv.valid, over.dispersion)
}

