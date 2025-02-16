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
  list(iv.sig, .divw(beta_exp_cor, beta_out, se_exp, se_out, iv.sig, over.dispersion))
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
    .divw(beta_exp_cor, beta_out, se_exp, se_out, iv.sig, over.dispersion.stg1)

  } else {
    .divw(beta_exp_cor, beta_out, se_exp, se_out, iv.valid, over.dispersion.stg2)
  }
}

