.stablility_selection_stg1 <- function(beta_exp, se_exp, rd, pi_thr) {
  n <- min(1/se_exp^2)
  nlam <- 10
  lambda <- exp(seq(log(quantile(abs(beta_exp)/se_exp^2, 1-sqrt(0.5/length(beta_exp)))/n), log(max(abs(beta_exp)/se_exp^2)/n), length = nlam))
  stab_sel <- matrix(NA, length(beta_exp), rd)
  for (i in 1:rd) {
    beta_exp_lasso <- matrix(NA, length(beta_exp), nlam)
    W <- runif(length(beta_exp), 0.5, 1)
    beta_exp_dp <- beta_exp + rnorm(length(beta_exp))*se_exp
    for (j in 1:nlam) {
      beta_exp_lasso[, j] <- (abs(beta_exp)>n*lambda[j]*se_exp^2/W)*1
    }
    stab_sel[, i] <- apply(beta_exp_lasso, 1, function(x) (1%in%x))
  }
  which(rowMeans(stab_sel) > pi_thr)
}

.g <- function(x, t) {
  (dnorm(t-x)-dnorm(t+x))/(pnorm(-t-x)+pnorm(-t+x))
}

.fix_point <- function(beta_exp, se_exp, sig_iv, t, maxiter = 500) {
  ratio_hat <- beta_exp[sig_iv]/se_exp[sig_iv]
  ratio_new <- ratio_hat
  for (i in 1:maxiter) {
    ratio_old <- ratio_new
    ratio_new <- ratio_hat - .g(ratio_old, t)
    if (norm(ratio_new-ratio_old, "2")/norm(ratio_old, "2") < 1e-5) {
      break
    }
  }
  beta_exp[sig_iv] <- ratio_new*se_exp[sig_iv]
  beta_exp
}

.stablility_selection_stg2 <- function(beta_exp, beta_out, se_exp, se_out, sig_iv,
                                       dp, pi_thr, maxiter = 500) {
  m <- length(sig_iv)
  n <- min(c(1/se_exp^2, 1/se_out^2))
  beta_exp <- beta_exp[sig_iv]
  beta_out <- beta_out[sig_iv]
  se_exp <- se_exp[sig_iv]
  se_out <- se_out[sig_iv]
  invalid_set <- matrix(NA, m, dp)
  for (i in 1:dp) {
    beta_exp_dp <- beta_exp + rnorm(m)*se_exp
    beta_out_dp <- beta_out + rnorm(m)*se_out
    alpha_cml <- matrix(NA, m, m-1)
    bic <- NULL
    for (K in 1:(m-1)) {
      beta.hat <- 0
      gamma <- rep(0, m)
      for (j in 1:maxiter) {
        beta.old <- beta.hat
        Q <- (beta_out_dp - beta.old*gamma)^2/se_out^2
        nonzero_index <- order(Q, decreasing = T)[1:K]
        alpha <- rep(0, m)
        alpha[nonzero_index] <- (beta_out_dp - beta.old*gamma)[nonzero_index]
        gamma <- ((beta_out_dp - alpha)*beta.old/se_out^2+beta_exp_dp/se_exp^2)/
          (beta.old^2/se_out^2+1/se_exp^2)
        beta.hat <- sum((beta_out_dp - alpha)*gamma/se_out^2)/sum(gamma^2/se_out^2)
        if (abs(beta.hat - beta.old)/abs(beta.old+1e-16) < 1e-7) {
          break
        }
      }
      gamma[nonzero_index] <- beta_exp_dp[nonzero_index]
      alpha[nonzero_index] <- (beta_out_dp - beta.hat*gamma)[nonzero_index]
      alpha_cml[, K] <- alpha
      bic <- c(bic, -2*sum(dnorm(beta_out_dp, beta.hat*gamma+alpha, se_out, log = T)+
                             dnorm(beta_exp_dp, gamma, se_exp, log = T)) + log(n)*K)
    }
    invalid_set[, i] <- (alpha_cml[, which.min(bic)]!=0)*1
  }
  sig_iv[which(rowMeans(invalid_set) < pi_thr)]
}

.divw <- function(beta_exp, beta_out, se_exp, se_out, ind, over.dispersion) {
  beta.hat <- sum(beta_exp[ind] * beta_out[ind]/(se_out[ind])^2)/
    sum((beta_exp[ind]^2 - se_exp[ind]^2)/(se_out[ind])^2)
  se.ratio <- se_exp[ind]/se_out[ind]
  mu <- beta_exp[ind]/se_exp[ind]
  tau.square <- 0
  if (over.dispersion) {
    tau.square <- sum(((beta_out[ind] - beta.hat * beta_exp[ind])^2 - se_exp[ind]^2 - beta.hat^2 * se_exp[ind]^2)/se_out[ind]^2)/
      sum(se_out[ind]^(-2))
  }
  V1 <- sum(se.ratio^2 * mu^2 + beta.hat^2 * se.ratio^4 *
              (mu^2 + 1) + tau.square * se.ratio^2/se_out[ind]^2 *
              mu^2)
  V2 <- sum(se.ratio^2 * (mu^2 - 1))
  beta.var <- V1/V2^2
  beta.se <- sqrt(beta.var)
  list(beta.hat = beta.hat,
       beta.se = beta.se,
       beta.p.value = 2*pnorm(-abs(beta.hat)/beta.se))
}
