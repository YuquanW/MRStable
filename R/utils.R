.stablility_selection_stg1 <- function(beta, se, rd, pi_thr) {
  n <- min(1/se^2)
  nlam <- 20
  lambda <- exp(seq(log(quantile(abs(beta)/se^2, 1-sqrt(0.5/length(beta)))/n), log(max(abs(beta)/se^2)/n), length = nlam))
  stab_sel <- matrix(NA, length(beta), rd)
  for (i in 1:rd) {
    beta_lasso <- matrix(NA, length(beta), nlam)
    W <- runif(length(beta), 0.75, 1)
    beta_dp <- beta + rnorm(length(beta))*se
    for (j in 1:nlam) {
      beta_lasso[, j] <- (abs(beta_dp)>n*lambda[j]*se^2/W)*1
    }
    stab_sel[, i] <- apply(beta_lasso, 1, max)
  }
  which(rowMeans(stab_sel) > pi_thr)
}

# .stablility_selection_stg1 <- function(beta_exp, se_exp, rd, pi_thr) {
#   n <- min(1/se_exp^2)
#   nlam <- 100
#   lambda <- exp(seq(log(quantile(abs(beta_exp)/se_exp^2, 0.1)/n), log(max(abs(beta_exp)/se_exp^2)/n), length = nlam))
#   print(lambda[1])
#   beta_exp_lasso <- matrix(NA, length(beta_exp), nlam)
#   bic <- NULL
#   for (j in 1:nlam) {
#     beta_exp_lasso[, j] <- (abs(beta_exp) - n*lambda[j]*se_exp^2)*
#       (abs(beta_exp)>n*lambda[j]*se_exp^2)*sign(beta_exp)
#     bic <- c(bic, sum((beta_exp-beta_exp_lasso[, j])^2/se_exp^2)+log(n)*sum(beta_exp_lasso[, j]!=0))
#   }
#   beta_exp_lasso1 <- beta_exp_lasso[, which.min(bic)]
#   plot(lambda, bic)
#   #print(quantile(abs(beta_exp)*(abs(beta_exp_lasso1)+1e-16)/se_exp^2, 0.1)/n)
#   lambda <- exp(seq(log(quantile(abs(beta_exp)*(abs(beta_exp_lasso1)+1e-16)/se_exp^2, 0.1)/n), log(max(abs(beta_exp)*(abs(beta_exp_lasso1)+1e-16)/se_exp^2)/n), length = nlam))
#
#   beta_exp_lasso <- matrix(NA, length(beta_exp), nlam)
#   bic <- NULL
#   for (j in 1:nlam) {
#     beta_exp_lasso[, j] <- (abs(beta_exp) - n*lambda[j]*se_exp^2/(abs(beta_exp_lasso1)+1e-16))*
#       (abs(beta_exp)>n*lambda[j]*se_exp^2/(abs(beta_exp_lasso1)+1e-16))*sign(beta_exp)
#     bic <- c(bic, sum((beta_exp-beta_exp_lasso[, j])^2/se_exp^2)+log(n)*sum(beta_exp_lasso[, j]!=0))
#   }
#   plot(lambda, bic)
#   #print(lambda[which.min(bic)])
#   beta_exp_lasso2 <- beta_exp_lasso[, which.min(bic)]
#   which(beta_exp_lasso2!=0)
# }

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

# .stablility_selection_stg2 <- function(beta_exp, beta_out, se_exp, se_out, sig_iv,
#                                        dp, pi_thr, maxiter = 500) {
#   m <- length(sig_iv)
#   n <- min(c(1/se_exp^2, 1/se_out^2))
#   beta_exp <- beta_exp[sig_iv]
#   beta_out <- beta_out[sig_iv]
#   se_exp <- se_exp[sig_iv]
#   se_out <- se_out[sig_iv]
#
#   valid_set <- matrix(NA, m, dp)
#   for (i in 1:dp) {
#     beta_exp_dp <- beta_exp + rnorm(m, sd = 0.2)*se_exp
#     beta_out_dp <- beta_out + rnorm(m, sd = 0.2)*se_out
#     alpha_cml <- matrix(NA, m, m-1)
#     rbic <- NULL
#     W <- runif(m, 0.8, 1)
#     for (K in 1:(m-1)) {
#       beta.hat <- 0
#       gamma <- rep(0, m)
#       for (j in 1:maxiter) {
#         beta.old <- beta.hat
#         Q <- (beta_out_dp - beta.old*gamma)^2/se_out^2
#         nonzero_index <- order(Q, decreasing = T)[1:K]
#         alpha <- rep(0, m)
#         alpha[nonzero_index] <- (beta_out_dp - beta.old*gamma)[nonzero_index]
#         gamma <- ((beta_out_dp - alpha)*beta.old/se_out^2+beta_exp_dp/se_exp^2)/
#           (beta.old^2/se_out^2+1/se_exp^2)
#         beta.hat <- sum((beta_out_dp - alpha)*gamma/se_out^2)/sum(gamma^2/se_out^2)
#         if (abs(beta.hat - beta.old)/abs(beta.old+1e-16) < 1e-7) {
#           break
#         }
#       }
#       gamma[nonzero_index] <- beta_exp_dp[nonzero_index]
#       alpha[nonzero_index] <- (beta_out_dp - beta.hat*gamma)[nonzero_index]
#       alpha_cml[, K] <- alpha
#       rbic <- c(rbic, sum((beta_out_dp - beta.hat*gamma - alpha)^2/(1.04*se_out^2)+
#                             (beta_exp_dp - gamma)^2/(1.04*se_exp^2)+W*log(n)*(alpha!=0)))
#     }
#     valid_set[, i] <- (alpha_cml[, which.min(rbic)]==0)*1
#   }
#   sig_iv[which(rowMeans(valid_set) > pi_thr)]
# }

# .stablility_selection_stg2 <- function(beta_exp, beta_out, se_exp, se_out, sig_iv,
#                                        dp, pi_thr, maxiter = 500) {
#   m <- length(sig_iv)
#   n <- min(c(1/se_exp^2, 1/se_out^2))
#   beta_exp <- beta_exp[sig_iv]
#   beta_out <- beta_out[sig_iv]
#   se_exp <- se_exp[sig_iv]
#   se_out <- se_out[sig_iv]
#   invalid_set <- matrix(NA, m, dp)
#   nlam <- 10
#   beta.ratio <- beta_out/beta_exp
#   se.ratio <- sqrt(se_out^2/beta_exp^2+beta_out^2*se_exp^2/beta_exp^4)
#   dens <- .hetero_kde(beta.ratio, se.ratio)
#   beta.init <- dens$x[which.max(dens$y)]
#
#   lambda <- exp(seq(log(min(abs(beta_out-beta.init*beta_exp)/se_out^2)/n), log(quantile(abs(beta_out-beta.init*beta_exp)/se_out^2, 0.5)/n), length = nlam))
#   for (i in 1:dp) {
#     beta_exp_dp <- beta_exp + rnorm(m, sd = 0.2)*se_exp
#     beta_out_dp <- beta_out + rnorm(m, sd = 0.2)*se_out
#     alpha_lasso <- matrix(NA, m, nlam)
#     W <- runif(m, 1, 1.2)
#     for (j in 1:nlam) {
#       beta.hat <- beta.init
#       gamma <- beta_exp_dp
#       for (k in 1:maxiter) {
#         beta.old <- beta.hat
#         alpha <- (abs(beta_out_dp - beta.old*gamma) - n*lambda[j]*se_out^2/W)*
#           sign(beta_out_dp - beta.old*gamma)*
#           (abs(beta_out_dp - beta.old*gamma)>n*lambda[j]*se_out^2/W)
#         gamma <- ((beta_out_dp - alpha)*beta.old/se_out^2+beta_exp_dp/se_exp^2)/
#           (beta.old^2/se_out^2+1/se_exp^2)
#         beta.hat <- sum((beta_out_dp - alpha)*gamma/se_out^2)/sum(gamma^2/se_out^2)
#         if (abs(beta.hat - beta.old)/abs(beta.old+1e-16) < 1e-5) {
#           break
#         }
#       }
#       alpha_lasso[, j] <- (alpha!=0)*1
#       #dnorm(beta_exp_dp, gamma, se_exp, log = T)) + log(n)*K)
#     }
#     invalid_set[, i] <- apply(alpha_lasso, 1, min)
#   }
#   sig_iv[which(rowMeans(invalid_set) < pi_thr)]
# }

.stablility_selection_stg2_M <- function(beta_exp, beta_out, se_exp, se_out, sig_iv, maxiter = 500) {
  m <- length(sig_iv)
  n <- min(c(1/se_exp^2, 1/se_out^2))
  beta_exp <- beta_exp[sig_iv]
  beta_out <- beta_out[sig_iv]
  se_exp <- se_exp[sig_iv]
  se_out <- se_out[sig_iv]
  nlam <- 50
  beta.ratio <- beta_out/beta_exp
  beta.init <- median(beta.ratio)
  alpha.init <- beta_out - beta.init*beta_exp
  lambda <- exp(seq(log(quantile(abs(beta_out-beta.init*beta_exp)*abs(alpha.init)/se_out^2, 0.1)/n),
                    log(max(abs(beta_out-beta.init*beta_exp)*abs(alpha.init)/se_out^2)/n), length = nlam))
  bic <- NULL
  alpha_lasso <- matrix(NA, m, nlam)
  for (j in 1:nlam) {
    beta.hat <- 0
    gamma <- beta_exp
    for (k in 1:maxiter) {
      beta.old <- beta.hat
      alpha <- (abs(beta_out - beta.old*gamma) - n*lambda[j]*se_out^2/(abs(alpha.init)+1e-16))*
        sign(beta_out - beta.old*gamma)*
        (abs(beta_out - beta.old*gamma)>n*lambda[j]*se_out^2/abs(alpha.init))
      gamma <- ((beta_out - alpha)*beta.old/se_out^2+beta_exp/se_exp^2)/
        (beta.old^2/se_out^2+1/se_exp^2)
      beta.hat <- sum((beta_out - alpha)*gamma/se_out^2)/sum(gamma^2/se_out^2)
      if (abs(beta.hat - beta.old)/abs(beta.old+1e-16) < 1e-5) {
        break
      }
    }
    alpha_lasso[, j] <- alpha
    bic <- c(bic, sum((beta_out - beta.hat*gamma - alpha)^2/(se_out^2)+
                        (beta_exp - gamma)^2/(se_exp^2)+log(n)*(alpha!=0)))
  }
  #plot(lambda, bic)
  sig_iv[which(alpha_lasso[, which.min(bic)]==0)]
}

.stablility_selection_stg2_P <- function(beta_exp, beta_out, se_exp, se_out, sig_iv, maxiter = 500) {
  m <- length(sig_iv)
  n <- min(c(1/se_exp^2, 1/se_out^2))
  beta_exp <- beta_exp[sig_iv]
  beta_out <- beta_out[sig_iv]
  se_exp <- se_exp[sig_iv]
  se_out <- se_out[sig_iv]
  nlam <- 50
  beta.ratio <- beta_out/beta_exp
  se.ratio <- median(sqrt(se_out^2/beta_exp^2+beta_out^2*se_exp^2/beta_exp^4))
  dens <- .hetero_kde(beta.ratio, se.ratio)
  #plot(dens$x, dens$y)
  beta.init <- dens$x[which.max(dens$y)]
  #print(beta.init)
  alpha.init <- beta_out - beta.init*beta_exp
  lambda <- exp(seq(log(quantile(abs(beta_out-beta.init*beta_exp)*abs(alpha.init)/se_out^2, 0.1)/n),
                    log(max(abs(beta_out-beta.init*beta_exp)*abs(alpha.init)/se_out^2)/n), length = nlam))
  bic <- NULL
  alpha_lasso <- matrix(NA, m, nlam)
  for (j in 1:nlam) {
    beta.hat <- 0
    gamma <- beta_exp
    for (k in 1:maxiter) {
      beta.old <- beta.hat
      alpha <- (abs(beta_out - beta.old*gamma) - n*lambda[j]*se_out^2/(abs(alpha.init)+1e-16))*
        sign(beta_out - beta.old*gamma)*
        (abs(beta_out - beta.old*gamma)>n*lambda[j]*se_out^2/abs(alpha.init))
      gamma <- ((beta_out - alpha)*beta.old/se_out^2+beta_exp/se_exp^2)/
        (beta.old^2/se_out^2+1/se_exp^2)
      beta.hat <- sum((beta_out - alpha)*gamma/se_out^2)/sum(gamma^2/se_out^2)
      if (abs(beta.hat - beta.old)/abs(beta.old+1e-16) < 1e-5) {
        break
      }
    }
    alpha_lasso[, j] <- alpha
    bic <- c(bic, sum((beta_out - beta.hat*gamma - alpha)^2/(se_out^2)+
                        (beta_exp - gamma)^2/(se_exp^2)+log(n)*(alpha!=0)))
  }
  #plot(lambda, bic)
  sig_iv[which(alpha_lasso[, which.min(bic)]==0)]
}

.divw <- function(beta_exp, beta_out, se_exp, se_out, ind, over.dispersion) {
  beta.hat <- sum(beta_exp[ind] * beta_out[ind]/(se_out[ind])^2)/
    sum((beta_exp[ind]^2 - se_exp[ind]^2)/(se_out[ind])^2)
  se.ratio <- se_exp[ind]/se_out[ind]
  mu <- beta_exp[ind]/se_exp[ind]
  tau.square <- 0
  if (over.dispersion) {
    tau.square <- max(sum(((beta_out[ind] - beta.hat * beta_exp[ind])^2 - se_exp[ind]^2 - beta.hat^2 * se_exp[ind]^2)/se_out[ind]^2)/
      sum(se_out[ind]^(-2)), 0)
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

.hetero_kde <- function(beta, se, res = 1000) {
  x_eval <- seq(
    min(beta) - 3 * max(se),
    max(beta) + 3 * max(se),
    length.out = res
  )
  dens <- sapply(x_eval, function(x) {
    sum(dnorm(x, mean = beta, sd = se))
  }) / length(beta)
  list(x = x_eval, y = dens)
}
