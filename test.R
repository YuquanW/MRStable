dat.exposure <- read.table("clumped_HDL_exp.txt", sep = "\t", header = T)
dat.outcome <- read.table("clumped_HDL_out.txt", sep = "\t", header = T)
dat.exposure$exposure <- "None"
dat.exposure$id.exposure <- "None"
dat.outcome$outcome <- "None"
dat.outcome$id.outcome <- "None"
dat <- TwoSampleMR::harmonise_data(dat.exposure, dat.outcome)
plot(dat$beta.exposure, dat$beta.outcome)
dat <- subset(dat, pval.exposure <= 0.05)

png("fig3-1.png", width = 2800, height = 2100, res = 300)
par(mai=c(1, 1.2, 1, 1))
plot(dat$se.exposure, -log10(dat$pval.exposure),
     xlab = expression(sigma[X[j]]),
     ylab = expression("-"*log[10]*"p"),
     pch = 16,
     cex = 2,
     cex.axis = 2,
     cex.lab = 2)
dev.off()


mr.divw::mr.divw(dat$beta.exposure, dat$beta.outcome,
                 dat$se.exposure,
                 dat$se.outcome, over.dispersion = T)
.divw(dat$beta.exposure, dat$beta.outcome,
      dat$se.exposure,
      dat$se.outcome,T)
mr.raps::mr.raps.mle.all(dat$beta.exposure,
                         dat$beta.outcome,
                         sqrt(dat$scale.exposure*dat$se.exposure^2),
                         sqrt(dat$scale.outcome*dat$se.outcome^2))
ldsc_divw(dat$beta.exposure,
          dat$beta.outcome,
          dat$se.exposure,
          dat$se.outcome,
          dat$scale.exposure,
          dat$scale.outcome,
          over.dispersion = F)
ldsc_mcp_divw(dat$beta.exposure,
              dat$beta.outcome,
              dat$se.exposure,
              dat$se.outcome,
              dat$scale.exposure,
              dat$scale.outcome,
              dat$samplesize.exposure,
              dat$samplesize.outcome,
              over.dispersion = F)

gamma <- dat$beta.exposure
beta_exp_simu <- gamma + rnorm(nrow(dat), 0, 1)*dat$se.exposure
alpha <- c(rnorm(floor(0.3*nrow(dat)), 0, 10*mean(dat$se.outcome)),
           rep(0, nrow(dat) - floor(0.3*nrow(dat))))
#delta <- rnorm(nrow(dat), 0, 5*mean(dat$se.outcome))
beta_out_simu <- gamma*0 + alpha + rnorm(nrow(dat), 0, 1)*dat$se.outcome
plot(beta_exp_simu, beta_out_simu)
ldsc_divw(beta_exp_simu,
          beta_out_simu,
          dat$se.exposure,
          dat$se.outcome,
          1,
          1,
          T)
ldsc_mcp_divw(beta_exp_simu,
              beta_out_simu,
              dat$se.exposure,
              dat$se.outcome,
              1,
              1,
              dat$samplesize.exposure, dat$samplesize.outcome,
              over.dispersion = T)
fn <- function(beta) {
  sum((dat$beta.outcome-beta*dat$beta.exposure)^2/(beta^2*dat$se.exposure^2+dat$se.outcome^2))
}
fn <- function(beta) {
  sum((dat$beta.outcome-beta*dat$beta.exposure)^2/dat$se.outcome^2 - beta^2*dat$se.exposure^2/dat$se.outcome^2)
}

plot(seq(-10, 10, 0.01), Map(fn, seq(-10, 10, 0.01)))
