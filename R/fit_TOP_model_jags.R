## TOP model in BUGS code for training multiple TFs in multiple cell types
## It can be extended to model the heterogeneity across replicates
TOP.model <- function() {
  for(i in 1:N) {
    chip[i] ~ dnorm(mu[i],  tau[tf[i], cell_type[i]])

    mu[i] <- alpha[tf[i], cell_type[i]] +
      pwm[i] * beta1[tf[i], cell_type[i]] +
      dnase.left2_sum[i] * beta2[tf[i], cell_type[i]] +
      dnase.left1_sum[i] * beta3[tf[i], cell_type[i]] +
      dnase.motif_sum[i] * beta4[tf[i], cell_type[i]] +
      dnase.right1_sum[i] * beta5[tf[i], cell_type[i]] +
      dnase.right2_sum[i] * beta6[tf[i], cell_type[i]]
  }

  for(j in 1:n_tfs) {
    for(k in 1:n_cell_types) {

      alpha[j, k] ~ dnorm(Alpha[j], 1)
      beta1[j, k] ~ dnorm(Beta1[j], 1)
      beta2[j, k] ~ dnorm(Beta2[j], 1)
      beta3[j, k] ~ dnorm(Beta3[j], 1)
      beta4[j, k] ~ dnorm(Beta4[j], 1)
      beta5[j, k] ~ dnorm(Beta5[j], 1)
      beta6[j, k] ~ dnorm(Beta6[j], 1)
      tau[j, k] ~ dgamma(T[j]^2, T[j])

    }
  }

  for(i in 1:n_tfs) {
    Alpha[i] ~ dnorm(A, 1)
    Beta1[i] ~ dnorm(B1, 1)
    Beta2[i] ~ dnorm(B2, 1)
    Beta3[i] ~ dnorm(B3, 1)
    Beta4[i] ~ dnorm(B4, 1)
    Beta5[i] ~ dnorm(B5, 1)
    Beta6[i] ~ dnorm(B6, 1)
    T[i] ~ dgamma(TAU^2, TAU)
  }

  A ~ dnorm(0, 1)
  B1 ~ dnorm(0, 1)
  B2 ~ dnorm(0, 1)
  B3 ~ dnorm(0, 1)
  B4 ~ dnorm(0, 1)
  B5 ~ dnorm(0, 1)
  B6 ~ dnorm(0, 1)
  TAU ~ dgamma(1, 1)
}


#' @title Fit TOP model using JAGS with R2jags package
#'
#' @param data.train combined training data.
#' @param TOP.model TOP model written in BUGS code.
#' @param parameters.to.save character vector of the names of the parameters to
#' save which should be monitored.
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.thin thinning rate, must be a positive integer.
#' @param n.chains number of Markov chains.
#' @param DIC logical; if TRUE (default), compute deviance, pD, and DIC.
#'
#' @return
#' @export
fit_TOP_jags <- function(data.train, TOP.model, parameters.to.save,
                         n.iter=1e5, n.burnin=5000, n.thin=2, n.chains=1, DIC=TRUE) {

  attach(data.train)

  N <- nrow(data.train)
  n_tfs <- length(unique(data.train$tf_id))
  n_cell_types <- length(unique(data.train$cell_type_id))

  jags.data <- list('pwm',
                    'dnase.left2_sum', 'dnase.left1_sum', 'dnase.motif_sum', 'dnase.right1_sum', 'dnase.right2_sum',
                    'chip_count',
                    'cell_type_id', 'tf_id',
                    'N', 'n_tfs', 'n_cell_types')

  ## using R2jags
  jags_fit <- jags(data = jags.data, parameters.to.save, TOP.model,
                  n.iter, n.burnin, n.thin, n.chains, DIC)

  return(jags_fit)

}


