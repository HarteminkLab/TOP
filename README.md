
<!-- README.md is generated from README.Rmd. Please edit that file -->
TOP: Transcription factor Occupancy Prediction
==============================================

<!-- badges: start -->
<!-- badges: end -->
TOP is a Bayesian hierarchical model trained using DNase-seq and ChIP-seq data to Predict transcription factor (TF) occupancy for multiple TFs across multiple cell types.

Install `TOP` R package
-----------------------

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("kevinlkx/TOP")
library(TOP)
```

Predict TF occupancy using trained `TOP` model
----------------------------------------------

### Prepare motif and DNase data

For each TF in each cell type, prepare data for candidate binding sites, including: PWM scores for motif matches, and DNase bins (using MILLIPEDE binning).

We provided a [pipeline](../doc/preprocessing.html) to obtain candidate binding sites and get DNase bins using DNase-seq bam files. You can also use your own pipelines to obtain motif matches, and DNase matrices (and normalize by library size, etc.), and use our function `millipede_bin_dnase` or `MILLIPEDE` software to bin DNase data with `MILLIPEDE` binning scheme.

### Predict TF occupancy using trained `TOP` model

You can [download](../data/) our pre-trained `TOP` models for Duke or UW DNase data, or use models trained from your own data.

Option1: Predict TF occupancy using posterior samples of regression coefficients trained from the TOP model.

``` r
# data: data frame or matrix. Columns are motif score and DNase features. Rows are candidate sites.
# alpha_samples: posterior samples of alpha parameters
# beta_samples: posterior samples of beta parameters
# tau_samples: posterior samples of tau parameters
# sample: logicals. If TRUE, samples from posterior predictions and then take the mean of posterior prediction samples.
# average_parameters: logicals. If TRUE, uses the posterior mean of regression coefficients to make predictions.
# transform: method used to transform ChIP-seq counts when training the TOP model. Default: transform = "asinh".
predict_TOP(data, alpha_samples, beta_samples, tau_samples, sample = TRUE, average_parameters = FALSE, transform = 'asinh')
```

Option 2: Predict TF occupancy using posterior mean of regression coefficients trained from the TOP model.

``` r
# data: data frame or matrix. Columns are motif score and DNase features. Rows are candidate sites.
# coef_mean: numeric vector. The posterior mean of trained regression coefficients.
# transform: method used to transform ChIP-seq counts when training the TOP model. Default: transform = "asinh".
predict_TOP_coef_mean(data, top_coef_mean, transform = 'asinh')
```

Train `TOP` model using motif, DNase and ChIP data from all TF-cell type combinations
-------------------------------------------------------------------------------------

### Prepare motif, DNase and ChIP training data

You will need to prepare training data for each TF in each cell type, including: PWM scores, normalized DNase data (with MILLIPEDE binning), and normalized ChIP read counts, and then combine training datasets for all training TF-cell type combinations.

### Train `TOP` model using combined training data from all TF-cell type combinations

`TOP` uses MCMC to fit the Bayesian hierarchical model. The current version was implemented using [JAGS](http://mcmc-jags.sourceforge.net). You will need `rjags` R package to fit the model.

``` r
library(rjags)
library(TOP)

# data.train: combined training data.
# TOP.model: TOP model written in BUGS code.
# parameters.to.save: character vector of the names of the parameters to save which should be monitored.
# n.iter number of total iterations per chain (including burn in).
# n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
# n.thin thinning rate, must be a positive integer (default=2).
# n.chains number of Markov chains (default: 1).
fit_TOP_jags(data.train, TOP.model, parameters.to.save, n.iter, n.burnin, n.thin, n.chains)
```
