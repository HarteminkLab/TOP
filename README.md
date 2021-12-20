
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TOP: Transcription factor Occupancy Prediction

<!-- badges: start -->
<!-- badges: end -->

TOP fits a Bayesian hierarchical model using motif information,
DNase-seq and ChIP-seq data from multiple TFs in multiple cell types.
Then, it can predict the quantitative occupancy of multiple
transcription factors (TFs) using data from a single DNase-seq
experiment.

# Install TOP R package

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("kevinlkx/TOP")
library(TOP)
```

Transcription factor Occupancy Profiler (TOP) fits a Bayesian
hierarchical model using transcription factor (TF) motif information,
DNase or ATAC-seq and ChIP-seq data from multiple TFs across multiple
cell types or conditions.

It can be used to predict the quantitative occupancy (or TF binding
probability) for many TFs using data from a single DNase- or ATAC-seq
experiment.

# Predict TF occupancy using TOP model

Here are some instructions on how to [predict TF occupancy using TOP
model](articles/predict_TF_occupancy_with_trained_model.html).

We will provide pre-trained models using ENCODE data
[here](https://users.cs.duke.edu/~amink/software/)..

# Train your own TOP model

You can also train TOP models using your own data.

Here are some instructions on how to [train a TOP
model](articles/train_TOP_model_JAGS.html).

It takes some time to train the model for many TFs and cell types, but
the training only needs to be done once. Once trained, you can use the
model to make predictions for new data sets.

### Reference

-   Luo, K., Zhong, J., Safi, A., Hong, L., Tewari, A., Song, L., Reddy,
    T., Ma, L., Crawford, G., & Hartemink, A. (2020) “Quantitative
    occupancy of myriad transcription factors from one DNase experiment
    enables efficient comparisons across conditions.” bioRxiv,
    bioRxiv:2020.06.28.171587.
