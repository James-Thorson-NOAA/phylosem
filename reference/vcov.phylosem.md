# Extract Variance-Covariance Matrix

extract the covariance of fixed effects, or both fixed and random
effects.

## Usage

``` r
# S3 method for class 'phylosem'
vcov(object, which = c("fixed", "random", "both"), ...)
```

## Arguments

- object:

  output from `phylosem`

- which:

  whether to extract the covariance among fixed effects, random effects,
  or both

- ...:

  ignored, for method compatibility
