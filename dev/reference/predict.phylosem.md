# predict phylosem

extract the predictive mean and covariance for nodes or tips based on
names from tree, `x$tree$tip.label` and `x$tree$node.label`

## Usage

``` r
# S3 method for class 'phylosem'
predict(
  object,
  taxon,
  what = c("p_g"),
  se.fit = TRUE,
  include_fixed = TRUE,
  ...
)
```

## Arguments

- object:

  Output from
  [`phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/phylosem.md)

- taxon:

  vector of taxon names (matched against tree) to predict traits

- what:

  What output to calculate (currently only the linear predictor is
  supported)

- se.fit:

  Whether to report predictive covariance, where `sqrt(diag(pred$cov))`
  is the standard errors

- include_fixed:

  Whether to include uncertainty in fixed effects when calculating
  predictive covariance

- ...:

  Not used

## Value

When `se.fit = FALSE`, it returns a matrix of predicted traits (rows) by
taxon (column). When `se.fit = TRUE`, it returns a list:

- est:

  matrix of predicted traits

- cov:

  array of predictive covariance among traits (1st and 2nd dims) by
  taxon (3rd dim)
