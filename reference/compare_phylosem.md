# Compare phylogenetic structural equation models

Fits several phylogenetic structural equation model for further
comparison

## Usage

``` r
compare_phylosem(
  sem_set,
  tree,
  data,
  family = rep("fixed", ncol(data)),
  covs,
  estimate_ou = FALSE,
  estimate_lambda = FALSE,
  estimate_kappa = FALSE,
  control = phylosem_control(),
  ...
)
```

## Arguments

- sem_set:

  A named list of structural equation model specifications, where each
  element will be passed as argument `sem` to
  [`phylosem`](https://james-thorson-noaa.github.io/phylosem/reference/phylosem.md)

- tree:

  phylogenetic structure, using class
  [`as.phylo`](https://rdrr.io/pkg/ape/man/as.phylo.html)

- data:

  data-frame providing variables being modeled. Missing values are
  inputted as NA. If an SEM includes a latent variable (i.e., variable
  with no available measurements) then it still must be inputted as a
  column of `data` with entirely NA values.

- family:

  Character-vector listing the distribution used for each column of
  `data`, where each element must be `fixed`, `normal`, `binomial`, or
  `poisson`. `family="fixed"` is default behavior and assumes that a
  given variable is measured exactly. Other options correspond to
  different specifications of measurement error.

- covs:

  optional: a character vector of one or more elements, with each
  element giving a string of variable names, separated by commas.
  Variances and covariances among all variables in each such string are
  added to the model. For confirmatory factor analysis models specified
  via `cfa`, `covs` defaults to all of the factors in the model, thus
  specifying all variances and covariances among these factors.
  *Warning*: `covs="x1, x2"` and `covs=c("x1", "x2")` are *not*
  equivalent: `covs="x1, x2"` specifies the variance of `x1`, the
  variance of `x2`, *and* their covariance, while `covs=c("x1", "x2")`
  specifies the variance of `x1` and the variance of `x2` *but not*
  their covariance.

- estimate_ou:

  Boolean indicating whether to estimate an autoregressive
  (Ornstein-Uhlenbeck) process using additional parameter `lnalpha`,
  corresponding to the `model="OUrandomRoot"` parameterization from
  phylolm as listed in
  [doi:10.1093/sysbio/syu005](https://doi.org/10.1093/sysbio/syu005)

- estimate_lambda:

  Boolean indicating whether to estimate additional branch lengths for
  phylogenetic tips (a.k.a. the Pagel-lambda term) using additional
  parameter `logitlambda`

- estimate_kappa:

  Boolean indicating whether to estimate a nonlinear scaling of branch
  lengths (a.k.a. the Pagel-kappa term) using additional parameter
  `lnkappa`

- control:

  Output from
  [`phylosem_control`](https://james-thorson-noaa.github.io/phylosem/reference/phylosem_control.md),
  used to define user settings, and see documentation for that function
  for details.

- ...:

  Additional arguments passed to
  [`phylosem`](https://james-thorson-noaa.github.io/phylosem/reference/phylosem.md)

## Value

An object (list) of class \`compare_phylosem\`, containing a list of
output from
[`phylosem`](https://james-thorson-noaa.github.io/phylosem/reference/phylosem.md)
