# Compare phylogenetic structural equation models

Fits several phylogenetic structural equation model for further
comparison

## Usage

``` r
compare_phylosem(
  sem_set,
  tree,
  data,
  family = Map(function(.) fixed(), colnames(data)),
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
  [`phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/phylosem.md)

- tree:

  phylogenetic structure, using class
  [`as.phylo`](https://rdrr.io/pkg/ape/man/as.phylo.html)

- data:

  data-frame providing numeric values for variables being modeled.
  Missing values are inputted as NA. If an SEM includes a latent
  variable (i.e., variable with no available measurements) then it still
  must be inputted as a column of `data` with entirely NA values.
  Bernoulli variables must be coded as 0s or 1s, and factors are not
  allowed.

- family:

  A named list of families, each returning a class `family`, including
  \[fixed()\], \[gaussian()\], \[binomial()\], \[Gamma()\], and
  \[poisson()\], with names that match levels of `colnames(data)` to
  allow different families by variable. Family \[fixed()\] specifies
  that states are known (i.e., measurements for that variable have no
  error). Other families allow users to supply a link function including
  \`identity\`, \`log\`, \`logit\`, or \`cloglog\`. For example
  `family = list(y = binomial("logit"), x = fixed())` would specify
  logit-linked Bernoulli distribution for variable \`data\$y\` and a
  fixed (no measurement error) distribution for \`data\$x\`. For many
  variables, it is convenient to do e.g.,
  `family = Map(function(.) gaussian(), colnames(tsdata))` rather than
  writing them all manually.

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
  [`phylosem_control`](https://james-thorson-noaa.github.io/phylosem/dev/reference/phylosem_control.md),
  used to define user settings, and see documentation for that function
  for details.

- ...:

  Additional arguments passed to
  [`phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/phylosem.md)

## Value

An object (list) of class \`compare_phylosem\`, containing a list of
output from
[`phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/phylosem.md)
