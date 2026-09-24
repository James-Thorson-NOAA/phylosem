# Detailed control for phylosem structure

Define a list of control parameters. Note that the format of this input
is likely to change more rapidly than that of
[`phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/phylosem.md)

## Usage

``` r
phylosem_control(
  nlminb_loops = 1,
  newton_loops = 1,
  trace = 0,
  eval.max = 1000,
  iter.max = 1000,
  getsd = TRUE,
  quiet = FALSE,
  run_model = TRUE,
  getJointPrecision = FALSE,
  profile = c()
)
```

## Arguments

- nlminb_loops:

  Integer number of times to call
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- newton_loops:

  Integer number of Newton steps to do after running
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- trace:

  Parameter values are printed every \`trace\` iteration for the outer
  optimizer. Passed to \`control\` in
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- eval.max:

  Maximum number of evaluations of the objective function allowed.
  Passed to \`control\` in
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- iter.max:

  Maximum number of iterations allowed. Passed to \`control\` in
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- getsd:

  Boolean indicating whether to call
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html)

- quiet:

  Boolean indicating whether to run model printing messages to terminal
  or not;

- run_model:

  Boolean indicating whether to estimate parameters (the default), or
  instead to return the model inputs and compiled TMB object without
  running;

- getJointPrecision:

  whether to get the joint precision matrix. Passed to
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html).

- profile:

  Parameters to profile out of the likelihood (this subset will be
  appended to random with Laplace approximation disabled).

## Value

An S3 object of class "phylosem_control" that specifies detailed model
settings, allowing user specification while also specifying default
values
