# Choose model

Choose model

## Usage

``` r
average(x, cut_off, avg_method)
```

## Arguments

- x:

  output from `compare_phylosem`

- cut_off:

  threshold where any model with delta-AIC greater than this value is
  excluded from average

- avg_method:

  see
  [`average_DAGs`](https://ax3man.github.io/phylopath/reference/average_DAGs.html)

## Value

Returns an AIC-weighted average of fitted models from
[`compare_phylosem`](https://james-thorson-noaa.github.io/phylosem/reference/compare_phylosem.md)
after conversion to format from \[phylopath::est_DAG\]
