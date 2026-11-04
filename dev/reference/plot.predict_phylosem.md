# plot predictions from phylosem

plot the predictive distribution for one or two variables extracted
using
[`predict.phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/predict.phylosem.md)

## Usage

``` r
# S3 method for class 'predict_phylosem'
plot(x, y, variables = 1, data, data_labels, distance_threshold = 0, tree, ...)
```

## Arguments

- x:

  Output from
  [`predict.phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/predict.phylosem.md)

- y:

  Not used

- variables:

  which variables to plot (can only plot one or two at a time)

- data:

  optional data.frame of fitted data e.g., `fit$data`

- data_labels:

  For each row of `data`, listing the corresponding name from
  `tree$tip.label`. Default pulls `data_labels` from `rownames(data)`

- distance_threshold:

  distance from taxa to use when plotting data

- tree:

  phylogenetic structure, using class
  [`as.phylo`](https://rdrr.io/pkg/ape/man/as.phylo.html)

- ...:

  Not used
