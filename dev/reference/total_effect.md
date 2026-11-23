# Calculate total effects

Calculate a data frame of total effects representing the estimated
effect of a change in any variable on every other variable.

## Usage

``` r
total_effect(object)
```

## Arguments

- object:

  Output from
  [`phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/phylosem.md)

## Value

A data frame listing the time-lag (lag), variable that is undergoing
some exogenous change (from), and the variable being impacted (to),
along with the total effect (total_effect) including direct and indirect
pathways, and the partial "direct" effect (direct_effect)

## Details

Total effects are taken from the Leontief matrix
\\\mathbf{(I-P)^{-1}}\\, where \\\mathbf{P}\\ is the path matrix across
variables. \\\mathbf{(I-P)}^{-1} \mathbf{\delta} \\ calculates the
effect of a perturbation represented by vector \\\mathbf{\delta}\\
