# Convert phylosem to phylo4d

Convert output from package phylosem to phylo4d object from package
phylobase

## Usage

``` r
as_phylo4d(object, what = c("Estimate", "Std. Error"))
```

## Arguments

- object:

  Output from
  [`phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/phylosem.md)

- what:

  Select what to convert (Estimate / Std. Error).

## Value

phylosem output to converted format supplied by
[`phylo4d`](https://rdrr.io/pkg/phylobase/man/phylo4d-methods.html)

## Details

This package is intended to for use in using plots assocaited with
package sem, e.g., using package plotSEM
[`semPlot::semPlotModel`](https://rdrr.io/pkg/semPlot/man/semPlotModel.html)
