# Organize stock-recruit time series

Provide observed stock biomass and recruitment time-series data, as well
as unfished spawning potential ratio and adult natural mortality rate,
to allow estimating maximum annual spawners per spawner (MASPS), which
can then be converted to steepness.

## Usage

``` r
beverton_holt(
  recruits,
  spawners,
  unfished_spawners_pre_recruit,
  natural_mortality,
  species
)
```

## Arguments

- recruits:

  Recruitment values for multiple stocks in a single named vector, with
  `names(recruits)` identifying the stock

- spawners:

  Spawning size (biomass or numbers) for multiple stocks in a named
  vector, where `spawners[1]` corresponds to the spawners that results
  in `recruits[1]`

- unfished_spawners_pre_recruit:

  SPR0 with units matching `recruits` and `spawners` for each stock,
  where `names(unfished_spawners_pre_recruit)` identifies the stock

- natural_mortality:

  annual adult natural mortality rate for each stock

- species:

  taxon label for each stock, matching `tree$tip.label` in the
  phylogenetic tree that will be used for analysis
