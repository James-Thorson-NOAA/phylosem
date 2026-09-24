# phylosem

Package *phylosem* combines features from structural equation models
(SEM), phylogenetic comparative methods (PCM), and generalized linear
mixed models (GLMM). By doing so, it incorporates a broad feature-set:

- Comparing multiple evolutionary models similar to *phylopath*
- Estimating trade-offs among multiple traits similar to *phylolm*
- Predicting missing trait values (and associated standard errors)
  similar to *Rphylopars*
- Estimating tradeoffs including recursive (cyclic) dependencies similar
  to package *sem*
- Applying ordination to multiple traits, similar to phylogenetic factor
  analysis in package *FishLife*

*phylosem* is specifically intended as a minimal implementation, and
uses standard packages for input/output formatting:

- Input: phylogenetic relatedness defined using class *phylo* in package
  *ape*
- Input: structural trade-offs specified using syntax defined by package
  *sem*
- Output: visualizing trade-offs using *semPlot*, *diagrammeR*, and
  *ggraph*
- Output: assembling trait predictions and standard errors using
  *phylobase*
- Output: plotting trait predictions using *phylosignal*

Please see package vignettes for more details regarding syntax and
features.

[![](https://www.r-pkg.org/badges/version/phylosem)](https://cran.r-project.org/package=phylosem)
[![](https://cranlogs.r-pkg.org/badges/phylosem)](https://cran.r-project.org/package=phylosem)
[![](https://cranlogs.r-pkg.org/badges/grand-total/phylosem)](https://cran.r-project.org/package=phylosem)
[![Codecov test
coverage](https://codecov.io/gh/James-Thorson-NOAA/phylosem/branch/remove-fit_tmb/graph/badge.svg)](https://app.codecov.io/gh/James-Thorson-NOAA/phylosem/tree/remove-fit_tmb)
[![Documentation](https://img.shields.io/badge/documentation-phylosem-orange.svg?colorB=E91E63)](https://james-thorson-noaa.github.io/phylosem/)
