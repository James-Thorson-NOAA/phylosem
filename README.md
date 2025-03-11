[![](https://www.r-pkg.org/badges/version/phylosem)](https://cran.r-project.org/package=phylosem)
[![](https://cranlogs.r-pkg.org/badges/phylosem)](https://cran.r-project.org/package=phylosem)
[![](https://cranlogs.r-pkg.org/badges/grand-total/phylosem)](https://cran.r-project.org/package=phylosem)
[![Codecov test coverage](https://codecov.io/gh/James-Thorson-NOAA/phylosem/branch/remove-fit_tmb/graph/badge.svg)](https://app.codecov.io/gh/James-Thorson-NOAA/phylosem/tree/remove-fit_tmb)
[![Documentation](https://img.shields.io/badge/documentation-phylosem-orange.svg?colorB=E91E63)](https://james-thorson-noaa.github.io/phylosem/)

## Phylogenetic structural equation models
Package _phylosem_ combines features from structural equation models (SEM), phylogenetic comparative methods (PCM), and generalized linear mixed models (GLMM).  By doing so, it incorporates a broad feature-set:

* Comparing multiple evolutionary models similar to _phylopath_
* Estimating trade-offs among multiple traits similar to _phylolm_
* Predicting missing trait values (and associated standard errors) similar to _Rphylopars_
* Estimating tradeoffs including recursive (cyclic) dependencies similar to package _sem_
* Applying ordination to multiple traits, similar to phylogenetic factor analysis in package _FishLife_

_phylosem_ is specifically intended as a minimal implementation, and uses standard packages for input/output formatting:

* Input: phylogenetic relatedness defined using class _phylo_ in package _ape_
* Input: structural trade-offs specified using syntax defined by package _sem_
* Output: visualizing trade-offs using _semPlot_, _diagrammeR_, and _ggraph_
* Output: assembling trait predictions and standard errors using _phylobase_
* Output: plotting trait predictions using _phylosignal_

Please see package vignettes for more details regarding syntax and features.

# NOAA Enterprise GitHub disclaimer
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

