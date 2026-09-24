# NOAA Enterprise GitHub disclaimer

## Phylogenetic structural equation models

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

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.
