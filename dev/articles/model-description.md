# phylosem model description

## Phylogenetic structural equation models

`phylosem` is an R package for fitting phylogenetic structural equation
models (PSEMs). The package generalizes features in existing R packages:

- `sem` for structural equation models (SEMs);
- `phylosem` for comparison among alternative path models;
- `phylolm` for fitting large linear models that arise as when
  specifying a SEM with one endogenous variable and multiple exogenous
  and independent variables;
- `Rphylopars` for interpolating missing values when specifying a SEM
  with an unstructured (full rank) covariance among variables;

In model configurations that can be fitted by both `phylosem` and these
other packages, we have confirmed that results are nearly identical or
otherwise identified reasons that results differ.

`phylosem` involves a simple user-interface that specifies the SEM using
notation from package `sem` and the phylogenetic tree using package
`ape`. It allows uers to specify common models for the covariance
including:

- Brownian motion (BM);
- Ornstein-Uhlenbeck (OU);
- Pagel’s lambda;
- Pagel’s kappa;

Output can be coerced to standard formats so that `phylosem` can use
plotting and summary functions form other packages. Available output
formats include:

- `sem`, for plotting the estimated SEM and summarizing direct and
  indirect effects;
- `phylopath`, for plotting and model comparison;
- `phylo4d` in R-package `phylobase` for plotting estimated traits;

Package *phylosem* (Thorson and Bijl 2023) involves specifying a
phylogenetic structural equation model (PSEM). This PSEM be viewed
either:

1.  *Weak interpretation*: as an expressive interface to parameterize
    the correlation among life-history traits, using as many or few
    parameters as might be appropriate; or
2.  *Strong interpretation*: as a structural causal model, allowing
    predictions about the consequence of counterfactual changes to the
    system.

We introduce PSEM first from the perspective of a software user (i.e.,
the interface) and then from the perspective of a statistician (i.e.,
the equations and their interpretation).

#### Viewpoint 1: Software interface

To specify a PSEM, the user uses *arrow notation*derived from package
`sem` (Fox 2006). For example, to specify that evolutionary chnages in
trait $`x`$ causes a change in trait $`y`$, we specify:

``` r
x -> y, slope
x <-> x, sd_x
y <-> y, sd_y
```

This then estimates a single parameter representing the impact of $`x`$
on $`y`$ (specified with one-headed arrows), as well as the square root
(Cholesky decomposition) of the exogenous covariance of of model
variables (specified with two-headed arrows). See
[`?phylosem`](https://james-thorson-noaa.github.io/phylosem/dev/reference/phylosem.md)
Details section for more details about syntax.

PSEM interactions can be as complicated or simple as desired, and can
include:

1.  Latent variables and loops (i.e., they are not restricted to
    directed acyclic graphs);
2.  Values that are fixed a priori, where the `parameter_name` is
    provided as `NA` and the starting value that follows is the fixed
    value;
3.  Values that are mirrored among path coefficients, where the same
    `parameter_name` is provided for multiple rows of the text file.

The user also specifies a distribution for measurement errors for each
variable using arguement `family`, and a dated phylogenetic tree `tree`
that is used to represent evolutionary correlations.

#### Viewpoint 2: Mathematical details

The PSEM defines a generalized linear mixed model (GLMM) for a
$`V \times J`$ matrix $`\mathbf{Y}`$, where $`y_{vj}`$ is the
measurement for vertex $`v`$ for variable $`j`$, where vertex $`v`$ can
be either an ancestral node or a tip (i.e., PSEM estimates ancestral
traits jointly with tips). This measurement matrix can include missing
values $`y_{tj} = \mathrm{NA}`$, and it will estimate a $`V \times J`$
matrix of latent states $`\mathbf{V}`$ for all modeled vertices and
variables, where $`\mathbf{X}_v`$ is the vector of traits for vertex
$`v`$ and $`\text{vec}(\mathbf{X})`$ is a $`VJ`$ length vector
constructed by stacking columns. PSEM also estimates a $`V \times J`$
matrix of process errors $`\mathbf{E}`$ where $`\mathbf{\epsilon}_v`$ is
the vector of process errors for vertex $`v`$, and
$`\text{vec}(\mathbf{E})`$ is the $`TJ`$ vector of stacked columns.

PSEM can be written in structural vector-autoregressive (SVAR) notation
as a lag-1 SVAR model :

``` math
\mathbf{x}_{c_e} = \underbrace{\mathbf{P} \mathbf{x}_{c_e}}_{\text{Relationships among traits}} + \underbrace{\rho \mathbf{x}_{p_e}}_{\text{Relationship among taxa}} + \mathbf{\epsilon}_{c_e}
```
where this expression is evaluated for each edge $`e`$ of the
phylogenetic tree connecting child $`c_e`$ with parent $`p_e`$.
$`\mathbf{P}`$ is the $`V \times V`$ matrix of relationships among
traits, and $`\mathbf{\epsilon}_v`$ is exogenous covariation for vertex
$`v`$:

``` math
\mathbf{\epsilon}_t \sim \text{MVN}(\mathbf{0,GG}^T)
```
and $`\mathbf{G}`$ is the square-root of exogenous covariance
$`\mathbf{GG}^T`$ that occurs in each time. $`\mathbf{G}`$ is then
estimated by PSEM, and is identifiable without constraints when
specified as a lower-triangle matrix (representing the Cholesky
decomposition of the covariance of exogenous process errors occurring in
each time).

This SVAR notation can also be rewritten as a joint simultaneous
equation model (SEM), although we do not elaborate here
((**thorson_graphical_2026?**))

### Measurement errors

PSEM includes multiple distribution for measurement errors. For example,
if the user specifies `family[[j] = fixed()` then:

``` math
y_{vj} = x_{vj} + d_{j}
```
for all vertices, where $`d_j`$ is an estimated mean parameter that is
only estimated when specifying an Ornstein-Uhlenbeck process (and
otherwise fixed at zero). Alternatively, if the user specifies
`family[[j]] = gaussian()` then:
``` math
y_{vj} \sim \mathrm{Normal}( x_{vj} + d_{j}, {\sigma_j}^2)
```
and $`{\sigma_j}^2`$ is then included as an estimated parameter. When
estimating missing values or masurement errors in $`\mathbf{Y}`$,
*phylosem* must then marginalize across the latent value of states
$`\mathbf{X}`$. It does this using the Laplace approximation (Skaug and
Fournier 2006), as implemented using the R-package TMB (Kristensen et
al. 2016). Computations involving sparse matrices are efficient using
the Matrix package (Bates et al. 2023) to interface with the Eigen
library (Guennebaud et al. 2010).

## Works cited

Bates, Douglas, Martin Maechler, and Mikael Jagan. 2023. *Matrix: Sparse
and Dense Matrix Classes and Methods*.
<https://CRAN.R-project.org/package=Matrix>.

Fox, John. 2006. “Structural Equation Modeling with the Sem Package in
R.” *Structural Equation Modeling-a Multidisciplinary Journal* 13:
465–86.

Guennebaud, G., B. Jacob, et al. 2010. *Eigen V3*.
<https://libeigen.gitlab.io/>.

Kristensen, Kasper, Anders Nielsen, Casper W. Berg, Hans Skaug, and
Bradley M. Bell. 2016. “TMB: Automatic Differentiation and Laplace
Approximation.” *Journal of Statistical Software* 70 (5): 1–21.
<https://doi.org/10.18637/jss.v070.i05>.

Skaug, Hans, and Dave Fournier. 2006. “Automatic Approximation of the
Marginal Likelihood in Non-Gaussian Hierarchical Models.” *Computational
Statistics & Data Analysis* 51 (2): 699–709.

Thorson, James T., and Wouter van der Bijl. 2023. “Phylosem: A Fast and
Simple R Package for Phylogenetic Inference and Trait Imputation Using
Phylogenetic Structural Equation Models.” *Journal of Evolutionary
Biology* 36 (10): 1357–64. <https://doi.org/10.1111/jeb.14234>.
