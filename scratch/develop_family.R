
setwd( R'(C:\Users\jtuth\Documents\GitHub\phylosem\src)' )
TMB::compile( "phylosem.cpp", framework = "TMBad" )

pak::local_install(R'(C:\Users\jtuth\Documents\GitHub\phylosem)')

library(phylosem)

# Settings
Ntree = 100
sd_x = 0.3
sd_y = 0.3
b0_x = 1
b0_y = 0
b_xy = 1

# Simulate tree
set.seed(1)
tree = ape::rtree(n=Ntree)

# Simulate data
x = b0_x + sd_x * phylolm::rTrait(n = 1, phy=tree)
ybar = b0_y + b_xy*x
y_normal = ybar + sd_y * phylolm::rTrait(n = 1, phy=tree)
y_pois = rpois( n=Ntree, lambda=exp(y_normal) )

# Construct, re-order, and reduce data
Data = data.frame(x=x,y=y_pois)

# Compare using phylolm::phyloglm
pglm = phylolm::phyloglm( y ~ 1 + x, data=Data, phy=tree, method="poisson_GEE" )
knitr::kable(summary(pglm)$coefficients, digits=3)

#
pglmm = phyr::pglmm_compare(
  y ~ 1 + x,
  family = "poisson",
  data = Data,
  phy = tree )
knitr::kable(summary(pglmm), digits=3)

#
pgsem = phylosem( sem = "x -> y, p",
          data = Data,
          #family = c("fixed","poisson"),
          family = list( x = fixed(), y = poisson("log") ),
          tree = tree,
          control = phylosem_control(quiet = TRUE) )

#
expect_equal(
  as.numeric(subset(summary(pgsem)$coefficients, Path == "x -> y" )['Estimate']),
  as.numeric(pglmm$B['x',]),
  tol = 1e-2
)
