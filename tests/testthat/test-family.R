
context("Testing cross platform and R version compatibility")

test_that("phylosem Poisson-PCGLMM is working ", {
  skip_on_cran()
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

  #
  pglmm = phyr::pglmm_compare(
    y ~ 1 + x,
    family = "poisson",
    data = Data,
    phy = tree )

  #
  pgsem = phylosem(
    sem = "x -> y, p",
    data = Data,
    family = list( x = fixed(), y = poisson("log") ),
    tree = tree,
    control = phylosem_control(quiet = TRUE)
  )

  #
  expect_equal(
    as.numeric(subset(summary(pgsem)$coefficients, Path == "x -> y" )['Estimate']),
    as.numeric(pglmm$B['x',]),
    tol = 1e-2
  )
})

test_that("phylosem Binomial-PCGLMM is working ", {
  skip_on_cran()
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
  y_binom = rbinom( n=Ntree, size=1, prob=plogis(y_normal) )

  # Construct, re-order, and reduce data
  Data = data.frame(x=x,y=y_binom)

  #
  pglmm = phyr::pglmm_compare(
    y ~ 1 + x,
    family = "binomial",
    data = Data,
    phy = tree )

  #
  pgsem = phylosem( sem = "x -> y, p",
            data = Data,
            family = list( x = fixed(), y = binomial("logit") ),
            tree = tree,
            control = phylosem_control(quiet = TRUE) )

  pgsem2 = phylosem( sem = "x -> y, p",
            data = Data,
            family = list( x = fixed(), y = categorical("binary") ),
            tree = tree,
            control = phylosem_control(quiet = TRUE) )


  #
  expect_equal(
    as.numeric(subset(summary(pgsem)$coefficients, Path == "x -> y" )['Estimate']),
    as.numeric(pglmm$B['x',]),
    tol = 5e-2
  )
  expect_equal(
    as.numeric(subset(summary(pgsem)$coefficients, Path == "x -> y" )['Estimate']),
    as.numeric(subset(summary(pgsem2)$coefficients, Path == "x -> y" )['Estimate']),
    tol = 5e-3
  )
})

test_that("phylosem categorical-PCGLMM is working ", {
  skip_on_cran()
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
  y_binom = rbinom( n=Ntree, size=1, prob=plogis(y_normal) )

  # Construct, re-order, and reduce data
  Data = data.frame(x=x,y=y_binom)

  #
  pglmm = phyr::pglmm_compare(
    y ~ 1 + x,
    family = "binomial",
    data = Data,
    phy = tree )

  #
              #
  expect_equal(
    as.numeric(subset(summary(pgsem)$coefficients, Path == "x -> y" )['Estimate']),
    as.numeric(pglmm$B['x',]),
    tol = 5e-2
  )
})
