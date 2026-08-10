
context("Testing ´random slopes and parameter mapping")

# Eastern Bering Sea pollcok
test_that("phylosem Poisson-PCGLMM is working ", {
  set.seed(123)

  library(phytools)
  library(ape)

  n_tips <- 100
  tree <- pbtree(n = n_tips, scale = 1)  # phytools; scale=1 sets total tree depth to 1

  # y _tilda_ 1 + w * x
  x <- rTraitCont( tree, sigma = 0.5, alpha = 1, model = "BM")
  w <- 1 + rTraitCont( tree, sigma = 0.5, alpha = 1, model = "BM")
  y = 1 + w * x + rTraitCont( tree, sigma = 0.5, alpha = 1, model = "BM")
  z = rTraitCont( tree, sigma = 0.5, alpha = 1, model = "BM")

  data = data.frame( x = x, y = y, w = NA, z = z)

  # Estimate random slope, but with extra mapping as well
  sem = "
    x -> y, w
    x <-> x, sigma1, 0.5
    w <-> w, sigma2, 0.1
    y <-> y, sigma1, 0.5
    z <-> z, NA, 1
  "

  fit = phylosem(
    data = data,
    tree = tree,
    sem = sem,
    estimate_xbar = c("x", "y", "w"),
    estimate_ou = FALSE,
    control = phylosem_control(
      trace = 1
    )
  )

  expect_equal( abs(as.numeric(fit$opt$par[1:2])), c(0.495, 1.238), tolerance = 0.01 )
})
