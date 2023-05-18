

context("Testing cross platform and R version compatibility")

# Eastern Bering Sea pollcok
test_that("phylosem example is working ", {
  #skip_on_ci()
  library(phylopath)
  library(phylobase)
  library(phylosignal)

  # Run phylosem
  model = "
    DD -> RS, p1
    BM -> LS, p2
    BM -> NL, p3
    NL -> DD, p4
  "
  psem = phylosem( sem = model,
          data = rhino[,c("BM","NL","DD","RS","LS")],
          tree = rhino_tree )
  # Check objective function
  expect_equal( as.numeric(psem$opt$obj), 1087.686, tolerance=1e-2 )

  # Convert and plot using phylopath
  coef_plot( as(psem,"fitted_DAG") )
  plot( as(psem,"fitted_DAG") )

  # Convet and plot using sem
  mysem = as(psem,"sem")
  sem::pathDiagram( model = mysem,
                  style = "traditional",
                  edge.labels = "values" )
  myplot = semPlot::semPlotModel( as(psem,"sem") )
  semPlot::semPaths( myplot,
                   nodeLabels = myplot@Vars$name )
  effects( as(psem,"sem") )

  # Convert and plot using phylobase / phylosignal
  plot( as(psem,"phylo4d") )
  barplot( as(psem,"phylo4d") )
  dotplot( as(psem,"phylo4d") )
  gridplot( as(psem,"phylo4d") )
})

