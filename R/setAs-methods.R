

#' Convert phylosem to other package formats
#'
#' Convert output from package phylosem to package-phylopath or package-sem
#'
#' @name setAs
#' @docType methods
#' @section Usage: \code{as(object, "fitted_DAG")}
#' @section Usage: \code{as(object, "sem")}
#' @seealso generic \code{\link[methods]{as}}, \code{\link[phylopath]{est_DAG}} from the
#' \code{phylopath} package, and \code{\link[sem]{sem}} from the
#' \code{sem} package.
#' @keywords methods
#' @aliases as as-method as,fitted_DAG,fitted_DAG-method,sem,sem-method
#'
#' @examples
#' \dontrun{
#' model = "
#'   DD -> RS, p1
#'   BM -> LS, p2
#'   BM -> NL, p3
#'   NL -> DD, p4
#' "
#' psem = phylosem( sem = model,
#'           data = rhino[,c("BM","NL","DD","RS","LS")],
#'           tree = rhino_tree )
#'
#' # Convert and plot using phylopath
#' library(phylopath)
#' coef_plot( as(psem,"fitted_DAG") )
#' plot( as(psem,"fitted_DAG") )
#'
#' # Convet and plot using sem
#' mysem = as(psem,"sem")
#' sem::pathDiagram( model = mysem,
#'                   style = "traditional",
#'                   edge.labels = "values" )
#' myplot = semPlot::semPlotModel( as(psem,"sem") )
#' semPlot::semPaths( myplot,
#'                    nodeLabels = Plot@Vars$name )
#' effects( as(psem,"sem") )
#' }
setAs("phylosem", "fitted_DAG", function(from, to) {

  # Base upon: https://github.com/fmichonneau/phylobase/blob/master/R/setAs-methods.R#L165-L228
  # extract and name identical to output from fitted_DAG
  out = list(
    coef = t(from$report$Rho_jj),
    se = t(as.list(from$opt$SD,what="Std. Error",report=TRUE)$Rho_jj)
  )
  dimnames(out$coef) = dimnames(out$se) = list( colnames(from$data), colnames(from$data) )

  # pass out
  class(out) = "fitted_DAG"
  return(out)

})

setAs("phylosem", "sem", function(from, to) {

  Sprime = from$report$V_jj
    rownames(Sprime) = colnames(Sprime) = colnames(from$data)
  out = sem::sem( from$SEM_model,
             S = Sprime,
             N = nrow(from$data) )   # data=data.frame(Database$Y_ij),

  # pass out
  return(out)

})


