#' Convert phylosem to sem
#'
#' Convert output from package phylosem to package sem
#'
#' @param x Output from \code{\link{phylosem}}
#'
#' @examples
#' \dontrun{
#' library(phylopath)
#' model = "
#'   DD -> RS, p1
#'   BM -> LS, p2
#'   BM -> NL, p3
#'   NL -> DD, p4
#' "
#' psem = phylosem( sem = model,
#'           data = rhino[,c("BM","NL","DD","RS","LS")],
#'           tree = rhino_tree )
#' sem::pathDiagram( model = phylosem2sem(psem),
#'                   style = "traditional",
#'                   edge.labels = "values" )
#' Plot = semPlot::semPlotModel( phylosem2sem(psem) )
#' semPlot::semPaths( Plot,
#'                    nodeLabels = Plot@Vars$name )
#' effects( phylosem2sem(psem) )
#' }
#'
#' @export
phylosem2sem <-
function( x ){

  Sprime = psem$report$V_jj
    rownames(Sprime) = colnames(Sprime) = colnames(psem$data)
  out = sem::sem( psem$SEM_model,
             S = Sprime,
             N = nrow(psem$data) )   # data=data.frame(Database$Y_ij),

  # pass out
  return(out)
}

