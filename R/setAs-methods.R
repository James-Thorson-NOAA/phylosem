

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
#'
#' @aliases as as-method as,fitted_DAG,fitted_DAG-method,sem,sem-method
setAs("phylosem", "fitted_DAG", function(from, to) {

  # Base upon: https://github.com/fmichonneau/phylobase/blob/master/R/setAs-methods.R#L165-L228
  # extract and name identical to output from fitted_DAG
  out = list(
    coef = t(from$report$Rho_jj),
    se = t(as.list(from$opt$SD, what="Std. Error", report=TRUE)$Rho_jj)
  )
  dimnames(out$coef) = dimnames(out$se) = list( colnames(from$data), colnames(from$data) )

  # pass out
  class(out) = "fitted_DAG"
  return(out)

})

setAs("phylosem", "sem", function(from, to) {

  Sprime = from$report$V_jj
    rownames(Sprime) = colnames(Sprime) = colnames(from$data)
  out = sem( from$SEM_model,
             S = Sprime,
             N = nrow(from$data) )   # data=data.frame(Database$Y_ij),

  # pass out
  return(out)

})

setAs("phylosem", "phylo4d", function(from, to) {

  #
  traits = from$report$x_vj
  colnames(traits) = colnames(from$data)
  out = phylo4d( x=from$tree, all.data=traits )

  # pass out
  return(out)

})

