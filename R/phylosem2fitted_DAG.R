#' Convert phylosem to fitted_DAG
#'
#' Convert output from package phylosem to package phylopath
#'
#' @param x Output from \code{\link{phylosem}}
#'
#' @export
phylosem2fitted_DAG <-
function( x ){

  # extract and name identical to output from fitted_DAG
  out = list(
    coef = t(psem$report$Rho_jj),
    se = t(as.list(psem$opt$SD,what="Std. Error",report=TRUE)$Rho_jj)
  )
  dimnames(out$coef) = dimnames(out$coef) = list( colnames(psem$data), colnames(psem$data) )

  # pass out
  class(out) = "fitted_DAG"
  return(out)
}



