#' Compare phylogenetic structural equation models
#'
#' Fits several phylogenetic structural equation model for further comparison
#'
#' @inheritParams phylosem
#'
#' @param sem_set A named list of structural equation model specifications,
#'        where each element will be passed as argument \code{sem} to
#'        \code{\link{phylosem}}
#'
#' @export
compare_phylosem <-
function( sem_set,
          tree,
          data,
          family = rep("fixed", ncol(data)),
          estimate_theta = FALSE,
          estimate_lambda = FALSE,
          estimate_kappa = FALSE,
          ... ){

  # Initialize object of results
  out = vector( "list", length=length(sem_set) )
  names(out) = names(sem_set)

  # Loop through models
  for( modelI in seq_along(sem_set) ){
    fit = tryCatch( phylosem( sem = sem_set[modelI],
                              tree = tree,
                              data = data,
                              family = family,
                              estimate_theta = estimate_theta,
                              estimate_lambda = estimate_lambda,
                              estimate_kappa = estimate_kappa,
                              ... ),
                    error = function(e) e )
    if( "phylosem" %in% class(fit) ){
      out[[modelI]] = fit
    }
  }

  # pass out
  return(out)
}

