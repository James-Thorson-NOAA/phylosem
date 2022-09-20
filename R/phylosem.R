#' Fit phylogenetic structural equation model
#'
#' Fits a phylogenetic structural equation model
#'
#' @param sem structural equation model structure, passed to either \code{\link[sem]{specifyModel}}
#'        or \code{\link[sem]{specifyEquations}} and then parsed to control
#'        the set of path coefficients and variance-covariance parameters
#' @param tree phylogenetic structure, using class \code{\link[ape]{as.phylo}}
#' @param data data-frame providing variables being modeled.  Missing values are inputted
#'        as NA.  If an SEM includes a latent variable (i.e., variable with no available measurements)
#'        then it still must be inputted as a column of \code{data} with entirely NA values.
#' @param family Character-vector listing the distribution used for each column of \code{data}, where
#'        each element must be \code{fixed}, \code{normal}, \code{binomial}, or \code{poisson}.
#' @param estimate_theta Boolean indicating whether to estimate an autoregressive (Ornstein-Uhlenbeck)
#'        process
#' @param estimate_lambda Boolean indicating whether to estimate additional branch lengths for
#'        phylogenetic tips (a.k.a. the Pagel-lambda term)
#' @param estimate_kappa Boolean indicating whether to estimate a nonlinear scaling of branch
#'        lengths (a.k.a. the Pagel-kappa term)
#'
#' @useDynLib phylosem
#' @export
phylosem <-
function( sem,
          tree,
          data,
          family = rep("fixed", ncol(data)),
          estimate_theta = FALSE,
          estimate_lambda = FALSE,
          estimate_kappa = FALSE,
          ... ){

  # Function that converts SEM model to a RAM, see `?sem` for more context
  build_ram = function( model, vars ){
    vars = sapply( vars, FUN=function(char){gsub("-", "", gsub(" ", "", char))} )
    n.paths = nrow(model)
    par.names = model[, 2]
    startvalues = model[,3]

    # EXCERPT FROM `getAnywhere("sem.semmod")`
    heads = from = to = rep(0, n.paths)
    for (p in 1:n.paths) {
      path = sem:::parse.path(model[p, 1])
      heads[p] = abs(path$direction)
      to[p] = path$second
      from[p] = path$first
      if (path$direction == -1) {
        to[p] = path$first
        from[p] = path$second
      }
    }
    missing_vars = setdiff( c(from,to), vars )
    if( length(missing_vars) > 0 ) stop( "Check `build_ram`:", paste0(missing_vars,sep=", ") )

    ram = data.frame(matrix(0, nrow=p, ncol=5))
    pars = na.omit(unique(par.names))
    ram[, 1] = heads
    ram[, 2] = apply(outer(vars, to, "=="), 2, which)
    ram[, 3] = apply(outer(vars, from, "=="), 2, which)
    par.nos = apply(outer(pars, par.names, "=="), 2, which)
    if(length(par.nos) > 0){
      ram[, 4] = unlist(lapply(par.nos, function(x) if (length(x)==0){0}else{x}))
    }
    ram[, 5] = startvalues
    colnames(ram) = c("heads", "to", "from", "parameter", "start")
    return(ram)
  }

  # Errors / warnings;  not using `assertthat::assert_that` to avoid dependency
  #if( estimate_theta==FALSE & fixed_root==TRUE ){
  #  stop("`fixed_root=TRUE` is only applicable when `estimate_theta=TRUE`")
  #}
  #if( measurement_error==TRUE & estimate_lambda==TRUE ){
  #  stop("measurement errors using `measurement_error=TRUE` are confounded with `estimate_lambda=TRUE`")
  #}
  if( isFALSE("phylo" %in% class(tree)) ){
    stop("Check `tree` input")
  }
  familycode_j = sapply( tolower(family), FUN=switch, "fixed"=0, "normal"=1, "norm"=1, "binomial"=2, "binom"=2, "poisson"=3, "pois"=3, NA )
  if( any(is.na(familycode_j)) ) stop("Check `family`")

  #
  SEM_model = tryCatch(
    sem::specifyModel( text=sem, exog.variances=TRUE, endog.variances=TRUE, covs=colnames(data) ),
    error = function(e) e
  )
  if( isFALSE("semmod" %in% class(SEM_model)) ){
    SEM_model = tryCatch(
      sem::specifyEquations( text=sem, exog.variances=TRUE, endog.variances=TRUE, covs=colnames(data) ),
      error = function(e) e
    )
  }
  if( isFALSE("semmod" %in% class(SEM_model)) ){
    stop("Must supply either input for `sem::specifyModel` or `sem::specifyEquations`")
  }
  RAM = build_ram( SEM_model, colnames(data) )

  #
  n_tip = ape::Ntip(tree)
  vroot = n_tip + 1
  edge_ez = tree$edge  # parent, child
  length_e = tree$edge.length
  height_v = ape::node.depth.edgelength(tree)
  if(vroot %in% edge_ez[,2]) stop("Check for problems")

  # Build data
  data_list = list( "n_tip" = n_tip,
                    "edge_ez" = edge_ez - 1,
                    "length_e" = length_e,
                    "RAM" = as.matrix(RAM[,1:4]),
                    "RAMstart" = as.numeric(RAM[,5]),
                    "estimate_theta" = estimate_theta,
                    "estimate_lambda" = estimate_lambda,
                    "estimate_kappa" = estimate_kappa,
                    "height_v" = height_v,
                    "y_ij" = as.matrix(data),
                    "familycode_j" = familycode_j )

  # Build parameters
  rmatrix = function(nrow, ncol) matrix(rnorm(nrow*ncol),nrow=nrow,ncol=ncol)
  parameters_list = list( "beta_z" = rep(0.1, max(RAM[,4])),
                          "lnsigma_j" = rep(0,ncol(data)),
                          "lntheta" = log(1),
                          "logitlambda" = plogis(2),
                          "lnkappa" = log(1),
                          "x_vj" = 0.1 * rmatrix( nrow=nrow(edge_ez)+1, ncol=ncol(data) ),
                          "xbar_j" = rep(0,ncol(data)) )
  # Build map
  map_list = list()
  # Start off map_list$x_vj, which has multiple constraints
  map_list$x_vj = array( 1:prod(dim(parameters_list$x_vj)), dim=dim(parameters_list$x_vj) )
  map_list$x_vj[data_list$n_tip+1,] = ifelse(colSums(!is.na(data))==0, NA, 1:ncol(data))

  # Build random
  random = c("x_vj")

  # Settings
  map_list$lnsigma_j = 1:length(parameters_list$lnsigma_j)
  for( j in 1:ncol(data)){
    # Turn off SD of measurement error
    if( familycode_j[j] %in% c(0,2,3) ){
      map_list$lnsigma_j[j] = NA
    }
    if( familycode_j[j] == 0 ){
      # Fix random-effects for tips at their observed values
      parameters_list$x_vj[1:n_tip,j] = ifelse( is.na(data_list$y_ij[,j]), parameters_list$x_vj[1:n_tip,j], data_list$y_ij[,j] )
      map_list$x_vj[1:n_tip,j] = ifelse( is.na(data_list$y_ij[,j]), map_list$x_vj[1:n_tip,j], NA )
    }
  }
  if( estimate_theta==FALSE ){
    map_list$lntheta = factor(NA)
  }
  if( estimate_theta==FALSE ){
    map_list$xbar_j = factor(rep( NA, length(parameters_list$xbar_j) ))
  }else{
    map_list$xbar_j = factor( ifelse(colSums(!is.na(data))==0, NA, 1:ncol(data)) )
  }
  if( estimate_lambda==FALSE ){
    map_list$logitlambda = factor(NA)
  }
  if( estimate_kappa==FALSE ){
    map_list$lnkappa = factor(NA)
  }

  # wrap up map_list$x_vj
  map_list$x_vj = factor(map_list$x_vj)
  map_list$lnsigma_j = factor(map_list$lnsigma_j)

  # Hardwire TMB using local path
  if( FALSE ){
    #dyn.unload( TMB::dynlib("phylosem") )          #
    setwd( "C:/Users/James.Thorson/Desktop/Git/phylosem/src/" )
    TMB::compile( "phylosem.cpp", flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign" )
    dyn.load( TMB::dynlib("phylosem") )          #
  }

  # Build TMB object
  obj = TMB::MakeADFun( data = data_list,
                        parameters = parameters_list,
                        map = map_list,
                        random = random,
                        DLL = "phylosem" )
  list_parameters(obj)
  obj$env$beSilent()       # if(!is.null(Random))
  results = list( data=data, SEM_model=SEM_model, obj=obj, call=match.call(),
                  map_list=map_list, parameters_list=parameters_list, data_list=data_list )
  #return(results)

  #
  results$opt = TMBhelper::fit_tmb( obj, ... )
  results$report = obj$report()
  results$parhat = obj$env$parList()
  class(results) = "phylosem"
  return( results )
}

#' Print parameter estimates and standard errors.
#'
#' @title Print parameter estimates
#' @param x Output from \code{\link{phylosem}}
#' @param ... Not used
#' @return NULL
#' @method print phylosem
#' @export
print.phylosem <- function(x, ...)
{
  cat("phylosem(.) result\n")
  if( "opt" %in% names(x) ){
    print( x$opt )
  }else{
    cat("`parameter_estimates` not available in `phylosem`\n")
  }
  invisible(x$opt)
}

#' Extract path coefficients.
#'
#' @title Extract path coefficients
#' @param x Output from \code{\link{phylosem}}
#' @param standardize Whether to standardize regression coefficients
#' @param ... Not used
#' @return NULL
#' @method coef phylosem
#' @export
coef.phylosem = function( x, standardized=FALSE ){
  beta_z = x$opt$par[names(x$opt$par)=="beta_z"]
  RAM = x$obj$env$data$RAM
  if(nrow(RAM) != nrow(x$SEM_model)) stop("Check assumptions")
  for( i in which(RAM[,1]==1) ){
    if( standardized==TRUE ){
      beta_z[i] = beta_z[i] * abs(beta_z[which( RAM[,'from']==RAM[i,'from'] & RAM[,'to']==RAM[i,'from'] )])
      beta_z[i] = beta_z[i] / abs(beta_z[which( RAM[,'from']==RAM[i,'to'] & RAM[,'to']==RAM[i,'to'] )])
    }
  }
  # Report variances
  for( i in which(RAM[,1]==2) ){
    if( RAM[i,'from'] == RAM[i,'to'] ){
      beta_z[i] = beta_z[i]^2
    }
  }
  SEM_params = beta_z[ifelse(RAM[,4]==0, NA, RAM[,4])]
  SEM_params = ifelse( is.na(SEM_params), as.numeric(x$SEM_model[,3]), SEM_params )
  return( data.frame(Path=x$SEM_model[,1], Parameter=x$SEM_model[,2], Estimate=SEM_params ) )
}