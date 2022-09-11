#' Fit phylogenetic structural equation model
#'
#' Fits a phylogenetic structural equation model
#'
#' @useDynLib phylosem
#' @export
phylosem <-
function( text,
          tree,
          data,
          measurement_error = FALSE,
          estimate_theta = TRUE,
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

  #
  SEM_model = sem::specifyModel( text=text, exog.variances=TRUE, endog.variances=TRUE, covs=colnames(data) )
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
                    "estimate_theta" = estimate_theta,
                    "estimate_lambda" = estimate_lambda,
                    "estimate_kappa" = estimate_kappa,
                    "height_v" = height_v,
                    "y_ij" = as.matrix(data) )

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

  # Build random
  random = c("x_vj")

  # Settings
  if( measurement_error==FALSE ){
    map_list$lnsigma_j = factor(rep( NA, length(parameters_list$lnsigma_j) ))
    parameters_list$lnsigma_j = rep( log(0.001), length(parameters_list$lnsigma_j) )
  }
  if( estimate_theta==FALSE ){
    map_list$lntheta = factor(NA)
    map_list$xbar_j = factor(rep( NA, length(parameters_list$xbar_j) ))
  }
  if( estimate_lambda==FALSE ){
    map_list$logitlambda = factor(NA)
  }
  if( estimate_kappa==FALSE ){
    map_list$lnkappa = factor(NA)
  }

  # Compile TMB
  if( TRUE ){
    #dyn.unload( TMB::dynlib("phylosem") )          #
    setwd( "C:/Users/James.Thorson/Desktop/Git/phylosem/src/" )
    TMB::compile( "phylosem.cpp", flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign" )
    dyn.load( TMB::dynlib("phylosem") )          #
  }

  # Build TMB object
  Obj = TMB::MakeADFun( data = data_list,
                        parameters = parameters_list,
                        map = map_list,
                        random = random,
                        DLL = "phylosem" )
  ThorsonUtilities::list_parameters(Obj)
  Obj$env$beSilent()       # if(!is.null(Random))

  #
  Opt = TMBhelper::fit_tmb( Obj, ... )
  Report = Obj$report()
  ParHat = Obj$env$parList()

  return( list(data=data, SEM_model=SEM_model, Opt=Opt, Report=Report, ParHat=ParHat) )
}
