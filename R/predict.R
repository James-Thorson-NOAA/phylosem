#' @title predict phylosem
#'
#' @description extract the predictive mean and covariance
#'        for nodes or tips based on names from tree,
#'        \code{x$tree$tip.label} and \code{x$tree$node.label}
#'
#' @param object Output from \code{\link{phylosem}}
#' @param taxon vector of taxon names (matched against tree) to predict traits
#' @param what What output to calculate
#'        (currently only the linear predictor is supported)
#' @param se.fit Whether to report predictive covariance, where
#'        \code{sqrt(diag(pred$cov))} is the standard errors
#' @param include_fixed Whether to include uncertainty in fixed effects
#'        when calculating predictive covariance
#' @param ... Not used
#'
#' @return
#' When \code{se.fit = FALSE}, it returns a matrix of predicted traits (rows) by taxon (column).
#' When \code{se.fit = TRUE}, it returns a list:
#' \describe{
#' \item{est}{matrix of predicted traits}
#' \item{cov}{array of predictive covariance among traits (1st and 2nd dims) by taxon (3rd dim)}
#' }
#'
#' @method predict phylosem
#' @export
predict.phylosem <-
function( object,
          taxon,
          what = c("p_g"),
          se.fit = TRUE,
          include_fixed = TRUE,
          ... ){

  what = match.arg(what)
  vi = match( taxon, c(object$tree$tip.label,object$tree$node.label) )
  ni = length(vi)
  nj = ncol(object$data)
  if(length(vi) != length(taxon)) stop("`taxon` not found in `tree$tip.label` or `tree$node.label`")

  # Return existing prediction
  out = t(object$report$x_vj[vi,,drop=FALSE])
  dimnames(out) = list( variable = colnames(object$data), taxon = taxon )

  #
  if( isTRUE(se.fit) ){
    if( isTRUE(include_fixed) ){
      if( is.null(object$sdrep$jointPrecision) ){
        stop("rerun with `getJointPrecision = TRUE`")
      }
      Q_nn = object$sdrep$jointPrecision
      L_nn = object$L_jointPrecision
    }else{
      stop("Option not added yet")
    }

    # Extract cells of joint covariance including NA for fixed values that are missing
    xmap_vj = matrix( as.numeric(object$tmb_inputs$map_list$x_vj), ncol = nj, byrow = FALSE )
    xmap_j = as.vector(xmap_vj[vi,])

    # Indexing
    p = na.omit(xmap_j)
    n = nrow(Q_nn)
    nk = length(p)

    # Do one forward solve
    E_nk = Matrix::sparseMatrix(i = p, j = seq_len(nk), x = 1, dims = c(n, nk))
    PE_nk = Matrix::solve(L_nn, E_nk, system = "P")
    W_nk  = Matrix::solve(L_nn, PE_nk, system = "L")
    V_kk  = as.matrix(Matrix::crossprod(W_nk))

    # Add zeros for fixed variables
    pfull = as.numeric(factor(xmap_j))
    Vfull = V_kk[pfull, pfull]
    Vfull[is.na(Vfull)] = 0

    # Make array
    V = array(
      sapply(seq_len(ni), function(i) {
        idx = seq( nj*(i - 1) + 1, i * nj )
        Vfull[idx, idx]
      }),
      dim = c(nj, nj, ni)
    )
    dimnames(V) = list( variable = colnames(object$data), variable = colnames(object$data), taxon = taxon )

    # Add to output
    out = list(
      est = out,
      cov = V
    )
  }

  class(out) = "predict_phylosem"
  return( out )
}


node_to_tip_dist <- function(phy, node) {
  if(is.character(node)){
    node = match( node, c(phy$tip.label, c(phy$node.label)) )
  }

  Nt   <- Ntip(phy)
  Nn   <- Nt + phy$Nnode
  anc  <- phy$edge[, 1]
  desc <- phy$edge[, 2]
  el   <- phy$edge.length

  # Put edges in an order where parents are processed before children
  phy <- ape::reorder(phy, "cladewise")
  anc <- phy$edge[, 1]; desc <- phy$edge[, 2]; el <- phy$edge.length

  parent  <- integer(Nn)
  parent[desc] <- anc
  elen_up <- numeric(Nn)
  elen_up[desc] <- el      # length of edge above each node

  # Mark the path from `node` up to the root
  on_path <- logical(Nn)
  path_dist <- numeric(Nn)                          # distance from `node` to each ancestor
  x <- node
  d <- 0
  on_path[x] <- TRUE
  while (parent[x] != 0L) {
    d <- d + elen_up[x]
    x <- parent[x]
    on_path[x] <- TRUE
    path_dist[x] <- d
  }

  # Distance from `node` to every node:
  #  - ancestors: path_dist
  #  - others: distance to nearest on-path ancestor (via preorder sweep) + path_dist of that ancestor
  dist <- numeric(Nn)
  dist[on_path] <- path_dist[on_path]
  for (i in seq_along(desc)) {                      # cladewise = parents before children
    ch <- desc[i]
    if (!on_path[ch]) dist[ch] <- dist[anc[i]] + el[i]
  }
  dist[1:Nt]
}

#' @title plot predictions from phylosem
#'
#' @description plot the predictive distribution for one or two variables
#'   extracted using \code{\link{predict.phylosem}}
#'
#' @inheritParams phylosem
#'
#' @param x Output from \code{\link{predict.phylosem}}
#' @param y Not used
#' @param variables which variables to plot (can only plot one or two at a time)
#' @param data optional data.frame of fitted data e.g., \code{fit$data}
#' @param distance_threshold distance from taxa to use when plotting data
#' @param ... Not used
#'
#' @method plot predict_phylosem
#' @export
plot.predict_phylosem <-
function( x,
          y,
          variables = 1,
          data,
          data_labels,
          distance_threshold = 0,
          tree,
          ... ){

  #
  if( all(is.numeric(variables)) ){
    variables = rownames(x$est)[variables]
  }
  assertCharacter(variables, any.missing = FALSE, min.len = 1, max.len = 2 )
  assertSubset(variables, rownames(x$est) )
  taxa = dimnames(x$cov)[[3]]

  #
  if( !missing(data) & !missing(data_labels) ){
    if( missing(tree) ){
      which_rows = which(data_labels %in% taxa)
    }else{
      which_rows = lapply(
        taxa,
        FUN = \(taxon){
          D = node_to_tip_dist( tree, taxon )
          near = tree$tip.label[which(D <= distance_threshold)]
          which(data_labels %in% near)
        }
      )
    }
    obs = data[ unlist(which_rows), variables ]
    obs_labels = rep( taxa, each = sapply(which_rows,length) )
  }else{
    obs = data.frame(Map(function(.) numeric(0), variables))
    obs_labels = vector(length=0)
  }


  # Avoid numerical singularity
  nj = dim(x$cov)[[1]]
  ni = dim(x$cov)[[3]]
  cov = x$cov + 1e-12 * diag(nj) %o% rep(1,ni)

  if( length(variables) == 1 ){
    from = x$est[variables,] - 3 * sqrt(x$cov[variables,variables,])
    to = x$est[variables,] + 3 * sqrt(x$cov[variables,variables,])
    xval = seq( min(c(from,obs[[1]]),na.rm=TRUE), max(c(to,obs[[1]]),na.rm=TRUE), length = 1000 )
    yval = sapply(
      seq_len(ni),
      function(i) dnorm(xval, x$est[variables,i], sqrt(x$cov[variables,variables,i]) )
    )

    # plot densities
    matplot(
      x = xval,
      y = yval,
      type = "l", lty = "solid",
      col = viridisLite::viridis(n = ni),
      ylab = "Distribution for mean",
      xlab = variables
    )
    legend( "topright", fill = viridisLite::viridis(n = ni), legend = taxa, bty = "n" )

    # Plot data
    points(
      x = obs, y = obs*0,
      col = viridisLite::viridis(n = ni)[match(obs_labels,taxa)]
    )
  }else if( length(variables) == 2 ){
    make_ellipse = function(mu, Sigma, n = 200) {   # p = 0.95,
      theta = seq(0, 2*pi, length.out = n)
      circle = rbind(cos(theta), sin(theta))  # * sqrt(qchisq(p, df = 2))
      t(mu + t(chol(Sigma)) %*% circle)   # n x 2 matrix of points
    }

    # set up the plot region to fit the largest ellipse
    ellipse_i = array(NA, dim = c(200,2,ni))
    for(i in seq_len(ni)){
      ellipse_i[,,i] = make_ellipse(x$est[variables,i], cov[variables,variables,i])
    }
    xylim = apply( ellipse_i, MARGIN = 2, FUN = range ) # rows: min, max.  cols: x, y

    # Plot stuff
    plot( 1, type = "n", xlim = xylim[,1], ylim = xylim[,2], xlab = variables[1], ylab = variables[2] )
    for(i in seq_len(ni)){
      lines(
        ellipse_i[,,i],
        col = viridisLite::viridis(n = ni)[i]
      )
      points(
        x$est[variables[1],i],
        x$est[variables[2],i],
        pch = 3, col = viridisLite::viridis(n = ni)[i]
      )
    }
    legend( "topright", fill = viridisLite::viridis(n = ni), legend = taxa, bty = "n" )

    # Plot data
    points(
      x = rep(xylim[1,1], nrow(obs)), y = obs[,2],
      col = viridisLite::viridis(n = ni)[match(obs_labels,taxa)]
    )
    points(
      x = obs[,1], y = rep(xylim[1,2], nrow(obs)),
      col = viridisLite::viridis(n = ni)[match(obs_labels,taxa)]
    )
  }
}


#' @title predict values for new tip
#'
#' @description Predict phylosem
#'
#' @inheritParams phytools::bind.tip
#'
#' @param x Output from \code{\link{phylosem}}
#' @param ... passed to \code{\link[phytools]{bind.tip}}
#'
#' @return predict new values
#'
#' @method predict phylosem
#' @export
#predict.phylosem <-
#function( x,
#          tip.label = "new_tip",
#          edge.length = NULL,
#          where = NULL,
#          ... ){
#
#  if( is.character(where) ){
#    where = which( c(x$tree$tip.label,x$tree$node.label) == where )
#    if(length(where)!=1) stop("`where` not found in `tree$tip.label` or `tree$node.label`")
#  }
#
#  # Return existing prediction OR rebuild
#  if( where <= Ntip(x$tree) ){
#    out = x$report$x_vj[where,]
#  }else{
#    # Add default edge.length
#    if(is.null(edge.length)){
#      if(is.ultrametric(x$tree)){
#        # node.depth, node.height, node.depth.edgelength
#        node_depth = node.depth.edgelength(x$tree)
#        edge_length = x$tree$edge.length
#        node_edgeout_length = edge_length[ match(1:(Ntip(x$tree)+Nnode(x$tree)),x$tree$edge[,1]) ]
#        node_and_edgeout_depth = node_depth + ifelse( is.na(node_edgeout_length), 0, node_edgeout_length )
#        tree_depth = max(node_and_edgeout_depth)
#        edge.length = tree_depth - node_depth[where]
#      }else{
#        stop("Must supply `edge.length`")
#      }
#    }
#
#    # build new tree
#    tree_new = phytools::bind.tip( x$tree,
#                                   tip.label = tip.label,
#                                   edge.length = edge.length,
#                                   where = where,
#                                   position = 0,
#                                   ... )
#
#    # refit
#    Args = as.list(x$call[-1])
#    Args$tree = tree_new
#    Args$startpar = x$opt$par
#    Args$quiet = TRUE
#    psem_new = do.call("phylosem", Args )
#
#    # extract
#    out = psem_new$report$x_vj[ which(psem_new$tree$tip.label == tip.label), ]
#  }
#
#  names(out) = colnames(x$data)
#  return( out )
#}
