
library(phylosem)
library(ape)
#library(checkmate)

######################
# Morphometrics module
######################

data(FishBase_and_Morphometrics, package = "FishLife")

data = FishBase_and_Morphometrics$Y_ij[,1:17]
tree = FishBase_and_Morphometrics$tree

#sem = paste( FishBase_and_Morphometrics$SEM_model[1:21,1], collapse = "\n" )
sem = "
temperature -> log(length_infinity), b1
temperature -> log(growth_coefficient), b2
temperature -> log(natural_mortality), b3
temperature -> log(weight_infinity), b4
log(length_infinity) -> log(growth_coefficient), varying
log(length_infinity) -> log(natural_mortality), b6
log(length_infinity) -> log(length_max), b7
log(length_infinity) -> log(weight_infinity), b8
log(natural_mortality) -> log(length_maturity), b9
log(natural_mortality) -> log(age_maturity), b10
log(natural_mortality) -> log(age_max), b11
log(growth_coefficient) -> log(length_maturity), b12
log(growth_coefficient) -> log(age_maturity), b13
log(weight_infinity) -> trophic_level, b14
log(weight_infinity) -> log(fecundity), b15
log(weight_infinity) -> log(offspring_size), b16
log(length_infinity) -> log(aspect_ratio), b17
log(aspect_ratio) -> log(max_body_width), b18
log(aspect_ratio) -> log(max_body_depth), b19
log(aspect_ratio) -> log(lower_jaw_length), b20
log(aspect_ratio) -> log(min_caudal_pedoncule_depth), b21
"

# Check family
sd_j = apply(
  data,
  MARGIN = 2,
  FUN = \(x) max(tapply(x, INDEX = rownames(data), FUN = sd), na.rm=TRUE)
)
family = ifelse( sd_j == -Inf, "fixed", "normal")

fit = phylosem(
  sem = sem,
  tree = tree,
  data = cbind(data, varying=NA),
  family = family,
  control = phylosem_control(
    trace = 1,
    newton_loops = 0
  )
)

##############################
# SR module
##############################

#source( R'(C:\Users\jtuth\Documents\GitHub\phylosem\R\experiments.R)' )
#source( R'(C:\Users\jtuth\Documents\GitHub\phylosem\R\parse_path.R)' )

data( FishBase_and_RAM, package = "FishLife")
FishBase_and_RAM$StockData

Z_ik = FishBase_and_RAM$Z_ik
Z_ik$Species = paste( Z_ik$Genus, Z_ik$Species )
for( taxon in colnames(Z_ik) ) Z_ik[,taxon] = factor(Z_ik[,taxon])
tree = ape::as.phylo(
  ~Class/Order/Family/Genus/Species,
  data = Z_ik,
  collapse=FALSE
)
tree$edge.length = rep(1,nrow(tree$edge))
tree = collapse.singles(tree)
tree = root(tree, node=ape::Ntip(tree)+1 )

#
#ln_MASPS = FishBase_and_RAM$beta_gv[,'ln_MASPS']
#ln_M = FishBase_and_RAM$beta_gv[,'M']
#log_MLSPS = FishBase_and_RAM$beta_gv[,'ln_MASPS'] - log(1 - exp(-exp(ln_M)))
#h = exp(log_MLSPS) / (4 + exp(log_MLSPS) )
#
#exp( log_MLSPS[FishBase_and_RAM$StockData[,'Stock_to_i']] )

# Stock attributes
natural_mortality = FishBase_and_RAM$StockData[,'M']
unfished_spawners_pre_recruit = FishBase_and_RAM$StockData[,'SPRF0']
species = paste( FishBase_and_RAM$Z_ik[,'Genus'], FishBase_and_RAM$Z_ik[,'Species'] )[FishBase_and_RAM$StockData[,'Stock_to_i']]
names(natural_mortality) = names(unfished_spawners_pre_recruit) = names(species) = paste0( "stock_", seq_len(nrow(FishBase_and_RAM$StockData)) )

# SR observations
recruits = FishBase_and_RAM$SR_obs[,'R_obs']
spawners = FishBase_and_RAM$SR_obs[,'SSB_obs']
names(recruits) = names(spawners) = paste0( "stock_", FishBase_and_RAM$SR_obs[,'StockNum'] )

experiments = beverton_holt(
  recruits,
  spawners,
  unfished_spawners_pre_recruit,
  natural_mortality,
  species
)

data = FishBase_and_RAM$Y_ij[,c(1:8,11)]
data = setNames(data, c("logL_inf", "logK", "logW_inf", "logA_max", "logA_mat", "logM", "logL_mat", "T", "logMASPS") )

sem = "
  T -> logL_inf, b1

  logL_inf -> logM, b2
  logL_inf -> logK, b3
  logL_inf -> logW_inf, b4
  logL_inf -> logMASPS, b5

  logL_inf -> logL_mat, b6
  logM -> logL_mat, b7

  logM -> logA_max, b8

  logA_max -> logA_mat, b9
  logM -> logA_mat, b10
"

# Check family
#sd_j = apply(
#  data,
#  MARGIN = 2,
#  FUN = \(x) max(tapply(x, INDEX = Z_ik$Species, FUN = sd), na.rm=TRUE)
#)
#family = ifelse( sd_j == -Inf, "fixed", "normal")
family = rep("normal", ncol(data))

fit = phylosem(
  sem = sem,
  tree = tree,
  data = data,
  family = family,
  data_labels = Z_ik$Species,
  experiments = experiments,
  control = phylosem_control(
    trace = 1,
    getsd = TRUE,
    newton_loops = 0
  )
)


#
if(FALSE){
  #sem
  #tree
  #data
  family = rep("normal", ncol(data))
  covs = colnames(data)
  estimate_ou = FALSE
  estimate_lambda = FALSE
  estimate_kappa = FALSE
  data_labels = Z_ik$Species
  tmb_inputs = NULL
  estimate_xbar = NULL
  #experiments = NULL
  control = phylosem_control()

  setwd( R'(C:\Users\jtuth\Documents\GitHub\phylosem\src)' )
  dyn.unload( "phylosem" )
  TMB::compile("phylosem.cpp", framework = "TMBad" )
  dyn.load( "phylosem" )

  opt = nlminb(
    obj$par, obj$fn, obj$gr,
    control = list(trace = 1, iter.max = 1e5, eval.max = 1e5)
  )

  rep = obj$report()
  parhat = obj$env$parList()
  colnames(parhat$x_vj) = colnames(data)

  #
  logMLSPS = 1 + parhat$x_vj[,9] - log( 1 - exp(-1 * exp(parhat$x_vj[,6])) )
  h = exp(logMLSPS) / ( 4 + exp(logMLSPS) )

  #
  plot( x = rep$logmu_k, y = log(recruits) )

  # Plot old against new
  which_col = 9
  old = colnames(FishBase_and_RAM$beta_gv)[c(1:8,11)][which_col]
  new = colnames(data)[which_col]
  row_names = rownames(FishBase_and_RAM$beta_gv)
  row_names = sapply(
    row_names,
    \(x){
      y = strsplit(x, split = "_", fixed = TRUE)[[1]][4:5]
      paste( y[1], y[2] )
    }
  )
  match_rows = match( row_names, c(tree$tip.label,tree$node.label) )
  match_table = na.omit( cbind(seq_along(match_rows), match_rows) )
  plot(
    x = FishBase_and_RAM$beta_gv[match_table[,1],old],
    y = parhat$x_vj[match_table[,2],new]
  )
  abline( a = 0, b = 1, lty = "dotted" )
}


# Modeled traits:
# log_MASPS = model variable
# log_one_minus_expnegM = model variable
# log_SPR0 = model variable
# b = model variable
# mean(B/R) = model variable
# mean(B) = model variable

# Derived traits
# log_MLSPS = log( MASPS / (1-exp(-M)) ) = log_MASPS - log(1 - exp(-M)) = log_MASPS - log_one_minus_expnegM
# h = MLSPS / (4 + MLSPS)
  # a = MLSPS / SPR0
# log_a = log_MLSPS - log_SPR0
# logR = log_a + log(SSB / (1 + SSB/b )

# Input observations
# B_over_S = time-series data
# B = time-series data

# Estimated traits
# (ahat = 1 / (mean(B_t/R_t) - parlist$beta_z[1] * mean(B_t)))
# log_a = -log( mean(B/R) - beta * mean(B) )


