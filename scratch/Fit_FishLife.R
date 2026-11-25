
#pak::local_install( R'(C:\Users\jtuth\Documents\GitHub\phylosem)' )

library(phylosem)
library(ape)
#library(checkmate)

######################
# Morphometrics module
######################

data(FishBase_and_Morphometrics, package = "FishLife")

data = FishBase_and_Morphometrics$Y_ij
tree = FishBase_and_Morphometrics$tree

#sem = paste( FishBase_and_Morphometrics$SEM_model[1:21,1], collapse = "\n" )
#sem = "
#temperature -> log(length_infinity), b1
#temperature -> log(growth_coefficient), b2
#temperature -> log(natural_mortality), b3
#temperature -> log(weight_infinity), b4
#log(length_infinity) -> log(growth_coefficient), b5
#log(length_infinity) -> log(natural_mortality), b6
#log(length_infinity) -> log(length_max), b7
#log(length_infinity) -> log(weight_infinity), b8
#log(natural_mortality) -> log(length_maturity), b9
#log(natural_mortality) -> log(age_maturity), b10
#log(natural_mortality) -> log(age_max), b11
#log(growth_coefficient) -> log(length_maturity), b12
#log(growth_coefficient) -> log(age_maturity), b13
#log(weight_infinity) -> trophic_level, b14
#log(weight_infinity) -> log(fecundity), b15
#log(weight_infinity) -> log(offspring_size), b16
#log(length_infinity) -> log(aspect_ratio), b17
#log(aspect_ratio) -> log(max_body_width), b18
#log(aspect_ratio) -> log(max_body_depth), b19
#log(aspect_ratio) -> log(lower_jaw_length), b20
#log(aspect_ratio) -> log(min_caudal_pedoncule_depth), b21
#
#log(length_infinity) -> spawning_typeguarders, b22
#log(length_infinity) -> spawning_typebearers, b23
#spawning_typebearers <-> spawning_typeguarders, cov1
#
#"

sem = FishBase_and_Morphometrics$text

# Check family
sd_j = apply(
  data,
  MARGIN = 2,
  FUN = \(x) max(tapply(x, INDEX = rownames(data), FUN = sd), na.rm=TRUE)
)
#family = ifelse( sd_j == -Inf, "fixed", "normal")
family = list(fixed(), gaussian())[ifelse( sd_j == -Inf, 1, 2 )]
names(family) = colnames(data)
for( i in 18:19 ){
  family[[i]] = categorical( group = "spawning_type" )
}
for( i in 20:23 ){
  family[[i]] = categorical( group = "habitat" )
}
for( i in 24:25 ){
  family[[i]] = categorical( group = "feeding_mode" )
}
for( i in 26:29 ){
  family[[i]] = categorical( group = "body_shape" )
}

fit = phylosem(
  sem = sem,
  tree = tree,
  data = cbind(data), #", varying=NA),
  family = family,
  control = phylosem_control(
    trace = 1,
    newton_loops = 0,
    nlminb_loops = 1,
    getsd = TRUE
  )
)

fit$parhat$x_vj

data.frame( FishBase_and_Morphometrics$SEM_model[1:76,1:2], L_z = FishBase_and_Morphometrics$ParHat$L_z )

#############################
# Morphometrics (small)
#############################

data(FishBase_and_Morphometrics, package = "FishLife")

#vars = c( "temperature", "log(length_infinity)", "feeding_modeplanktivorous_or_other", "feeding_modemacrofauna" )
vars = c(
  "temperature", "log(length_infinity)",
  #"spawning_typeguarders", "spawning_typebearers",
  #"habitatbathymetric", "habitatbenthopelagic", "habitatreefassociated", "habitatpelagic",
  #"feeding_modeplanktivorous_or_other", "feeding_modemacrofauna",
  "body_shapeelongated", "body_shapeshort_and_or_deep", "body_shapeeellike", "body_shapeother"
)
data = FishBase_and_Morphometrics$Y_ij[,vars]
tree = FishBase_and_Morphometrics$tree

# prune tree for speed
#tree = ape::collapse.singles(tree)
#tree = keep.tip(
#  phy = tree,
#  tip = intersect( tree$tip.label, rownames(data) )
#)

#sem = paste( FishBase_and_Morphometrics$SEM_model[1:21,1], collapse = "\n" )
sem = "
temperature -> log(length_infinity), b1
#temperature -> log(growth_coefficient), b2
#temperature -> log(natural_mortality), b3
#temperature -> log(weight_infinity), b4
#log(length_infinity) -> log(growth_coefficient), b5
#log(length_infinity) -> log(natural_mortality), b6
#log(length_infinity) -> log(length_max), b7
#log(length_infinity) -> log(weight_infinity), b8
#log(natural_mortality) -> log(length_maturity), b9
#log(natural_mortality) -> log(age_maturity), b10
#log(natural_mortality) -> log(age_max), b11
#log(growth_coefficient) -> log(length_maturity), b12
#log(growth_coefficient) -> log(age_maturity), b13
#log(weight_infinity) -> trophic_level, b14
#log(weight_infinity) -> log(fecundity), b15
#log(weight_infinity) -> log(offspring_size), b16
#log(length_infinity) -> log(aspect_ratio), b17
#log(aspect_ratio) -> log(max_body_width), b18
#log(aspect_ratio) -> log(max_body_depth), b19
#log(aspect_ratio) -> log(lower_jaw_length), b20
#log(aspect_ratio) -> log(min_caudal_pedoncule_depth), b21

#log(length_infinity) -> spawning_typeguarders, b22
#log(length_infinity) -> spawning_typebearers, b23
#spawning_typebearers <-> spawning_typeguarders, cov1

#log(length_infinity) -> habitatpelagic, c1
#log(length_infinity) -> habitatbenthopelagic, c2
#log(length_infinity) -> habitatreefassociated, c3
#log(length_infinity) -> habitatbathymetric, c4
#habitatbenthopelagic <-> habitatpelagic, c5
#habitatreefassociated <-> habitatpelagic, c6
#habitatbathymetric <-> habitatpelagic, c7
#habitatreefassociated <-> habitatbenthopelagic, c8
#habitatbathymetric <-> habitatbenthopelagic, c9
#habitatbathymetric <-> habitatreefassociated, c10

#log(length_infinity) -> feeding_modeplanktivorous_or_other, d1
#log(length_infinity) -> feeding_modemacrofauna, d2
#feeding_modemacrofauna <-> feeding_modeplanktivorous_or_other, d3

log(length_infinity) -> body_shapeelongated, e1
log(length_infinity) -> body_shapeother, e2
log(length_infinity) -> body_shapeshort_and_or_deep, e3
log(length_infinity) -> body_shapeeellike, e4
body_shapeother <-> body_shapeelongated, e5
body_shapeshort_and_or_deep <-> body_shapeelongated, e6
body_shapeeellike <-> body_shapeelongated, e7
body_shapeshort_and_or_deep <-> body_shapeother, e8
body_shapeeellike <-> body_shapeother, e9
body_shapeeellike <-> body_shapeshort_and_or_deep, e10
"

#sem = FishBase_and_Morphometrics$text

# Check family
n_continuous = 2
sd_j = apply(
  data[,1:n_continuous],
  MARGIN = 2,
  FUN = \(x) max(tapply(x, INDEX = rownames(data), FUN = sd), na.rm=TRUE)
)
#family = ifelse( sd_j == -Inf, "fixed", "normal")
family = list(fixed(), gaussian())[ifelse( sd_j == -Inf, 1, 2 )]
#for( i in (length(family)+1:2) ){
#  family[[i]] = categorical( group = "spawning_type" )
#}
#for( i in (length(family)+1:4) ){
#  family[[i]] = categorical( group = "habitat" )
#}
#for( i in (length(family)+1:2) ){
#  family[[i]] = categorical( group = "feeding_mode" )
#}
for( i in (length(family)+1:4) ){
  family[[i]] = categorical( group = "body_shape" )
}
names(family) = colnames(data)

fit = phylosem(
  sem = sem,
  tree = tree,
  data = cbind(data), #", varying=NA),
  family = family,
  control = phylosem_control(
    trace = 1,
    newton_loops = 0,
    nlminb_loops = 0,
    getsd = FALSE
  )
)

obj = fit$obj
opt = nlminb( obj$par, obj$fn, obj$gr, control = list(trace = 1))
rep = obj$report( obj$env$last.par.best )

H = obj$env$spHess( obj$env$last.par.best, random = TRUE )

fit$parhat$x_vj

data.frame( FishBase_and_Morphometrics$SEM_model[1:76,1:2], L_z = FishBase_and_Morphometrics$ParHat$L_z )

head( cbind(data, fit$rep$jnll_ig), 100 )


#############################
# SR module
#  Version saved in FishLife ...
#    D:\Work files (primary backup)\Collab-2017\2017 -- FishLife recruitment\2019-09-22\Taxon_v2_14_0\Nfac=-6_Nobs=0
#  Easier version ...
#    D:\Work files (primary backup)\Collab-2017\2017 -- FishLife recruitment\2019-09-22\Taxon_v2_14_0\Nfac=0_Nobs=0
##############################

#source( R'(C:\Users\jtuth\Documents\GitHub\phylosem\R\experiments.R)' )
#source( R'(C:\Users\jtuth\Documents\GitHub\phylosem\R\parse_path.R)' )

library(phylosem)
library(ape)
#library(phylosem)

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
unfished_spawners_per_recruit = FishBase_and_RAM$StockData[,'SPRF0']
species = paste( FishBase_and_RAM$Z_ik[,'Genus'], FishBase_and_RAM$Z_ik[,'Species'] )[FishBase_and_RAM$StockData[,'Stock_to_i']]
names(natural_mortality) = names(unfished_spawners_per_recruit) = names(species) = paste0( "stock_", seq_len(nrow(FishBase_and_RAM$StockData)) )

# SR observations
recruits = FishBase_and_RAM$SR_obs[,'R_obs']
spawners = FishBase_and_RAM$SR_obs[,'SSB_obs']
names(recruits) = names(spawners) = paste0( "stock_", FishBase_and_RAM$SR_obs[,'StockNum'] )

experiments = beverton_holt(
  species = species,
  recruits = recruits,
  spawners = spawners,
  unfished_spawners_per_recruit = unfished_spawners_per_recruit,
  natural_mortality = natural_mortality,
  variation_in_compensation = TRUE,      # if TRUE, variation goes to zero
  form = c("MAER", "MLER")[2],
  SD_logb = NA
)

data = FishBase_and_RAM$Y_ij[,c(1:8,11)]
# Same order as Collab-2017\2017 -- FishLife recruitment
data = setNames(data, c("logL_inf", "logK", "logW_inf", "logA_max", "logA_mat", "logM", "logL_mat", "T", "phi") )

# Simple SEM
sem = "
  T -> logL_inf, b1

  logL_inf -> logM, b2
  logL_inf -> logK, b3
  logL_inf -> logW_inf, b4
  logL_inf -> phi, b5

  logL_inf -> logL_mat, b6
  logM -> logL_mat, b7

  logM -> logA_max, b8

  logA_max -> logA_mat, b9
  logM -> logA_mat, b10
  logM -> phi, b11
"

# 3-factors with same order as FishLife
#sem = "
#  F1 -> logL_inf, F11
#  F2 -> logL_inf, NA, 0
#  F3 -> logL_inf, NA, 0
#
#  F1 -> logK, F12
#  F2 -> logK, F22
#  F3 -> logK, NA, 0
#
#  F1 -> logW_inf, F13
#  F2 -> logW_inf, F23
#  F3 -> logW_inf, F33
#
#  F1 -> logA_max, F14
#  F2 -> logA_max, F24
#  F3 -> logA_max, F34
#
#  F1 -> logA_mat, F15
#  F2 -> logA_mat, F25
#  F3 -> logA_mat, F35
#
#  F1 -> logM, F16
#  F2 -> logM, F26
#  F3 -> logM, F36
#
#  F1 -> logL_mat, F17
#  F2 -> logL_mat, F27
#  F3 -> logL_mat, F37
#
#  F1 -> T, F18
#  F2 -> T, F28
#  F3 -> T, F38
#
#  F1 -> logMASPSAR, F19
#  F2 -> logMASPSAR, F29
#  F3 -> logMASPSAR, F39
#
#  F1 <-> F1, NA, 1
#  F2 <-> F2, NA, 1
#  F3 <-> F3, NA, 1
#"
#data = cbind( data, F1 = NA, F2 = NA, F3 = NA )

# No factors, independent by variable
#sem = ""

# Check family
#sd_j = apply(
#  data,
#  MARGIN = 2,
#  FUN = \(x) max(tapply(x, INDEX = Z_ik$Species, FUN = sd), na.rm=TRUE)
#)
#family = ifelse( sd_j == -Inf, "fixed", "normal")
#family = ifelse( colnames(data) == "logMASPS", "fixed", "normal" )
#family = ifelse( colnames(data) == "logMASPS", "fixed", "normal" )
family = Map(function(.) gaussian(), colnames(data))
# logMASPS has variation among stocks within species in FishLife v2
family$phi = fixed()
# Factors are fixed
#family$F1 = fixed()
#family$F2 = fixed()
#family$F3 = fixed()

fit = phylosem(
  sem = sem,
  tree = tree,
  data = data,
  family = family,
  data_labels = Z_ik$Species,
  experiments = experiments,
  control = phylosem_control(
    trace = 1,
    getsd = TRUE,     # allow summary(.)
    newton_loops = 0,
    nlminb_loops = 1
  )
)

rep = fit$obj$report()
parhat = fit$obj$env$parList()
colnames(parhat$x_vj) = colnames(data)

# Inspect simple old
if(FALSE){
  load( R'(D:\Work files (primary backup)\Collab-2017\2017 -- FishLife recruitment\2019-09-22\Taxon_v2_14_0\Nfac=0_Nobs=0\Return.RData)' )
  cols = match( c("ln_MASPS", "M"), colnames(Return$Y_ij) )
  MLSPS = 1 + exp(Return$ParHat$beta_gj[,cols[1]]) / ( 1 - exp(-1 * exp(Return$ParHat$beta_gj[,cols[2]])) )
  h = MLSPS / ( 4 + MLSPS )
}

# Inspect old estimates
if(FALSE){
  load( R'(D:\Work files (primary backup)\Collab-2017\2017 -- FishLife recruitment\2019-09-22\Taxon_v2_14_0\Nfac=-6_Nobs=0\Return.RData)' )
  cols = match( c("ln_MASPS", "M"), colnames(Return$Y_ij) )
  MLSPS = 1 + exp(Return$ParHat$beta_gj[,cols[1]]) / ( 1 - exp(-1 * exp(Return$ParHat$beta_gj[,cols[2]])) )
  h = MLSPS / ( 4 + MLSPS )

  #
  summary(Return$ParHat$beta_gj[,cols])

  # Compare intercepts
  Return$ParHat$alpha_j
  parhat$x_vj[ ape::Ntip(tree)+1, ]

  # Compare stock-specific stuff
  Return$ParHat$bparam_stock
  fit$parhat$logb_p
}

#
cols = match( c("phi", "logM"), colnames(fit$data) )
MLSPS = 1 + exp(fit$parhat$x_vj[,cols[1]]) / ( 1 - exp(-1 * exp(fit$parhat$x_vj[,cols[2]])) )
h = MLSPS / ( 4 + MLSPS )
summary(fit$parhat$x_vj[,cols])

#
h_p = rep$MLRR_p / (4 + rep$MLRR_p)
mean(h_p)

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

########################
# Alternative
########################


library(phylosem)
library(ape)
library(rfishbase)
fishtable <- as.data.frame(load_taxa())

# Old tree
#Z_ik = unique(FishBase_and_RAM$Z_ik)
#Z_ik$Species = paste( Z_ik$Genus, Z_ik$Species )
##for( taxon in colnames(Z_ik) ) Z_ik[,taxon] = factor(Z_ik[,taxon])
##tree = ape::as.phylo(
##  ~Class/Order/Family/Genus/Species,
##  data = Z_ik,
##  collapse=FALSE
##)
##tree$edge.length = rep(1,nrow(tree$edge))
##tree = collapse.singles(tree)
##tree = root(tree, node=ape::Ntip(tree)+1 )
#syn <- synonyms(Z_ik$Species)
#species = syn$Species[match( Z_ik$Species, syn$synonym )]
#species_df = fishtable[ match(species, fishtable$Species), c("Class","Order","Family","Genus","Species") ]
#species_df = unique(species_df)
#for( ci in 1:5 ) species_df[,ci] = factor(species_df[,ci])
#tree = ape::as.phylo( ~Class/Order/Family/Genus/Species, data=species_df, collapse=FALSE)
#tree$edge.length = rep(1,nrow(tree$edge))
#tree = collapse.singles(tree)
#tree = root(tree, node=ape::Ntip(tree)+1 )

# Make tree for rfishbase taxonomy
colnames(fishtable) <- tolower(colnames(fishtable))
for( ci in c(2,3,5,6,7) ){
  fishtable[,ci] = ifelse(
    !is.na(as.numeric(fishtable[,ci])),
    NA,
    fishtable[,ci]
  )
}
fishtable = na.omit(fishtable[,c("species","genus","family","class","order")])
for( ci in 1:5 ) fishtable[,ci] = factor(fishtable[,ci])
tree = ape::as.phylo( ~class/order/family/genus/species, data=fishtable, collapse=FALSE)
tree$edge.length = rep(1,nrow(tree$edge))
tree = collapse.singles(tree)
tree = root(tree, node=ape::Ntip(tree)+1 )

# Stock attributes for SR analysis
data( FishBase_and_RAM, package = "FishLife")
natural_mortality = FishBase_and_RAM$StockData[,'M']
unfished_spawners_per_recruit = FishBase_and_RAM$StockData[,'SPRF0']
species = paste( FishBase_and_RAM$Z_ik[,'Genus'], FishBase_and_RAM$Z_ik[,'Species'] )[FishBase_and_RAM$StockData[,'Stock_to_i']]

# Fix names in SR database
syn <- synonyms(species)
species = syn$Species[match( species, syn$synonym )]
names(natural_mortality) = names(unfished_spawners_per_recruit) = names(species) = paste0( "stock_", seq_len(nrow(FishBase_and_RAM$StockData)) )

# SR observations
recruits = FishBase_and_RAM$SR_obs[,'R_obs']
spawners = FishBase_and_RAM$SR_obs[,'SSB_obs']
names(recruits) = names(spawners) = paste0( "stock_", FishBase_and_RAM$SR_obs[,'StockNum'] )

experiments = beverton_holt(
  species = species,
  recruits = recruits,
  spawners = spawners,
  unfished_spawners_per_recruit = unfished_spawners_per_recruit,
  natural_mortality = natural_mortality,
  variation_in_compensation = FALSE,      # if TRUE, variation goes to zero
  form = c("MAER", "MLER")[2],
  SD_logb = 10
)

data = cbind( species = Z_ik$Species, phi = unname(FishBase_and_RAM$Y_ij[,c(11),drop=FALSE]) )
# Same order as Collab-2017\2017 -- FishLife recruitment
family = list(phi = fixed())

# Subset
data = subset( data, data$species %in% tree$tip.label )

# Simple SEM
sem = "
"

fit = phylosem(
  sem = sem,
  tree = tree,
  data = data[,2,drop=FALSE],
  family = family,
  data_labels = data$species,
  experiments = experiments,
  estimate_ou = FALSE,
  estimate_lambda = TRUE,
  control = phylosem_control(
    trace = 1,
    getsd = FALSE,     # allow summary(.)
    newton_loops = 0,
    nlminb_loops = 1
  )
)

rep = fit$obj$report()
parhat = fit$obj$env$parList()
#colnames(parhat$x_vj) = colnames(data)

#
cols = match( c("phi"), colnames(fit$data) )
MLSPS = 1 + exp(fit$parhat$x_vj[,cols[1]])
h = MLSPS / ( 4 + MLSPS )
mean(h)

table(names(fit$opt$par))
fit$opt$par[c("logitlambda","logalpha")]

Match = match( Z_ik$Species, tree$tip.label )
mean(h[Match], na.rm=TRUE)
