

library(RTMB)
library(Matrix)
library(mvtnorm)
library(ape)
library(sem)
library(ggplot2)
library(gridExtra)

root_dir = R'(C:\Users\jtuth\Documents\GitHub\graphical_mixed_model)'
data_dir = file.path( root_dir, "data" )

# Functions
source( file.path(root_dir, "trait_functions.R") )

# Load data
tree = read.tree( file.path( data_dir, "VertTree_mammals.tre" ) )
max_edge = max(tree$edge.length)
tree$edge.length = tree$edge.length / max_edge
s_root = Ntip(tree) + 1
n_nodes = Nnode(tree)
n_tips = Ntip(tree)

traits = read.table( file.path(data_dir,"PanTHERIA_1-0_WR05_Aug2008.txt"),
            row.names = NULL )

tree$tip.label = tolower( tree$tip.label )
traits$traits_binom = paste0( tolower(traits$MSW05_Species), "_", traits$MSW05_Binomial )

trait_to_tip = match( traits$traits_binom, tree$tip.label )
traits = traits[ which(!is.na(trait_to_tip)), ]

add_NA = function(vec) ifelse( as.numeric(vec) == -999, NA, as.numeric(vec) )
data = data.frame(
  ln_metabolism = log(add_NA(traits[,'X18.1_BasalMetRate_mLO2hr'])),
  ln_range = log(add_NA(traits[,'X22.1_HomeRange_km2'])),
  ln_size = log(add_NA(traits[,'X5.1_AdultBodyMass_g']))
)
rownames(data) = traits$traits_binom

#
sem = "
  ln_size -> ln_metabolism, b1
  ln_size -> ln_range, b2
"
ou_j = c(TRUE, TRUE, TRUE)

#
plm = phylolm::phylolm(
  ln_metabolism ~ ln_size,
  data = data,
  phy = tree,
  model = "OUrandomRoot"
)

#############
# Fit in RTMB
#############

# Whether to drop root for ISAR
drop_bm_root = TRUE
which_drop = ((n_tips+n_nodes) * (seq_along(ou_j)-1) + (n_tips+1))[(ou_j==FALSE)]
assemble_version = 3

#
if( is.null(tree$node.label) & (n_nodes > 0) ){
  tree$node.label = paste0("node_",seq_len(n_nodes))
}

#
SEM_model = specifyModel( text=sem,
                          exog.variances=TRUE,
                          endog.variances=TRUE,
                          covs=colnames(data),
                          quiet=TRUE )
model = build_ram( model = SEM_model,
           vars = colnames(data) )

y_sj = as.matrix(data[match(c(tree$tip.label,tree$node.label),rownames(data)),,drop=FALSE])
parlist = list(
  y_sj = y_sj,
  beta_p = rep(1,max(model$parameter)),
  ln_theta = rep(0,ncol(data)),
  xbar = rep(0,ncol(data))
)
map = list(
  #ln_theta = factor( rep(1,ncol(data)) ),
  y_sj = ifelse( is.na(parlist$y_sj), seq_len(prod(dim(parlist$y_sj))), NA ),
  ln_theta = factor(ifelse(ou_j, seq_len(ncol(data)), NA))
)
if(isTRUE(drop_bm_root)) map$y_sj[which_drop] = NA
map$y_sj = factor(map$y_sj)

parlist$y_sj = ifelse( is.na(parlist$y_sj), 0, parlist$y_sj )

#method = "GMRF"
obj = MakeADFun( func = get_nll,
                  random = "y_sj",
                  map = map,
                  parameters = parlist,
                  silent = TRUE )
opt = nlminb( obj$par, obj$fn, obj$gr,
              control = list() )
rep = obj$report()
sdrep = sdreport(obj)

# Compare with phylosem ... logLik matches exactly when shared ln_theta
# map$ln_theta = factor(c(1,1,1))
psem = phylosem::phylosem(
  data = data,
  sem = sem,
  estimate_ou = all(ou_j),
  tree = tree
)

# Compare SDs
c( opt$par['beta_p'], psem$opt$par['beta_z'] )


#####################
# Experiment with moderated slopes
#####################

library(phylosem)

#
sem = "
  ln_size -> ln_metabolism, slope
  ln_size -> ln_range, b2
"

data$slope = NA

#sem
#tree
#data
family = rep("fixed", ncol(data))
covs = colnames(data)
estimate_ou = FALSE
estimate_lambda = FALSE
estimate_kappa = FALSE
data_labels = rownames(data)
tmb_inputs = NULL
control = phylosem_control()

