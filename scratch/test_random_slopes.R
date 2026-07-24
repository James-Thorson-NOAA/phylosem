

#pak::pak("james-thorson-NOAA/phylosem@dev")

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

# map$ln_theta = factor(c(1,1,1))
psem = phylosem::phylosem(
  data = data,
  sem = sem,
  estimate_ou = TRUE,
  tree = tree,
  estimate_xbar = colnames(data)
)

if( FALSE ){
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
}

#############
# Simulated example
#############

pak::pak( "James-Thorson-NOAA/phylosem@dev" )
set.seed(123)

library(phylosem)
library(phytools)
library(ape)

n_tips <- 100
tree <- pbtree(n = n_tips, scale = 1)  # phytools; scale=1 sets total tree depth to 1

x <- rTraitCont( tree, sigma = 0.5, alpha = 1, model = "BM")
w <- 1 + rTraitCont( tree, sigma = 0.5, alpha = 1, model = "BM")
y = 1 + w * x + rTraitCont( tree, sigma = 0.5, alpha = 1, model = "BM")

data = data.frame( x = x, y = y, w = NA )
sem = "
  x -> y, w
"

fit = phylosem(
  data = data,
  tree = tree,
  sem = sem,
  estimate_xbar = c("x", "y", "w"),
  estimate_ou = FALSE,
  control = phylosem_control(
    trace = 1
  )
)

cor( w, fit$parhat$x_vj[seq_len(n_tips),3] )
plot( w, fit$parhat$x_vj[seq_len(n_tips),3] )

