
library(ape)
library(phylosem)

setwd( R'(C:\Users\jtuth\Documents\GitHub\phylosem)' )

# https://github.com/tessington/metabolic_index/tree/main/data
all.dat = readRDS( file.path("data-raw", "alldata_taxonomy.RDS") )

# https://github.com/tessington/metabolic_index/blob/main/code/02-fit_tmb_model.R
kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W/wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$Pcrit_atm<- all.dat$Pcrit / 101.325 # convert from kPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit)

### Create new ParentChild matrix for reduced taxonomic structure ####
taxa.list <- c("Phylum", "Class","Order", "Family", "Genus", "Species")
#taxa.info <- make_taxa_tree(all.dat, taxa.list)
for( taxon in taxa.list ) all.dat[,taxon] = factor(all.dat[,taxon])
tree = ape::as.phylo(
  ~Phylum/Class/Order/Family/Genus/Species,
  data = all.dat,
  collapse=FALSE
)
tree$edge.length = rep(1,nrow(tree$edge))
tree = collapse.singles(tree)
tmp = root(tree, node=ape::Ntip(tree)+1 )

invtemp = all.dat$inv.temp
logW = log(all.dat$W)
minuslogpo2 = -log(all.dat$Pcrit)


#        if (j == 0) {
#          spc_ij(i,j) = exp(beta_gj(spc_in_PCgz( i ), j ));
#        }
#        if (j >0) {
#          spc_ij(i,j) = beta_gj(spc_in_PCgz( i ), j );
#        }
# V = spc_ij.col(0);
#    n_pow = spc_ij.col(1);
#    Eo = spc_ij.col(2);

# https://github.com/tessington/metabolic_index/blob/main/code/TMB/hierarchical_mi_base.cpp#L50-L55
# mu( id ) =  Eo( taxa_id( id ) ) * invtemp( id ) + n_pow( taxa_id( id ) ) * logW( id  ) - log(V( taxa_id( id ) ))
# jnll_comp( 1 ) = -sum( dnorm( minuslogpo2, mu, sigma, true) );

data = data.frame(
  invtemp = invtemp,
  logW = logW,
  minuslogpo2 = minuslogpo2,
  n_pow = NA,
  Eo = NA
)

sem = "
  logW -> minuslogpo2, n_pow
  invtemp -> minuslogpo2, Eo
"

fit = phylosem(
  data = data,
  sem = sem,
  tree = tree,
  data_labels = as.character(all.dat$Species),
  family = c("normal", "normal", "normal", "fixed", "fixed")
)

