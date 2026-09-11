
#pak::pak( "James-Thorson-NOAA/phylosem@dev" )
#remotes::install_github( "James-Thorson-NOAA/phylosem@dev" )
#set.seed(123)

library(phylosem)
library(phytools)
library(ape)

n_tips <- 100
tree <- pbtree(n = n_tips, scale = 1)  # phytools; scale=1 sets total tree depth to 1

# stocks_sim ... no length species -> stock (used for estimation)
stocks_est = tree
for( ti in seq_len(n_tips) ){
  stocks_est = bind.tip(
    stocks_est,
    tip.label = paste0("t",ti,".1"),
    edge.length = 0,
    where = which(stocks_est$tip.label == paste0("t",ti))
  )
}

# stocks_est ... yes length species -> stock (used for simulation)
stock_sd = 0.1
stocks_sim = stocks_est
stocks_sim$edge.length = ifelse( stocks_sim$edge.length==0, stock_sd, stocks_sim$edge.length )
node.depth.edgelength(stocks_sim)

# Visualize
setNames(
  data.frame(stocks_sim$edge,stocks_sim$edge.length,stocks_est$edge.length),
  c("parent", "child", "length_sim", "length_est")
)

x <- rTraitCont( stocks_sim, sigma = 0.5, alpha = 1, model = "BM")
y = 1 + x + rTraitCont( stocks_sim, sigma = 0.5, alpha = 1, model = "BM")

data = data.frame( x = x, y = y )
sem = "
  x -> y, w
"

fit = phylosem(
  data = data,
  tree = stocks_est,
  sem = sem,
  estimate_xbar = c("x", "y"),
  estimate_lambda = TRUE,
  #estimate_ou = TRUE,
  control = phylosem_control(
    trace = 1
  )
)

1 - plogis(fit$opt$par['logitlambda'])


