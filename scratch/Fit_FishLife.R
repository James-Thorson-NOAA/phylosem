
library(phylosem)
data(FishBase_and_Morphometrics, package = "FishLife")

data = FishBase_and_Morphometrics$Y_ij[,1:17]
tree = FishBase_and_Morphometrics$tree

#sem = paste( FishBase_and_Morphometrics$SEM_model[1:21,1], collapse = "\n" )
sem = "
temperature -> log(length_infinity), b1
temperature -> log(growth_coefficient), b2
temperature -> log(natural_mortality), b3
temperature -> log(weight_infinity), b4
log(length_infinity) -> log(growth_coefficient), b5
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
  data = data,
  family = family,
  control = phylosem_control(
    trace = 1,
    newton_loops = 0
  )
)

