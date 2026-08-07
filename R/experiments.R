
beverton_holt =
function( recruits,
          spawners,
          unfished_spawners_pre_recruit,
          natural_mortality,
          species ){

  assertNumeric( recruits, names = "named", finite = TRUE )
  assertNumeric( spawners, names = "named", len = length(recruits), finite = TRUE )
  assertNumeric( unfished_spawners_pre_recruit, names = "named", finite = TRUE )
  assertNumeric( natural_mortality, names = "named", len = length(unfished_spawners_pre_recruit), finite = TRUE )
  assertCharacter( species, names = "named", len = length(unfished_spawners_pre_recruit) )

  time_data = data.frame(
    population = names(recruits),
    R = recruits,
    S = spawners
  )

  population_data = merge(
    data.frame( population = names(unfished_spawners_pre_recruit), SPR0 = unfished_spawners_pre_recruit),
    data.frame( population = names(natural_mortality), M = natural_mortality)
  )
  population_data = merge(
    data.frame( population = names(species), species = species),
    population_data
  )

  out = list(
    type = "BH",
    time_data = time_data,
    population_data = population_data
  )
  return(out)
}

#linear =
#function( formula,
#          data ){
#
#  #
#  model.matrix( formula, data = data )
#
#}
