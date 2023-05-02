#' Convert variable names to a factor-model specification
#'
variables_to_factors <-
function( var_names,
          n_factors,
          method,
          small_var = 0.001,
          concatenate = TRUE,
          estimate_diagonal = FALSE ){

  text = NULL

  # Method-1
  if(method=="covariance"){
    for(i in seq_len(n_factors) ){
    for(j in i:length(var_names)){
      string = paste0(
        var_names[i],
        " <-> ",
        var_names[j],
        ", cov_",
        var_names[i],
        "_",
        var_names[j],
        ", 1"
      )
      #if(text==""){
      #  text = string
      #}else{
      #  text = paste0( text, " \n ", string)
      #}
      text = c( text, string )
    }}
    for(i in (n_factors+seq_len(length(var_names)-n_factors)) ){
      string = paste0(
        var_names[i],
        " <-> ",
        var_names[i],
        #", proces_var",
        #", ", 0.1
        ", ", NA,
        ", ", small_var
      )
      #if(text==""){
      #  text = string
      #}else{
      #  text = paste0( text, " \n ", string)
      #}
      text = c( text, string )
    }
  }

  # Method-2
  if(method=="factors"){
    for(i in seq_len(n_factors) ){
    for(j in i:length(var_names)){
      string = paste0(
        "factor",i,
        " -> ",
        var_names[j],
        ", f",i,"_",
        var_names[j],
        ", 1"
      )
      #if(text==""){
      #  text = string
      #}else{
      #  text = paste0( text, " \n ", string)
      #}
       text = c( text, string )
    }}
    for(i in seq_len(n_factors) ){
      string = paste0(
        "factor",i,
        " <-> ",
        "factor",i,
        ", ", NA,
        ", 1"
      )
      text = c( text, string )
    }
    for(i in seq_along(var_names) ){
      if( estimate_diagonal==FALSE ){
        string = paste0(
          var_names[i],
          " <-> ",
          var_names[i],
          ", ", NA,
          ", ", small_var
        )
      }else{
        string = paste0(
          var_names[i],
          " <-> ",
          var_names[i],
          ", sd", i
        )
      }
      #if(text==""){
      #  text = string
      #}else{
      #  text = paste0( text, " \n ", string)
      #}
      text = c( text, string )
    }
  }
  if(is.null(text)) text = ""
  if( concatenate==TRUE ){
    text = paste0( text, collapse=" \n ")
  }

  return(text)
}
