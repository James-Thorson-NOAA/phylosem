
#' @title
#' Family for data that are known without error
#'
#' @description
#' Allows using \code{family = fixed()} to specify data that have no measurement error
#'
#' @export
fixed <- function() {
  link1 = "identity"
  l1 <- substitute(link1)
  if (!is.character(l1)) l1 <- deparse(l1)
  structure(
    list(
      link = l1,
      type = "fixed",
      family = "fixed",
      clean_name = "fixed"
    ),
    class = "family"
  )
}

#' @title
#' Family for categorical data
#'
#' @description
#' Allows using \code{family = categorical()} to specify factor-valued traits
#'
#' @param group string indicating which trait a given category belongs to
#'
#' @export
categorical <- function( group ) {
  link1 = "mlogit"
  l1 <- substitute(link1)
  if (!is.character(l1)) l1 <- deparse(l1)
  structure(
    list(
      link = l1,
      group = as.character(group),
      type = "categorical",
      family = "categorical",
      clean_name = "categorical"
    ),
    class = "family"
  )
}



