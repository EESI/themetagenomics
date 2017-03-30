#' Prevent x renaming
#'
#' @param x x
#' @param value value
#'
#' @name prev_rename
NULL

#' @rdname prev_rename
#' @export
`names<-.topics` <- function(x,value){
  warning('topics-class xs cannot be renamed.')
  return(x)
}

#' @rdname prev_rename
#' @export
`attributes<-.topics` <- function(x,value){
  warning('topics-class attributes cannot be renamed.')
  return(x)
}

#' @rdname prev_rename
#' @export
`names<-.themetadata` <- function(x,value){
  warning('themetadata-class xs cannot be renamed.')
  return(x)
}

#' @rdname prev_rename
#' @export
`attributes<-.themetadata` <- function(x,value){
  warning('themetadata-class attributes cannot be renamed.')
  return(x)
}

#' @rdname prev_rename
#' @export
`names<-.functions` <- function(x,value){
  warning('functions-class xs cannot be renamed.')
  return(x)
}

#' @rdname prev_rename
#' @export
`attributes<-.functions` <- function(x,value){
  warning('functions-class attributes cannot be renamed.')
  return(x)
}

#' @rdname prev_rename
#' @export
`names<-.effects` <- function(x,value){
  warning('effects-class xs cannot be renamed.')
  return(x)
}

#' @rdname prev_rename
#' @export
`attributes<-.effects` <- function(x,value){
  warning('effects-class attributes cannot be renamed.')
  return(x)
}
