#' Prevent object renaming in class topics
#' @export
`names<-.topics` <- function(object,value){
  warning('topics-class objects cannot be renamed.')
  return(object)
}

#' Prevent attribute renaming in class topics
#' @export
`attributes<-.topics` <- function(object,value){
  warning('topics-class attributes cannot be renamed.')
  return(object)
}

#' Prevent object renaming in class themetadata
#' @export
`names<-.themetadata` <- function(object,value){
  warning('themetadata-class objects cannot be renamed.')
  return(object)
}

#' Prevent attribute renaming in class themetadata
#' @export
`attributes<-.themetadata` <- function(object,value){
  warning('themetadata-class attributes cannot be renamed.')
  return(object)
}

#' Prevent object renaming in class functions
#' @export
`names<-.functions` <- function(object,value){
  warning('functions-class objects cannot be renamed.')
  return(object)
}

#' Prevent attribute renaming in class functions
#' @export
`attributes<-.functions` <- function(object,value){
  warning('functions-class attributes cannot be renamed.')
  return(object)
}

#' Prevent object renaming in class effects
#' @export
`names<-.effects` <- function(object,value){
  warning('effects-class objects cannot be renamed.')
  return(object)
}

#' Prevent attribute renaming in class topics
#' @export
`attributes<-.effects` <- function(object,value){
  warning('effects-class attributes cannot be renamed.')
  return(object)
}
