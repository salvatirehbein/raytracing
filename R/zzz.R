.onAttach <- function(libname, pkgname) {
  msg <- paste0 ("\n")
  packageStartupMessage (msg)
}
