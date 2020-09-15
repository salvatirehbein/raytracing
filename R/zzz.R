.onAttach <- function(libname, pkgname) {
  msg <- paste0 ("Thanks for using the raytracing package!\n\n",
           " Please type: citation(\"raytracing\")\n for details on how to cite this package.\n")
  packageStartupMessage (msg)
}
