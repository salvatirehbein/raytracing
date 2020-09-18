.onAttach <- function(libname, pkgname) {
  msg <- paste0 ("Thanks for using the raytracing package!\n\n",
           "Though times through COVID-19 pandemy trying to develop this tool. So, please I would really appreciate any citation when using this package. Many thanks! Tip: citation(\"raytracing\")")
  packageStartupMessage (msg)
}
