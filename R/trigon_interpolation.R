#' Trigonometric interpolation
#'
#' This function performs trigonometric interpolation for the
#' passed basic state variable and the requested latitude
#'
#' @param betamz meridional gradient of the absolute vorticity
#'  in Mercator coordinates
# #' @param umz zonal wind in mercator coordinates
#' @param y Numeric. The latitude where the interpolation is
#' required
#' @param yk Numeric vector of the data to be interpolated.
#' For instance, umz or betam
#' @param nlat Numeric. The total latitudes of the vector yk.
#' @return Numeric value
#' @export
#' @examples {
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' umz <- colMeans(b$um, na.rm = TRUE)
#' y0 <- -17
#' it(y = -17, yk = umz, nlat = length(b$lat))
#' }
it <- function(y,
               yk,
               nlat) {
  n <- nlat/2
  # Calculate the A parameter
  yk <- ifelse(yk >= 1e4, 0, yk)
  A <- yk[nlat]

  for (i in 1:n) {
    A <- sum(A, yk[i], yk[i + n], na.rm = TRUE)
  }

  # Calculate the ak and bk parameters
  l_ak <- list()
  l_bk <- list()
  for (i in 1:n) {
    ak <- 0
    bk <- 0
    for (j in 1:nlat) {
      mx <- 2*pi*i*(j-1)/(2*n + 1)
      ak <- sum(ak, yk[j]*cos(mx), na.rm = TRUE)
      bk <- sum(bk, yk[j]*sin(mx), na.rm = TRUE)
    }
    l_ak[[i]] <- ak
    l_bk[[i]] <- bk
  }

  # Great sum
  soma <- 0
  for (i in 1:n) {
    kx <- 2 * i * pi * (90 - y)/180
    soma <- sum(soma,
                l_ak[[i]] * cos(kx),
                l_bk[[i]] * sin(kx),
                na.rm = TRUE)
  }

  # Calculate the final interpolated value T(x)
  tx <- (2 * soma + A) / (2*n + 1)
  return(tx)
}
