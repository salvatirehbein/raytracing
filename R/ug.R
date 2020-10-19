#' Calculate the zonal group velocity
#'
#' This function calculates the zonal group velocity.
#'
#' @param beta meridional gradient of the absolute vorticity
#'  in Mercator coordinates
# #' @param umz zonal wind in mercator coordinates
#' @param lat Numeric vector of latitudes from minor to major
#' @param K Numeric. Total Rossby wavenumber
#' @param Ks2 Numeric. Stationary Rossby wavenumber
#' @param a Earth's radio
#' @return Numeric.
#' @export
#' @examples {
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(u = input)
#' y0 <- -30
#' lat <- rev(b$lat)
#' phirad <- lat*pi/180
#' betamz <- rev(colMeans(b$betam, na.rm = TRUE))
#' umz <- rev(colMeans(b$u, na.rm = TRUE))*cos(phirad)
#' beta_y0 <- trin(y0, yk = betamz)
#' u_y0 <- trin(y0, yk = umz)
#' Ks2_y0 <- beta_y0/u_y0
#'
#' calcUg(beta = beta_y0,
#'        Ks2 = Ks2_y0,
#'        lat = lat)
#' }
calcUg <- function(beta,
                   Ks2,
                   lat,
                   K = 3,
                   a = 6371000) {

  k <- K/a

  ug <-  (2 * beta * k^2 / Ks2^2 ) * 180 / (a*pi)

  return (ug)
}
