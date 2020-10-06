#' Calculate the zonal group velocity
#'
#' This function calculates the zonal group velocity.
#'
#' @param betamz meridional gradient of the absolute vorticity
#'  in Mercator coordinates
# #' @param umz zonal wind in mercator coordinates
#' @param y A latitude
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
#' b <- betaks(ifile = input)
#' betam_z <- colMeans(b$betam, na.rm = TRUE)
#' um_z <- colMeans(b$um, na.rm = TRUE)
#' y0 <- -30
#' calcUg(betamz = betam_z,
#'        Ks2 = betam_z[ypos(y0, b$lat)]/um_z[ypos(y0, b$lat)],
#'        y = y0,
#'        lat = b$lat) # 7.039724e-05
#' }
calcUg <- function(betamz,
                   # umz,
                   Ks2,
                   y,
                   lat,
                   K = 3,
                   a = 6371000) {
  k <- K/a
  # Ks2 <- betamz[ypos(y, lat)]/umz[ypos(y, lat)]

  ug <-  (2 * betamz[ypos(y, lat)] * k^2 / (Ks2^2) ) * 180 / (a*pi)

  return (ug)
}
