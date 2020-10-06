#' Calculate the meridional group velocity
#'
#' This function calculates the meridional group velocity.
#'
#' @param betamz meridional gradient of the absolute vorticity
#'  in Mercator coordinates
# #' @param umz zonal wind in mercator coordinates
#' @param y A latitude
#' @param lat Numeric vector of latitudes from minor to major
#' @param K Numeric. Total Rossby wavenumber
#' @param Ks2 Numeric. Stationary Rossby wavenumber
#' @param tl Numeric. Turning latitude. It will always start with a
#'  positive tl and automatically change to negative
#'  after the turning latitude or before a critical value.
#' @param direction Controls the wave displacement.
#' The wave goes to the north of the source using 1 and
#' to the south of the source using -1.
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
#' calcVg(betamz = betam_z,
#'        Ks2 = betam_z[ypos(y0, b$lat)]/um_z[ypos(y0, b$lat)],
#'        y = y0,
#'        lat = b$lat)
#' }
calcVg <- function(betamz,
                   # umz,
                   Ks2,
                   y,
                   lat,
                   K = 3,
                   direction = -1,
                   tl = 1,
                   a = 6371000) {
  # Ks2 <- betamz[ypos(y, lat = lat)]/umz[ypos(y, lat = lat)]
  k <- K/a
  l <- sqrt(abs(Ks2 - k^2))
  phirad <- lat*pi/180
  return (
    (direction*
      tl*
      2*
      betamz[ypos(y, lat = lat)]*
      k*
      l*
      cos(phirad[ypos(y, lat = lat)]) /
      Ks2^2) * (180/(pi*a))
  )
}
