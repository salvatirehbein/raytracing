#' Calculate the meridional group velocity
#'
#' This function calculates the meridional group velocity.
#'
#' @param betamz meridional gradient of the absolute vorticity
#'  in Mercator coordinates
#' @param umz zonal wind in mercator coordinates
#' @param y A latitude
#' @param lat numeric vector of latitudes from minor to major
#' @param K Total Rossby wavenumber
#' @param tl Turning latitude. It will always start with a
#'  positive tl and automatically change to negative
#'  after the turning latitude or before a critical value.
#' @param direction Controls the wave displacement.
#' The wave goes to the north of the source using 1 and
#' to the south of the source using -1.
#' @param a Earth's radio
#' @return The meridional group velocity
#' @export
#' @examples {
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' y0 <- -30
#' calcVg(betamz = colMeans(b$betam, na.rm = TRUE),
#'        umz = colMeans(b$um, na.rm = TRUE),
#'        y = y0,
#'        lat = b$lat)
#' }
calcVg <- function(betamz,
                   umz,
                   y,
                   lat,
                   K = 3,
                   direction = -1,
                   tl = 1,
                   a = 6371000) {
  Ks2 <- betamz/umz
  k <- K/a
  l <- sqrt(abs(Ks2 - k^2))
  phirad <- lat*pi/180
  return (
    (direction*
      tl*
      2*
      betamz[ypos(y, lat = lat)]*
      k*
      l[ypos(y, lat = lat)]*
      cos(phirad[ypos(y, lat = lat)]) /
       Ks2[ypos(y, lat = lat)]^2) * (180/(pi*a))
  )
}
