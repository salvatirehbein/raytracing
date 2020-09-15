#' Calculate the zonal group velocity
#'
#' This function calculates the zonal group velocity.
#'
#' @param betamz meridional gradient of the absolute vorticity
#'  in Mercator coordinates
#' @param umz zonal wind in mercator coordinates
#' @param y A latitude
#' @param lat numeric vector of latitudes from minor to major
#' @param K Total Rossby wavenumber
#' @param w wave frequency
#' @param a Earth's radio
#' @return The zonal group velocity
#' @export
#' @examples {
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' y0 <- -30
#' calcUg(betamz = colMeans(b$betam, na.rm = TRUE),
#'        umz = colMeans(b$um, na.rm = TRUE),
#'        y = y0,
#'        lat = b$lat) # 7.039724e-05
#' }
calcUg <- function(betamz,
                   umz,
                   y,
                   lat,
                   K = 3,
                   w = 0,
                   a = 6371000) {
  Ks2 <- betamz/umz
  k <- K/a
  cg <- w/k
  return (
    (cg + 2 * betamz[ypos(y, lat)] * k^2 /
       (Ks2[ypos(y, lat)]^2) ) * 180 / (a*pi)
  )
}
