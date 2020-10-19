#' Calculate the meridional group velocity
#'
#' This function calculates the meridional group velocity.
#'
#' @param beta meridional gradient of the absolute vorticity
#'  in Mercator coordinates
# #' @param umz zonal wind in mercator coordinates
#' @param y The latitude y0 or y1
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
#' calcVg(beta = beta_y0,
#'        Ks2 = Ks2_y0,
#'        y = y0,
#'        lat = lat)
#' }
calcVg <- function(beta,
                   Ks2,
                   y,
                   lat,
                   K = 3,
                   direction = -1,
                   tl = 1,
                   a = 6371000) {
  k <- K/a
  l2 <- Ks2 - k^2
  l <- sqrt(abs(l2))

    vg <- (direction * tl*
             2 * beta * k * sqrt(abs(l2)) * cos(y*pi/180) / Ks2^2) *
          (180/(pi*a))

    return (vg)
}
