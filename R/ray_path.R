#' Calculate the ray paths / great circles
#'
#' This function calculates the great circles using the (x0, y0) obtained
#' with `ray` or `ray_source`. It returns a LINESTRING geometry that is ready
#' for plot.
#'
#' @param x vector with the longitude obtained with `ray` or `ray_source`
#' @param y vector with the latitude obtained with `ray` or `ray_source`
#' @return  sf data.frame
#' @importFrom sf st_linestring st_as_sfc st_segmentize st_shift_longitude
#' @importFrom units set_units
#' @export
#' @examples {
#' # For Magana and Ambrizzi (1998):
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_300hPa_1997-98DJF.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' ma98 <- ray(betamz = colMeans(b$betam, na.rm = TRUE),
#'                 umz = colMeans(b$um, na.rm = TRUE),
#'                 lat = b$lat,
#'                 K = 3,
#'                 itime = 10,
#'                 x0 = -130,
#'                 y0 = -17,
#'                 dt = 12 * 60 * 60,
#'                 direction = -1)
#'
#' # For Coelho et al. (2015):
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' co2015 <- ray(betamz = colMeans(b$betam, na.rm = TRUE),
#'           umz = colMeans(b$um, na.rm = TRUE),
#'           lat = b$lat,
#'           K = 3,
#'           itime = 30,
#'           x0 = -135 + 360,
#'           y0 = -30,
#'           dt = 6 * 60 * 60,
#'           direction = -1)
#'
#' r1 <- co2015
#' # We want to obtain all the linestrings except for the turning points
#' rp <- ray_path(x = r1$x0, y = r1$y0)
#' plot(rp, axes  = TRUE, col = "red")
#'
#'
#' ### ray source
#' r1 <- ray_source(betamz = colMeans(b$betam, na.rm = TRUE),
#'                 umz = colMeans(b$um, na.rm = TRUE),
#'                 lat = b$lat,
#'                 direction = -1,
#'                 x0 = -(130) + 360,
#'                 y0 = -(17),
#'                 K = c(3, 4, 5),
#'                 dt = 12 * 60 * 60,
#'                 itime = 30)
#' }
ray_path <- function(x,
                     y) {
  # For calculating Great Circle
  # TODO: remove sf dependence;


  geo <- cbind(x, y)

  dl <- lapply(1:(nrow(geo) - 1), function(i){
    sf::st_linestring(rbind(geo[i, ], geo[i+1, ]))
  })

  dl[[length(dl) + 1]] <- NA
  dls <- sf::st_as_sfc(dl, crs = 4326)

  # Calculate the Great Circle:
  # Entra com longitudes 0 a 360, mas apos aplicar st_segmentize
  # a longitude vai de -180 a 180 (centralizando em 0).
  dfl <- sf::st_segmentize(dls, units::set_units(2, "km"))

  # Voltar as longitudes de 0 a 360 para centralizar o plot no Pacifico
  # dfl <- sf::st_shift_longitude(x = dfl)
  # # plot(dls, axes  = T, graticule = T, lwd = 2)
  # plot(dfl, add  = T, col = "red")

  return(dfl)

}
