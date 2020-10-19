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
#' b <- betaks(u = input)
#' ma98 <- ray(betam = b$betam,
#'             u = b$u,
#'             lat = b$lat,
#'             K = 3,
#'             itime = 10,
#'             x0 = -130,
#'             y0 = -17,
#'             dt = 12 * 60 * 60,
#'             direction = -1)
#' plot(ma98$geometry, axes = T)
#' rp <- ray_path(x = ma98$lon, y = ma98$lat)
#' plot(rp, add = TRUE)
#'
#' # For Coelho et al. (2015):
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(u = input)
#' co2015 <- ray(betam = b$betam,
#'               u = b$u,
#'               lat = b$lat,
#'               K = 3,
#'               itime = 30,
#'               x0 = -135,
#'               y0 = -30,
#'               dt = 6 * 60 * 60,
#'               direction = -1)
#' plot(co2015$geometry, axes = T)
#' rp <- ray_path(x = co2015$lon, y = co2015$lat)
#' plot(rp, add = TRUE)
#' }
ray_path <- function(x,
                     y) {
  geo <- cbind(x, y)

  dl <- lapply(1:(nrow(geo) - 1), function(i){
    sf::st_linestring(rbind(geo[i, ], geo[i+1, ]))
  })

  dl[[length(dl) + 1]] <- NA
  dls <- sf::st_as_sfc(dl, crs = 4326)

  # Calculate the Great Circle:
  dfl <- sf::st_segmentize(dls, units::set_units(2, "km"))

  return(dfl)

}
