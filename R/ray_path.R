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
#' # Coelho et al. (2015):
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(u = input)
#' rt <- ray(betam = b$betam,
#'           u = b$u,
#'           lat = b$lat,
#'           K = 3,
#'           itime = 30,
#'           x0 = -135,
#'           y0 = -30,
#'           dt = 6 * 60 * 60,
#'           direction = -1)
#' plot(rt$geometry, axes = TRUE)
#' rp <- ray_path(x = rt$lon, y = rt$lat)
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
