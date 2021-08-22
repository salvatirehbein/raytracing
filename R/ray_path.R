#' Calculate the ray paths / segment of great circles
#'
#' This function calculates the segments great circles using the (lat, lon)
#' coordinates obtained with `ray` or `ray_source`.
#' It returns a LINESTRING geometry that is ready for plot.
#'
#' @param x vector with the longitude obtained with `ray` or `ray_source`
#' @param y vector with the latitude obtained with `ray` or `ray_source`
#' @return  sfc_LINESTRING sfc
#' @importFrom sf st_linestring st_as_sfc st_segmentize st_sfc
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
#'           dt = 6,
#'           direction = -1)
#' rp <- ray_path(x = rt$lon, y = rt$lat)
#' plot(rp,  axes = TRUE, graticule = TRUE)
#' }
ray_path <- function(x,
                     y) {
  geo <- cbind(x, y)

  dl <- lapply(1:(nrow(geo) - 1), function(i){
    sf::st_linestring(rbind(geo[i, ], geo[i+1, ]))
  })

  dl[[length(dl) + 1]] <- NA
  dls <- sf::st_as_sfc(dl, crs = 4326)

  # Calculate the segment of Great Circle
  dfl <- sf::st_segmentize(dls, units::set_units(20, "km"))
  dfl <- sf::st_wrap_dateline(dfl)

  # # Obtain multilinestring for west to east hemispheric waves
  # co <- as.data.frame(sf::st_coordinates(dfl))
  # names(co) <- c("xini", "yini", "id")
  # co$xend <- co$xini[c(2:length(co$xini), NA)]
  # co$yend <- co$yini[c(2:length(co$yini), NA)]
  # co2 <- co[-nrow(co),]
  #
  # # Detect where the start long is greater than end long:
  # l <- lapply(1:nrow(co2), function(i){
  #   ifelse(co2$xini[i] < co2$xend[i],
  #          paste0(co2$xini[i],
  #                 " ",
  #                 co2$yini[i],
  #                 ", ",
  #                 co2$xend[i],
  #                 " ",
  #                 co2$yend[i]),
  #          paste0("AQUI ")
  #   )
  # })
  #
  # # create the multilinestring feature:
  # ll <- c(paste0(unlist(l)[1:nrow(co2)-1], ","),
  #         paste0(l[nrow(co2)], ")"))
  # ll <- paste(unlist(ll), collapse=" ")
  # ll <- gsub(pattern = ", AQUI , ", replacement = "), (", x = ll)
  #
  # sdf <- sf::st_sfc(geometry = sf::st_as_sfc(
  #                    paste0("MULTILINESTRING ((",ll, ")")),
  #                  crs = 4326)
  #
  return(dfl)

}
