#' Filter the ray paths that arrives in an area of interest
#'
#' \code{wave_arrival} ingests the ray paths to filter by determined area of
#' interest. Default CRS 4326.
#'
#' @param x sf data.frame object with the LINESTRINGS to be filtered.
#' @param aoi String giving the path and the filename of the area of interest.
#' By default is NULL. If no \strong{aoi} is not provided, the xmin, xmax, ymin,
#' and ymax must be provided.
#' @param xmin Numeric. Indicates the western longitude to be used in the
#' range -180 to 180.
#' @param xmax Numeric. Indicates the eastern longitude to be used in the
#' range -180 to 180.
#' @param ymin Numeric. Indicates the southern longitude to be used in the
#' range -90 to 90.
#' @param ymax Numeric. Indicates the northern longitude to be used in the
#' range -90 to 90.
#' @param ofile Character;  Output file name with .csv extension,
#' for instance, "/user/aoi_ray.csv"
#' @return sf data.frame
#' @importFrom sf st_read st_polygon st_as_sfc st_as_text st_sf
#' st_intersects
#' @importFrom utils write.csv
#' @export
#' @examples {

#' }
wave_arrival <- function(x,
                         aoi = NULL,
                         xmin,
                         xmax,
                         ymin,
                         ymax,
                         ofile){

  if(!is.null(aoi)){
    region <- sf::st_read(aoi)
    plot(region$geometry, axes = TRUE)
  } else {
    # define area of interest:
    caixa_p <- sf::st_polygon(
      list(rbind(c(xmin, ymin),
                 c(xmin, ymax),
                 c(xmax, ymax),
                 c(xmax, ymin),
                 c(xmin, ymin)))
    )
    caixa_sf <- sf::st_as_sfc(
      sf::st_as_text(caixa_p),
      crs = 4326
    )
    region <- sf::st_as_sf(caixa_sf)
    plot(region$x, axes = TRUE)
  }

  # Inicio do filtro ####
  arqi_sf <- sf::st_sf(x, geometry = x$geometry)

  # Passa ou nao pela region:
  xf <- sf::st_intersects(
    arqi_sf,
    region,
    sparse = FALSE
  )
  arqi_sf$is_in_region <- xf[, 1]
  ids_na_caixa <- unique(arqi_sf[arqi_sf$is_in_region == T, ]$id)

  dfinal <- arqi_sf[arqi_sf$id %in% ids_na_caixa, ]
  if(!missing(ofile)) {
    write.csv(x = dfinal, file = ofile)
  }
  return(dfinal)
}
