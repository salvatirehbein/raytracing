#' Get the position of a given latitute y
#'
#' This function get the position in a vector of a given latitute y.
#' @family Interpolation
#' @param y numeric value of one latitude
#' @param lat numeric vector of latitudes from minor to major
#' @param yk numeric vector to be approximated
#' @param mercator Logical. Is it require to transform the final data
#' in mercator coordinates? Default is FALSE.
#' @return The position where the latitude y has the minor
#'  difference with lat
#' @export
#' @examples {
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(u = input)
#' ykk <- rev(colMeans(b$betam))
#' ypos(y = -30, lat = seq(90, -90, -2.5), yk = ykk)
#'
#' }
ypos <- function(y, lat, yk, mercator = FALSE){
  if(length(lat) != length(yk)) {
    stop("Length of lat and yk must be equal ")
  }
  df <- data.frame(ilat = seq_along(lat),
                   lat = lat,
                   y = y,
                   diff = abs(lat - y))
  x <- df[df$diff == min(df$diff, na.rm = TRUE), ]$ilat
  yy <- yk[x]
  if(mercator) yy <- yk[x]/cos(y*pi/180)
  return(yy)
}
