#' Get the position of a given latitute y
#'
#' This function get the position in a vector of a given latitute y.
#'
#' @param y numeric value of one latitude
#' @param lat numeric vector of latitudes from minor to major
#' @return The position where the latitude y has the minor
#'  difference with lat
#' @export
#' @examples {
#' ypos(y = -30, lat = seq(-90, 90, 2.5))
#' }
ypos <- function(y, lat){
  df <- data.frame(ilat = seq_along(lat),
                   lat = lat,
                   y = y,
                   diff = abs(lat - y))
  df[df$diff == min(df$diff, na.rm = TRUE), ]$ilat
}
