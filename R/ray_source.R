#' Calculate the Rossby waves ray paths over a source region
#'
#' This function calculates the Rossby waves ray paths based on
#' Hoskin and Ambrizzi (1993) in a source area of interest.
#' The source area is given bi the combination  of latitudes y0,
#' longitudes x0 and wave number K. This means that these numeric
#' vectors can have kength bigger than one.
#'
#' @param betamz matrix of (longitude = rows x latitude from minor
#' to major = columns) meridional gradient of the absolute vorticity
#' in Mercator coordinates. This package includes the
#' function \code{\link{betaks}} that import NetCDF files
#' using \strong{ncdf4}.
#' @param umz matrix of (longitude = rows x latitude from
#' minor to major = columns) zonal wind in mercator coordinates.
#'  This package includes the  function \code{\link{betaks}}
#'  that import NetCDF files using \strong{ncdf4}.
#' @param x0 Numeric value; longitude (0 to 360 degrees)
#' @param y0 Numeric value; latitude
#' @param lat Numeric vector of latitudes from minor to major
#'  (ex: -90 to 90)
#' @param K Numeric value; Total Rossby wavenumber
#' @param dt Numeric value; Timestep for integration (seconds)
#' @param itime Numeric value; total integration time
#' @param direction Numeric value (possibilities: 1 or -1).
#' It controls the wave displacement:
#' If 1, the wave goes to the north of the source;
#' If -1, the wave goes to the south of the source.
#' @param tl Numeric value; Turning latitude. Do not change this!
#' It will always start with a positive tl (1) and automatically
#'  change to negative (-1) after the turning latitude.
#' @param a Earth's radio (m)
#' @param ofile Character;  Output file name with .csv extension,
#' for instance, "/user/ray.csv"
#' @param verbose Boolean; if TRUE (default) return messages
#' during compilation
#' @return Rossby wave ray paths
#' @importFrom utils write.csv
#' @export
#' @examples {
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' a <- ray_source(betamz = colMeans(b$betam, na.rm = TRUE),
#'                 umz = colMeans(b$um, na.rm = TRUE),
#'                 lat = b$lat,
#'                 direction = -1)
#'
#' # Simple plot:
#' a$K <- factor(a$K)
#' library(ggplot2)
#' ww <- map_data('world', ylim=c(-60,60))
#'
#' ggplot(a, aes(x = x0 - 360,
#'               y = y0)) +
#'   geom_point(data = a[!(a$tun_y0 == -1 | a$tun_y1 == -1 | a$id == 0), ],
#'              aes(x = x0 - 360,
#'                  y = y0,
#'                  colour = K),
#'                  size = 3) +
#'  theme_bw() +
#'   geom_polygon(data = ww,
#'                aes(x = long,
#'                    y = lat,
#'                    group = group),
#'                    alpha = 0.0,
#'                    col = "grey")
#'}
ray_source <- function(betamz,
                       umz,
                       lat,
                       x0 = -(130:132) + 360,
                       y0 = -(16:18),
                       K = c(3, 4, 5),
                       dt = 12 * 60 * 60,
                       itime = 10,
                       direction = c(-1, 1),
                       tl = 1,
                       a = 6371000,
                       verbose = FALSE,
                       ofile){

  df <- expand.grid(x0, y0)
  names(df) <- c("lon", "lat")
  dir <- direction
  wn <- K

  ddir <- lapply(seq_along(dir), function(i){
    ddf <- lapply(1:nrow(df), function(j){
      dwn <- lapply(1:length(wn), function(k){
        cbind(ray(betamz = betamz,
                  umz = umz,
                  lat = lat,
                  itime = itime,
                  K = wn[k],
                  dt = dt,
                  direction = dir[i],
                  x0 = df$lon[j],
                  y0 = df$lat[j],
                  verbose = verbose),
              par = j)
      })
      do.call("rbind", dwn)
    })
    do.call("rbind", ddf)
  })

  DF <- as.data.frame(do.call("rbind", ddir))

  if(!missing(ofile)) {
    utils::write.csv(x = DF, file = ofile, row.names = FALSE)
  }
  return(DF)
}
