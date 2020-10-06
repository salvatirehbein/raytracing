#' Calculate the Rossby waves ray paths over a source region
#'
#' \code{ray_source} must ingest the zonal means of \strong{betam} and
#' \strong{um}, plus the latitude vector (\strong{lat}). \code{ray_source}
#' returns the Rossby wave ray paths (lat/lon) triggered according a set of
#' experiments. The ray paths will be calculated for the combination
#' of initial sources/positions (x0, y0), total wavenumbers (K),
#' and directions, using the same basic state given by  \code{\link{betaks}}.
#'
#' @param betamz zonal mean of \strong{betam}. \strong{betam} is a
#' matrix (longitude = rows x latitude from minor to major = columns)
#' obtained with \code{\link{betaks}}. It is suggested to use
#' \code{\link{colMeans}} base R function in order to obtain the
#' zonal mean of \strong{betam}.
#' @param umz  zonal mean of \strong{um}. \strong{um} is a
#' matrix (longitude = rows x latitude from minor to major = columns)
#' obtained with \code{\link{betaks}}. It is suggested to use
#' \code{\link{colMeans}} base R function in order to obtain the
#' zonal mean of \strong{um}.
#' @param x0 Vector with longitudes (0 to 360 degrees)
#' @param y0 Vector with latitudes
#' @param lat Numeric vector of latitudes from minor to major
#'  (ex: -90 to 90)
#' @param K Vector; Total Rossby wavenumber
#' @param dt Numeric value; Timestep for integration (seconds)
#' @param itime Numeric value; total integration time
#' @param direction Vector with two possibilities: 1 or -1
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
#' @return  sf data.frame
#' @importFrom utils write.csv
#' @export
#' @examples {
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_300hPa_1997-98DJF.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' r1 <- ray_source(betamz = colMeans(b$betam, na.rm = TRUE),
#'                 umz = colMeans(b$um, na.rm = TRUE),
#'                 lat = b$lat,
#'                 direction = -1,
#'                 x0 = -130 + 360,
#'                 y0 = -17,
#'                 K = c(5),
#'                 dt = 12 * 60 * 60,
#'                 itime = 10*2)
#' # Simple plot:
#'
#'}
ray_source <- function(betamz,
                       umz,
                       lat,
                       x0,
                       y0,
                       K,
                       dt,
                       itime,
                       direction,
                       tl = 1,
                       a = 6371000,
                       verbose = FALSE,
                       ofile){
  # faz combinacao dos x0 e y0
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
              par = j # de coordenadas x0 y0
              )
      })
      do.call("rbind", dwn)
    })
    do.call("rbind", ddf)
  })

  DF <- do.call("rbind", ddir)

  if(!missing(ofile)) {
    utils::write.csv(x = DF, file = ofile, row.names = FALSE)
  }
  return(DF)
}
