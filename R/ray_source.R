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
#' @param interpolate Logical, use \code{\link{trin}}, TRUE or FALSE
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
#' b <- betaks(u = input)
#' r1 <- ray_source(betam = b$betam,
#'                  u = b$u,
#'                  lat = b$lat,
#'                  direction = -1,
#'                  x0 =  -130,#-seq(129,131, 0.5),
#'                  y0 = -17.06,#-seq(17,17.08, 0.01),
#'                  K = c(3),
#'                  dt = 12,
#'                  itime = 20*2, interpolation = "trin")
#' rl <- ray_path(r1$lon, r1$lat)
#' plot(r1$geometry, reset = F, axes = T)
#' plot(rl, col = "red", add = T)
#' # Simple plot:
#'# r1 <- r1[r1$iday %in% c(0,12),]
#' plot(r1["id"], axes = T, key.pos = 3)
#' plot(r1["id"],
#'      axes = TRUE,
#'      reset = FALSE,
#'      key.pos = 1,
#'      pal = c("red", "blue", "black"),
#'      pch = 16,
#'      key.length = 1)
#' lrs <- split(r1, r1$id)
#'
#' rl <- lapply(seq_along(lrs), function(i){
#'   x <- ray_path(x = lrs[[i]]$lon, y = lrs[[i]]$lat)
#' })
#' plot(rl[[1]], add = T, col = "red")
#' plot(rl[[2]], add = T, col = "blue")
#' plot(rl[[3]], add = T, col = "black")
#'
#'}
ray_source <- function(betam,
                       u,
                       lat,
                       x0,
                       y0,
                       K,
                       dt,
                       itime,
                       direction,
                       interpolation = "ypos",
                       tl = 1,
                       a = 6371000,
                       verbose = FALSE,
                       ofile){
  # combine the x0 and y0:
  df <- expand.grid(x0, y0)
  names(df) <- c("lon", "lat")
  dir <- direction
  wn <- K

  ddir <- lapply(seq_along(dir), function(i){
    ddf <- lapply(1:nrow(df), function(j){
      dwn <- lapply(1:length(wn), function(k){
        cbind(ray(betam = betam,
                  u = u,
                  lat = lat,
                  itime = itime,
                  K = wn[k],
                  dt = dt,
                  direction = dir[i],
                  interpolation = interpolation,
                  x0 = df$lon[j],
                  y0 = df$lat[j],
                  verbose = verbose),
              par = j, # de coordenadas x0 y0
              id = paste0("K:",wn[k],
                         ", (x0,y0):(",df$lon[j],",",df$lat[j],")"))
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
