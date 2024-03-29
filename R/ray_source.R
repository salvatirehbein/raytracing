#' Calculate the Rossby waves ray paths over a source region
#'
#' \code{ray_source} returns the Rossby wave ray paths (lat/lon) triggered from
#' one or more initial source/position (x0, y0), one or more total
#' wavenumber (K), and one or more direction set up when invoking the function.
#' \code{ray_source} must ingest the meridional gradient of the absolute
#' vorticity in mercator coordinates\strong{betam}, the zonal mean wind
#' \strong{u}, and the latitude vector (\strong{lat}). Those variables can be
#' obtained (recommended) using \code{\link{betaks}} function. The zonal means
#' of the basic state will be calculated along the \strong{ray} program, as well
#' as the conversion to mercator coordinates of \strong{u}.
#' The resultant output is a spatial feature object from a combination of
#' initial and final positions/sources, total wavenumbers (K), and directions.
#'
#' @param betam matrix (longitude = rows x latitude from minor to
#' major = columns) obtained with \code{\link{betaks}}. \strong{betam} is the
#' meridional gradient of the absolute vorticity in mercator coordinates
#' @param u matrix (longitude = rows x latitude from minor to
#' major = columns) obtained with \code{\link{betaks}}. Is the zonal wind speed
#' in the appropriate format for the \code{ray}. It will be converted in mercator
#' coordinates inside the \code{ray}
#' @param x0 Vector with the initial longitudes (choose between -180 to 180)
#' @param y0 Vector with the initial latitudes
#' @param lat Numeric vector of latitudes from minor to major
#'  (ex: -90 to 90). Obtained with \code{\link{betaks}}
#' @param K Vector; Total Rossby wavenumber
#' @param cx numeric. Indicates the zonal phase speed. The program is designed
#' for eastward propagation (cx > 0) and stationary waves (cx = 0, the default).
#' @param dt Numeric value; Timestep for integration (hours)
#' @param itime Numeric value; total integration time. For instance, 10 days
#' times 4 times per day
#' @param direction Vector with two possibilities: 1 or -1
#' It controls the wave displacement:
#' If 1, the wave goes to the north of the source;
#' If -1, the wave goes to the south of the source.
#' @param interpolation Character. Set the interpolation method to be used:
#' \code{\link{trin}} or \code{\link{ypos}}
#' @param tl Numeric value; Turning latitude. Do not change this!
#' It will always start with a positive tl (1) and automatically
#'  change to negative (-1) after the turning latitude.
#' @param a Earth's radio (m)
#' @param ofile Character;  Output file name with .csv extension,
#' for instance, "/user/ray.csv"
#' @param verbose Boolean; if TRUE (default) return messages
#' during compilation
#' @return  sf data.frame
#' @importFrom sf st_as_sf st_set_geometry
#' @importFrom utils write.csv
#' @export
#' @examples \dontrun{
#' #do not run
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(u = input)
#' rt <- ray_source(betam = b$betam,
#'                  u = b$u,
#'                  lat = b$lat,
#'                  K = 3,
#'                  itime = 10*4,
#'                  cx = 0,
#'                  x0 = -c(130, 135),
#'                  y0 = -30,
#'                  dt = 6,
#'                  direction = -1,
#'                  interpolation = "trin")
#'
#' # Plot:
#' data(coastlines)
#' plot(coastlines,
#'      reset = FALSE,
#'      axes = TRUE,
#'      graticule = TRUE,
#'      col = "grey",
#'      main = "Coelho et al. (2015): JFM/2014")
#' plot(rt[sf::st_is(rt, "LINESTRING"),]["lon_ini"],
#'      add = TRUE,
#'      lwd = 2,
#'      pal = colorRampPalette(c("black", "blue")))
#'}
ray_source <- function(betam,
                       u,
                       lat,
                       x0,
                       y0,
                       K,
                       cx,
                       dt,
                       itime,
                       direction,
                       interpolation = "trin",
                       tl = 1,
                       a = 6371000,
                       verbose = FALSE,
                       ofile){
  # combine the x0 and y0:
  # df <- expand.grid(x0, y0)
  # names(df) <- c("lon", "lat")
  df <- data.frame(lon = x0,
                   lat = y0)
  cat("\nInitial Poits at: \n")
  print(df)


  dir <- direction
  wn <- K

  ddir <- lapply(seq_along(dir), function(i){
    ddf <- lapply(1:nrow(df), function(j){
      dwn <- lapply(1:length(wn), function(k){
       dx <- cbind(ray(betam = betam,
                  u = u,
                  lat = lat,
                  itime = itime,
                  K = wn[k],
                  # cx = cx,
                  dt = dt,
                  direction = dir[i],
                  interpolation = interpolation,
                  x0 = df$lon[j],
                  y0 = df$lat[j],
                  verbose = verbose),
              id = paste0("K",wn[k],"_lati",df$lat[j],"_loni",df$lon[j]),
              direction = dir[i])

       cat(paste0("\nK = ",wn[k],"  lat = ", df$lat[j],"   lon = ",df$lon[j]),"\n")
      # Avoid stopping when the wave does not propagates, that is,
      # when "Error in geo[i + 1, ] : subscript out of bounds" happens
      # #ERROR HANDLING
      possibleError <- tryCatch(
      dl <- sf::st_as_sf(sf::st_set_geometry(dx, NULL),
                         geometry = ray_path(x = dx$lon,
                                              y = dx$lat)),
        error=function(e) e
      )

       if(!inherits(possibleError, "error")){
       #REAL WORK
         dl <- sf::st_as_sf(sf::st_set_geometry(dx, NULL),
                            geometry = ray_path(x = dx$lon,
                                                y = dx$lat))
         rbind(dx, dl)
       } else {
         message(possibleError)
         cat(paste0("\nSkipping to the next iteration... \n"))
       }
      })
        do.call("rbind", dwn)
    })
      do.call("rbind", ddf)
  })

    DF <- do.call("rbind", ddir)

    if(!missing(ofile)) {
      write.csv(x = DF, file = ofile)
    }
    return(DF)
}
