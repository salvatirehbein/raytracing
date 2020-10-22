#' Calculates the Rossby waves ray paths
#'
#' \code{ray} returns the Rossby wave ray paths (lat/lon) triggered from
#' one initial source/position (x0, y0), one total wavenumber (K), and one
#' direction set up when invoking the function.
#' \code{ray} must ingest the meridional gradient of the absolute vorticity
#' in mercator coordinates\strong{betam}, the zonal mean wind \strong{u},
#' and the latitude vector (\strong{lat}). Those variables can be obtained
#' (recommended) using \code{\link{betaks}} function. The zonal means of the
#' basic state will be calculated along the \strong{ray} program, as well as
#' the conversion to mercator coordinates of \strong{u}.
#'
#' @param betam matrix (longitude = rows x latitude from minor to
#' major = columns) obtained with \code{\link{betaks}}. \strong{betam} is the
#' meridional gradient of the absolute vorticity in mercator coordinates
#' @param u matrix (longitude = rows x latitude from minor to
#' major = columns) obtained with \code{\link{betaks}}. Is the zonal wind speed
#' in the appropriate format for the \code{ray}. It will be converted in mercator
#' coordinates inside the \code{ray}
#' @param x0 Numeric value. Initial longitude (choose between -180 to 180)
#' @param y0 Numeric value. Initial latitude
#' @param lat Numeric vector of latitudes from minor to major
#'  (ex: -90 to 90). Obtained with \code{\link{betaks}}
#' @param K Numeric value; Total Rossby wavenumber
#' @param dt Numeric value; Timestep for integration (hours)
#' @param itime Numeric value; total integration time. For instance, 10 days
#' times 4 times per day
#' @param direction Numeric value (possibilities: 1 or -1)
#' It controls the wave displacement:
#' If 1, the wave goes to the north of the source;
#' If -1, the wave goes to the south of the source.
#' @param interpolation Character. Set the interpolation method to be used:
#' \code{\link{trin}} or \code{\link{ypos}}
#' @param tl Numeric value; Turning latitude. Do not change this!
#' It will always start with a positive tl (1) and automatically
#'  change to negative (-1) after the turning latitude
#' @param a Earth's radio (m)
#' @param ofile Character;  Output file name with .csv extension,
#' for instance, "/user/ray.csv"
#' @param verbose Boolean; if TRUE (default) return messages
#' during compilation
#' @seealso \code{\link{ray_source}}
#' @return sf data.frame
#' @importFrom utils write.csv
#' @importFrom sf st_as_sf
#' @export
#' @examples {
#' # For Coelho et al. (2015):
#' rm(list = ls()); graphics.off(); gc()
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                      package = "raytracing")
#' b <- betaks(u = input)
#' rt <- ray(betam = b$betam,
#'           u = b$u,
#'           lat = b$lat,
#'           K = 3,
#'           itime = 10 * 4,
#'           x0 = -130,
#'           y0 = -30,
#'           dt = 6,
#'           direction = -1,
#'           interpolation = "trin")
#' rp <- ray_path(rt$lon, rt$lat)
#' plot(rp,
#'      main = "Coelho et al. (2015): JFM/2015",
#'      axes = TRUE,
#'      cex = 2,
#'      graticule = TRUE)
#'}
ray <- function(betam,
                u,
                lat,
                x0,
                y0,
                K,
                dt,
                itime,
                direction,
                interpolation = "trin",
                tl = 1,
                a = 6371000,
                verbose = FALSE,
                ofile){

  dt <- dt * 60 * 60

  # must be from north to south
  if(!is.unsorted(lat)) { # is sorted from south to north
    message("Sorting latitude from north to south")
    betamz <- rev(colMeans(betam, na.rm = TRUE)) # in mercator coord.
    uz <- rev(colMeans(u, na.rm = TRUE))         # in latlon coord.
    lat <- lat[order(lat, decreasing = TRUE)]
  } else {
    betamz <- colMeans(betam, na.rm = TRUE)      # in mercator coord.
    uz <- colMeans(u, na.rm = TRUE)              # in latlon coord.
  }

  phirad <- lat*pi/180

  k <- K/a
  k2 <- (K/a)^2
  x_ini <- x0
  y_ini <- y0

  l_time <- list()
  l_x0 <- list()
  l_y0 <- list()


  count <- 0
  while(TRUE) {
    count <- count + 1

    # Break 1
    if(abs(y0) > 90) break

    if(interpolation == "trin") {

      beta_y0 <- trin(y = y0, yk = betamz)
      u_y0 <- trin(y = y0, yk = uz, mercator = TRUE)

    } else if(interpolation  == "ypos"){

      beta_y0 <- ypos(y = y0, lat = lat, yk = betamz)
      u_y0 <- ypos(y = y0, lat = lat, yk = uz, mercator = TRUE)

    }

    Ks2_y0 <- beta_y0/u_y0
    l2_y0 <- Ks2_y0 - k2

    ## Skip to the next timestep
    if(count == 1 & l2_y0 < 0) next

     if(round(l2_y0, 13) <= 0) tl <- -1

    # Break 2
    if(Ks2_y0 < 0) break

    ug0 <-(2 * beta_y0 * k2 / Ks2_y0^2 ) * 180 / (a*pi)

    vg0 <- (direction * tl *
              2 * beta_y0 * k * sqrt(abs(l2_y0))*
              cos(y0*pi/180)/Ks2_y0^2)*(180/(pi*a))

    y1 <- y0 + dt*vg0

    if(interpolation == "trin") {

      beta_y1 <- trin(y = y1, yk = betamz)
      u_y1 <- trin(y = y1, yk = uz, mercator = TRUE)

    } else if(interpolation  == "ypos"){

      beta_y1 <-  ypos(y = y1, lat = lat, yk = betamz)
      u_y1 <- ypos(y = y1, lat = lat, yk = uz, mercator = TRUE)

    }

    Ks2_y1 <- beta_y1/u_y1
    l2_y1 <- Ks2_y1 - k2

    # Skip to the next timestep
    if(count == 1 & l2_y1 < 0) next

    if(round(l2_y1, 13) < 0) tl <- -1

    # Break 3
    if(Ks2_y1  < 0) break

    ug1 <-  (2 * beta_y1 * k2 / Ks2_y1^2 ) * 180 / (a*pi)

    vg1 <- (direction * tl *
              2 * beta_y1 * k * sqrt(abs(l2_y1))*
              cos(y1*pi/180)/Ks2_y1^2)*(180/(pi*a))

    x2 <- x0 + 0.5*dt*(ug0+ug1)
    y2 <- y0 + 0.5*dt*(vg0+vg1)

    x0 <- x2
    y0 <- y2

    l_time[[count]] <- count
    l_x0[[count]] <- x0
    l_y0[[count]] <- y0

    # Break 4
    if(count == itime) break
  }

  count <- c(0, unlist(l_time))
  df <- data.frame(
    K = K,
    lat_ini = y_ini,
    lon_ini = x_ini,
    time = dt*count/(60*60),
    iday = rep(seq(0, 23, dt/(60*60)), length(count))[1:length(count)],
    lat = c(y_ini, unlist(l_y0)),
    lon = c(x_ini, unlist(l_x0))
    )

  df$lon_shift <- ifelse(df$lon < 0, df$lon + 360, df$lon)
  df$y0 <- df$lat
  df$x0 <- df$lon
  pontos <- sf::st_as_sf(df,
                         coords = c("x0", "y0"),
                         crs = 4326)

  if(!missing(ofile)) {
    utils::write.csv(x = pontos, file = ofile, row.names = FALSE)
  }
  return(pontos)
}
