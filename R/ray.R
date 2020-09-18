#' Calculate the Rossby waves ray paths
#'
#' This function calculates the Rossby waves ray paths based on
#' Hoskin and Ambrizzi (1993).
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
#' # For Magana & Ambrizzi (1998):
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_300hPa_1997-98DJF.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' a <- ray(betamz = colMeans(b$betam, na.rm = TRUE),
#'                 umz = colMeans(b$um, na.rm = TRUE),
#'                 lat = b$lat,
#'                 K = 3,
#'                 itime = 10,
#'                 x0 = -130 + 360,
#'                 y0 = -17,
#'                 dt = 12 * 60 * 60,
#'                 direction = -1)
#'
#' # For Coelho et al. (2015):
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' a <- ray(betamz = colMeans(b$betam, na.rm = TRUE),
#'                 umz = colMeans(b$um, na.rm = TRUE),
#'                 lat = b$lat,
#'                 K = 3,
#'                 itime = 30,
#'                 x0 = -135 + 360,
#'                 y0 = -30,
#'                 dt = 6 * 60 * 60,
#'                 direction = -1)
#'
#' # Simple plot:
#' library(ggplot2)
#' ww <- map_data('world', ylim=c(-60,60))
#'
#' ggplot() +  theme_bw() +
#'   geom_polygon(data = ww,
#'                aes(x = long, y = lat, group = group),
#'                alpha = 0.0, col = "grey") +
#'   geom_point(data = a[!(a$tun_y0 == -1 |
#'                           a$tun_y1 == -1 | a$id == 0), ],
#'              aes(x = x0 - 360, y = y0), size = 3) +
#'   labs(x = NULL,
#'        y = NULL,
#'        title = NULL)
#'}
ray <- function(betamz,
                       umz,
                       lat,
                       x0 = -130 + 360,
                       y0 = -17,
                       K = 3,
                       dt = 12 * 60 * 60,
                       itime = 10,
                       direction = -1,
                       tl = 1,
                       a = 6371000,
                       verbose = FALSE,
                       ofile){
  k2 <- (K/a)^2
  x_ini <- x0
  y_ini <- y0
  # define lists for storing data
  l_time   <- list()
  l_lon_x0 <- list()
  l_lat_y0 <- list()
  l_y0     <- list()
  l_tun_y0 <- list()
  l_l2_y0 <- list()
  l_Ks2_y0  <- list()
  l_betamz_y0 <- list()
  l_umz_y0 <- list()
  l_lon_x1 <- list()
  l_lat_y1 <- list()
  l_y1     <- list()
  l_tun_y1 <- list()
  l_l2_y1 <- list()
  l_Ks2_y1  <- list()
  l_betamz_y1 <- list()
  l_umz_y1 <- list()

  # Beginning of the while loop parent 1 ####
  count <- 0
  while(TRUE) {
    count <- count + 1

    if(verbose) {
      cat("iteration: ", count,
          ". poslat: ", ypos(y0, lat = lat), "\n  lat y0:", y0, "\n")
    }

    ug0 <- calcUg(betamz = betamz, umz = umz, y = y0, lat = lat)
    vg0 <- calcVg(betamz = betamz, umz = umz, y = y0, lat = lat,
                  K = K, direction = direction)

    # Write the elements at y0 on the lists
    Ks2_y0 <- betamz[ypos(y0, lat = lat)] / umz[ypos(y0, lat = lat)]
    l2_y0 <- Ks2_y0 - k2

    l_time[[count]] <- count
    l_y0[[count]] <- ypos(y0, lat = lat)
    l_lon_x0[[count]] <- x0
    l_lat_y0[[count]] <- y0
    l_tun_y0[[count]] <- ifelse(l2_y0 < 0, -1, 1)
    l_l2_y0[[count]] <- l2_y0
    l_Ks2_y0[[count]] <- Ks2_y0
    l_betamz_y0[[count]] <- betamz[ypos(y0, lat = lat)]
    l_umz_y0[[count]] <- umz[ypos(y0, lat = lat)]

    y1 <- y0 + dt*vg0

    # Write the elements at y1 on the lists
    Ks2_y1 <- betamz[ypos(y1, lat = lat)] / umz[ypos(y1, lat = lat)]
    l2_y1 <- Ks2_y1 - k2

    l_y1[[count]] <- ypos(y1, lat = lat)
    l_lat_y1[[count]] <- y1
    l_tun_y1[[count]] <- ifelse(l2_y1 < 0, -1, 1)
    l_l2_y1[[count]] <- l2_y1
    l_Ks2_y1[[count]] <- Ks2_y1
    l_betamz_y1[[count]] <- betamz[ypos(y1, lat = lat)]
    l_umz_y1[[count]] <- umz[ypos(y1, lat = lat)]

    # condition for leaving the parent while loop 1: turning latitude (l2 < 0)
    if(l2_y0 < 0){
      if(verbose) cat("l2_y0 < 0")
      break
    }

    ug1 <- calcUg(betamz = betamz, umz = umz, y = y1, lat = lat)
    vg1 <- calcVg(betamz = betamz, umz = umz, y = y1, lat = lat,
                  K = K, direction = direction)

    x2 <- x0 + 0.5*dt*(ug0+ug1)
    y2 <- y0 + 0.5*dt*(vg0+vg1)

    x0 <- x2
    y0 <- y2

    if(count == itime) {
      break
    }
  }
  message("\nLeaving parent while loop 1\n")

  # Store data before turning latitude
  dff <- data.frame(K = rep(K, length(unlist(l_time))),
                    x_ini = rep(x_ini, length(unlist(l_time))),
                    y_ini = rep(y_ini, length(unlist(l_time))),
                    direction = rep(direction, length(unlist(l_time))),
                    time = unlist(l_time),
                    x0 = unlist(l_lon_x0),
                    y0 = unlist(l_lat_y0),
                    y1 = unlist(l_lat_y1),
                    ilat_y0 = unlist(l_y0),
                    ilat_y1 = unlist(l_y1),
                    tun_y0 = unlist(l_tun_y0),
                    tun_y1 = unlist(l_tun_y1),
                    l2_y0 = unlist(l_l2_y0),
                    l2_y1 = unlist(l_l2_y1),
                    Ks2_y0 = unlist(l_Ks2_y0),
                    Ks2_y1 = unlist(l_Ks2_y1),
                    betamz_y0 = unlist(l_betamz_y0),
                    betamz_y1 = unlist(l_betamz_y1),
                    umz_y0 = unlist(l_umz_y0),
                    umz_y1 = unlist(l_umz_y1))

  dff$id <- 1:nrow(dff)

  x0 <- dff[max(dff$id),]$x0
  y0 <- dff[max(dff$id),]$y0
  message("Restarting from the turning longitude: x0 = ", x0)
  message("Restarting from the turning latitude: y0 = ", y0)

  countold <- count

  # Beginning of the while loop parent 2 ####
  count <- 0
  while(TRUE) {
    count <- count + 1
    if(verbose) {
      cat("\niteration: ", count,
          ". poslat: ", ypos(y0, lat = lat), "\n")
    }

    ug0 <- calcUg(betamz = betamz, umz = umz, y = y0, lat = lat)
    vg0 <- calcVg(betamz = betamz, umz = umz, y = y0, lat = lat,
                  K = K, direction = direction, tl = -1)

    # Write the elements at y0 on the lists
    Ks2_y0 <- betamz[ypos(y0, lat = lat)] / umz[ypos(y0, lat = lat)]
    l2_y0 <- Ks2_y0 - k2

    l_time[[count]] <- count
    l_y0[[count]] <- ypos(y0, lat = lat)
    l_lon_x0[[count]] <- x0
    l_lat_y0[[count]] <- y0
    l_tun_y0[[count]] <- ifelse(l2_y0 < 0, -1, 1)
    l_l2_y0[[count]] <- l2_y0
    l_Ks2_y0[[count]] <- Ks2_y0
    l_betamz_y0[[count]] <- betamz[ypos(y0, lat = lat)]
    l_umz_y0[[count]] <- umz[ypos(y0, lat = lat)]

    y1 <- y0 + dt*vg0

    # Write the elements at y1 on the lists
    Ks2_y1 <- betamz[ypos(y1, lat = lat)] / umz[ypos(y1, lat = lat)]
    l2_y1 <- Ks2_y1 - k2

    l_y1[[count]] <- ypos(y1, lat = lat)
    l_lat_y1[[count]] <- y1
    l_tun_y1[[count]] <- ifelse(l2_y1 < 0, -1, 1)
    l_l2_y1[[count]] <- l2_y1
    l_Ks2_y1[[count]] <- Ks2_y1
    l_betamz_y1[[count]] <- betamz[ypos(y1, lat = lat)]
    l_umz_y1[[count]] <- umz[ypos(y1, lat = lat)]

    ug1 <- calcUg(betamz = betamz, umz = umz, y = y1, lat = lat)
    vg1 <- calcVg(betamz = betamz, umz = umz, y = y1, lat = lat,
                  K = K, direction = direction, tl = -1)

    x2 <- x0 + 0.5*dt*(ug0+ug1)
    y2 <- y0 + 0.5*dt*(vg0+vg1)

    if(abs(x0-x2)<2e-1 & abs(y0-y2)<2e-1) {
      message("x & y differences almost 0")
      break
    }
    x0 <- x2
    y0 <- y2


    if(count == itime) {
      message("itime reached")
      break
    }
  }

  # Store data before turning latitude
  dff2 <- data.frame(K = rep(K, length(unlist(l_time))),
                     x_ini = rep(x_ini, length(unlist(l_time))),
                     y_ini = rep(y_ini, length(unlist(l_time))),
                     direction = rep(direction, length(unlist(l_time))),
                     time = unlist(l_time),
                     x0 = unlist(l_lon_x0),
                     y0 = unlist(l_lat_y0),
                     y1 = unlist(l_lat_y1),
                     ilat_y0 = unlist(l_y0),
                     ilat_y1 = unlist(l_y1),
                     tun_y0 = unlist(l_tun_y0),
                     tun_y1 = unlist(l_tun_y1),
                     l2_y0 = unlist(l_l2_y0),
                     l2_y1 = unlist(l_l2_y1),
                     Ks2_y0 = unlist(l_Ks2_y0),
                     Ks2_y1 = unlist(l_Ks2_y1),
                     betamz_y0 = unlist(l_betamz_y0),
                     betamz_y1 = unlist(l_betamz_y1),
                     umz_y0 = unlist(l_umz_y0),
                     umz_y1 = unlist(l_umz_y1))
  dff2$id <- 0:(nrow(dff2)-1)
  DF <- rbind(dff, dff2)
  if(!missing(ofile)) {
    utils::write.csv(x = DF, file = ofile, row.names = FALSE)
  }
  return(DF)
}
