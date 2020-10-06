#' Calculates the Rossby waves ray paths
#'
#' \code{ray} must ingest the zonal means of \strong{betam} and \strong{um},
#' plus the latitude vector (\strong{lat}). \code{ray} returns the Rossby wave
#' ray paths (lat/lon) triggered from an initial source/position (x0, y0),
#' a total wavenumber (K), and direction set up when invoking the function.
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
#' @param x0 Numeric value; longitude (0 to 360 degrees)
#' @param y0 Numeric value; latitude
#' @param lat Numeric vector of latitudes from minor to major
#'  (ex: -90 to 90)
#' @param K Numeric value; Total Rossby wavenumber
#' @param dt Numeric value; Timestep for integration (seconds)
#' @param itime Numeric value; total integration time. For instance, 10 days
#' times 4 times per day.
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
#' @return sf data.frame
#' @importFrom utils write.csv
#' @importFrom sf st_sf
#' @export
#' @examples {
#' # For Magana and Ambrizzi (1998):
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_300hPa_1997-98DJF.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' ma98 <- ray(betamz = colMeans(b$betam, na.rm = TRUE),
#'                 umz = colMeans(b$um, na.rm = TRUE),
#'                 lat = b$lat,
#'                 K = 3,
#'                 itime = 10*4,
#'                 x0 = -130 + 360,
#'                 y0 = -17,
#'                 dt = 6 * 60 * 60,
#'                 direction = -1)
#' head(ma98)
#' # For Coelho et al. (2015):
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' co2015 <- ray(betamz = colMeans(b$betam, na.rm = TRUE),
#'                 umz = colMeans(b$um, na.rm = TRUE),
#'                 lat = b$lat,
#'                 K = 3,
#'                 itime = 10 * 4,
#'                 x0 = -130 + 360,
#'                 y0 = -30,
#'                 dt = 6 * 60 * 60,
#'                 direction = -1)
#' head(co2015)
#'}
ray <- function(betamz,
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
      cat(paste0("Iteration: ", count,
                 "   ypos: ", ypos(y0, lat = lat), " y0: ", round(y0, 2), "\n"))
    }


    Ks2_y0 <- betamz[ypos(y0, lat)]/umz[ypos(y0, lat)]
    ug0 <- calcUg(betamz = betamz, Ks2 = Ks2_y0, y = y0, lat = lat)
    vg0 <- calcVg(betamz = betamz, Ks2 = Ks2_y0, y = y0, lat = lat,
                  K = K, direction = direction)

    # Write the elements at y0 on the lists
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

    # Write the elements at y1 on the lists before test the turning latitude in y0
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
    # Turning latitude occurs when k = Ks --> l = 0
    if(round(l2_y0, 13) <= 0){
      message("l2_y0 <= 0\n",
              "Turning latitude reached \n",
              "Leaving parent while loop 1\n")
      break
    }

    ug1 <- calcUg(betamz = betamz, Ks2 = Ks2_y1, y = y1, lat = lat)
    vg1 <- calcVg(betamz = betamz, Ks2 = Ks2_y1, y = y1, lat = lat,
                  K = K, direction = direction)

    x2 <- x0 + 0.5*dt*(ug0+ug1)
    y2 <- y0 + 0.5*dt*(vg0+vg1)

    x0 <- x2
    y0 <- y2

    if(count == itime) {
      message("\nitime reached \n",
              "Leaving parent while loop 1\n")
      break
    }
  }


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
                    l2_y0_rounded = round(unlist(l_l2_y0), 13),
                    l2_y1_rounded = round(unlist(l_l2_y1), 13),
                    Ks2_y0 = unlist(l_Ks2_y0),
                    Ks2_y1 = unlist(l_Ks2_y1),
                    betamz_y0 = unlist(l_betamz_y0),
                    betamz_y1 = unlist(l_betamz_y1),
                    umz_y0 = unlist(l_umz_y0),
                    umz_y1 = unlist(l_umz_y1))


  # dff <- dff[dff$l2_y0_rounded >= 0.0e+00, ]
  dff$id <- 1:nrow(dff)

  x0 <- dff[max(dff$id),]$x0
  y0 <- dff[max(dff$id),]$y0
  message("Restarting from the turning longitude: x0 = ", round(x0, 3))
  message("Restarting from the turning latitude: y0 = ", round(y0, 3), "\n")


  # Beginning of the while loop parent 2 ####
  count <- 0
  while(TRUE) {
    count <- count + 1

    if(verbose) {
      cat(paste0("Iteration: ", count,
                 "   ypos: ", ypos(y0, lat = lat), " -->  lat y0:", round(y0, 2), "\n"))
    }


    # Stationary Rossby waves are only possible for beta >=0 & U > 0 (Hoskins & Ambrizzi, 1993):
    Ks2_y0 <- betamz[ypos(y0, lat)]/umz[ypos(y0, lat)]
    l2_y0 <- Ks2_y0 - k2

    ug0 <- calcUg(betamz = betamz, Ks2 = Ks2_y0, y = y0, lat = lat)
    vg0 <- calcVg(betamz = betamz, Ks2 = Ks2_y0, y = y0, lat = lat,
                  K = K, direction = direction, tl = -1)

    y1 <- y0 + dt*vg0
    Ks2_y1 <- betamz[ypos(y1, lat = lat)] / umz[ypos(y1, lat = lat)]
    l2_y1 <- Ks2_y1 - k2

    # Write the elements at y0 on the lists
    l_time[[count]] <- count
    l_y0[[count]] <- ypos(y0, lat = lat)
    l_lon_x0[[count]] <- x0
    l_lat_y0[[count]] <- y0
    l_tun_y0[[count]] <- ifelse(l2_y0 < 0, -1, 1)
    l_l2_y0[[count]] <- l2_y0
    l_Ks2_y0[[count]] <- Ks2_y0
    l_betamz_y0[[count]] <- betamz[ypos(y0, lat = lat)]
    l_umz_y0[[count]] <- umz[ypos(y0, lat = lat)]

    # Write the elements at y1 on the lists
    l_y1[[count]] <- ypos(y1, lat = lat)
    l_lat_y1[[count]] <- y1
    l_tun_y1[[count]] <- ifelse(l2_y1 < 0, -1, 1)
    l_l2_y1[[count]] <- l2_y1
    l_Ks2_y1[[count]] <- Ks2_y1
    l_betamz_y1[[count]] <- betamz[ypos(y1, lat = lat)]
    l_umz_y1[[count]] <- umz[ypos(y1, lat = lat)]


    ug1 <- calcUg(betamz = betamz, Ks2 = Ks2_y1, y = y1, lat = lat)
    vg1 <- calcVg(betamz = betamz, Ks2 = Ks2_y1, y = y1, lat = lat,
                  K = K, direction = direction, tl = -1)

    x2 <- x0 + 0.5*dt*(ug0+ug1)
    y2 <- y0 + 0.5*dt*(vg0+vg1)

    # Leave the while parent loop 2 when x & y differences are almost 0
    # The wave will stop propagating from this point
    if(abs(x0 - x2) < 1e-13 & abs(y0 - y2) < 1e-13) {
      message(paste0("|x0 - x2| & |y0 - y2| differences are almost 0 --> ",
                     "The wave is not propagating. \n",
                     "Leaving parent while loop 2\n"))
      break
    }
    x0 <- x2
    y0 <- y2


    if(count == itime) {
      message("itime reached \n",
              "Leaving parent while loop 2\n")
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
                     l2_y0_rounded = round(unlist(l_l2_y0), 13),
                     l2_y1_rounded = round(unlist(l_l2_y1), 13),
                     Ks2_y0 = unlist(l_Ks2_y0),
                     Ks2_y1 = unlist(l_Ks2_y1),
                     betamz_y0 = unlist(l_betamz_y0),
                     betamz_y1 = unlist(l_betamz_y1),
                     umz_y0 = unlist(l_umz_y0),
                     umz_y1 = unlist(l_umz_y1))

  # dff2 <- dff2[dff2$l2_y0_rounded >= 0.0e+00, ]

  dff2$id <- -(1:nrow(dff2))

  DF <- rbind(dff, dff2)

  # Filter turning latitudes & latlons where l2 ~ 0.
  DF <- DF[DF$tun_y0 != -1 &
             DF$tun_y1 != -1 &
             DF$l2_y0_rounded > 1e-13 &
             DF$l2_y1_rounded > 1e-13, ]

  DF$hora <- rep(1:2, nrow(DF))[1:nrow(DF)]
  DF <- DF[DF$hora == 1,]

  # Add points to the graph
  DF$x00 <- DF$x0 - 360
    pontos <- sf::st_as_sf(DF, coords = c("x00", "y0"), crs = 4326)

  # Calculate Great Circle for ray tracing
  r <- ray_path(DF$x0 - 360, DF$y0)
  DF <- sf::st_sf(DF, geometry = r)

  # Simple plot
  plot(DF$geometry, axes = TRUE)
  plot(pontos$geometry, add = TRUE, pch = 16, cex = 2)


  if(!missing(ofile)) {(DF)
    utils::write.csv(x = DF, file = ofile, row.names = FALSE)
  }
  return(DF)
}
