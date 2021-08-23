#' Calculates Beta and Ks
#'
#' \code{betaks} ingests the time-mean zonal wind (u), transform it in
#' mercator coordinates (um); calculates the meridional gradient of
#' the absolute vorticity (beta) in mercator coordinates (betam);
#' and, finally, calculates stationary wavenumber (Ks) in mercator coordinates
#' (ksm) (see: Hoskins and Ambrizzi, 1993). \code{betaks} returns the um, betam,
#' and lat, for being ingested in \code{\link{ray}} or
#' \code{\link{ray_source}}.
#'
#' @param u String indicating the input data filename. The file to be
#' passed consists in a netCDF file with only time-mean zonal wind at one
#' pressure level, latitude in ascending order (not a requisite), and longitude
#' from 0 to 360.
#' It is required that the read dimensions express
#' longitude (in rows) x latitude (in columns).
#' \strong{u} also can be a numerical matrix with time-mean zonal wind at one
#' pressure level, latitude in ascending order (not a requisite), and longitude
#' from 0 to 360. It is required that the read dimensions express longitude
#' (in rows) x latitude (in columns).
#' @param ofile String indicating the file name for store output data.
#' If missing, will not return a netCDF file
#' @param uname String indicating the variable name field
#' @param lat String indicating the name of the latitude field. If
#' \strong{u} is a matrix, \strong{lat} must be numeric.
#' @param lon String indicating the name of the longitude field.If
#' \strong{u} is a matrix, \strong{lon} must be numeric from 0 to 360.
#' @param a Numeric indicating the Earth's radio (m)
#' @param plots Logical, if TRUE will produce filled.countour plots
#' @param show.warnings Logical, if TRUE will warns about NaNs in sqrt(<0)
#' @return list with one vector (lat) and 3 matrices (um, betam, and ksm)
#' @importFrom ncdf4 nc_open ncvar_get nc_close ncdim_def ncvar_def
#'  nc_create ncatt_put
#' @importFrom graphics filled.contour
#' @export
#' @examples {
#' # u is NetCDF and lat and lon characters
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(u = input, plots = TRUE)
#' b$ksm[] <- ifelse(b$ksm[] >= 16 |
#'                   b$ksm[] <= 0, NA, b$ksm[])
#' cores <- c("#ff0000","#ff5a00","#ff9a00","#ffce00","#f0ff00")
#' graphics::filled.contour(b$ksm[, -c(1:5, 69:73)] ,
#'                          col = rev(colorRampPalette(cores, bias = 0.5)(20)),
#'                          main = "Ks")
#'
#' # u, lat and lon as numeric
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.bin",
#'                       package = "raytracing")
#' u <- readBin(input,
#'              what = numeric(),
#'              size = 4,
#'              n = 144*73*4)
#' lat <- seq(-90, 90, 2.5)
#' lon <- seq(-180, 180 - 1, 2.5)
#' u <- matrix(u,
#'             nrow = length(lon),
#'             ncol = length(lat))
#' graphics::filled.contour(u, main = "Zonal Wind Speed [m/s]")
#' b <- betaks(u, lat, lon)
#' b$ksm[] <- ifelse(b$ksm[] >= 16 |
#'                   b$ksm[] <= 0, NA, b$ksm[])
#' cores <- c("#ff0000","#ff5a00","#ff9a00","#ffce00","#f0ff00")
#' graphics::filled.contour(b$ksm[, -c(1:5, 69:73)] ,
#'                          col = rev(colorRampPalette(cores, bias = 0.5)(20)),
#'                          main = "Ks")
#' }
betaks <- function(u,
                   lat = "lat",    # lon,
                   lon = "lon",
                   uname = "uwnd",   # lat
                   ofile,
                   a = 6371000,
                   plots = FALSE,
                   show.warnings = FALSE
) {
  if(inherits(u, "character")){
    message("Detecting NetCDF")
    ncu      <- ncdf4::nc_open(filename = u)
    u        <- ncdf4::ncvar_get(nc = ncu, varid = uname)
    lat      <- ncdf4::ncvar_get(nc = ncu, varid = lat)
    lon      <- ncdf4::ncvar_get(nc = ncu, varid = lon)
    ncdf4::nc_close(ncu)

  }


  dx       <- abs(lon[2] - lon[1])
  dy       <- abs(lat[2] - lat[1])
  nlat     <- length(lat)
  nlon     <- length(lon)

  # Define parameters and constants
  if(is.unsorted(lat)) {
    message("Sorting latitude from south to north")
    phi <- lat[order(lat, decreasing = FALSE)]
    u <- u[, nlat:1]
  } else {
    phi <- lat # nocov
  }

  if(plots) graphics::filled.contour(u, main = "u")

  omega <- 2*pi/(24*60*60)
  dphi <- dy*pi/180
  dphia <- dphi*a
  phirad <- matrix(phi*(pi/180),
                   nrow = 1,
                   ncol = nlat)

  # Calculate beta terms
  # 1st term: df/dy ####
  dfdy <- matrix(NA,
                 nrow = 1,
                 ncol = nlat - 2)
  for(j in 2:(nlat - 1)) {
    dfdy[1, j-1] <- (
      2*omega*sin(phirad[1, j + 1]) - 2*omega*sin(phirad[1, j - 1]))/(2*dphia)
  }

  #  message("df/dy must be close to 0 at the poles and > 0 at the Equator\n")
  if(plots) graphics::plot(as.vector(dfdy),
                           main= "df/dy",
                           pch = 16)

  # 2nd term: d2u/dy2 ####
  m_phirad <- matrix(as.vector(phirad),
                     nrow = nlon,
                     ncol = nlat,
                     byrow = TRUE)
  uc <- u*cos(m_phirad)

  d2udy2 <- matrix(NA,
                   nrow = nlon,
                   ncol = nlat - 2)
  for(j in 2:(nlat - 1)) {
    for(i in 1:nlon) {
      d2udy2[i, j-1] <- (uc[i, j + 1] - 2*uc[i, j] + uc[i, j-1]) / (dphia^2)

    }
  }
  if(plots) graphics::filled.contour(d2udy2, main = expression("d2u/dy2"))

  # Calculate Beta ####
  m_dfdy <- matrix(as.vector(dfdy),
                   nrow = nlon,
                   ncol = nlat - 2,
                   byrow = TRUE)
  beta <- m_dfdy - d2udy2

  beta_f <- cbind(beta[, 1],
                  beta,
                  beta[, ncol(beta)])

  if(plots) graphics::filled.contour(beta_f[,], main = "Beta")

  # Calculate Beta Mercartor --> beta * cos(phi) #####
  beta_mercator <- beta_f*cos(m_phirad)

  if(plots) graphics::filled.contour(beta_mercator,
                                     main = "Beta Mercator")

  # Ks mercator #######
  ks_mercator <- matrix(NA, nrow = nlon, ncol = nlat)

  ks_mercator[] <- ifelse(
    beta_mercator[] < 0 & u[] != 0, -1,
    ifelse(
      beta_mercator[] > 0 & u[] < 0, 30,
      ifelse(
        beta_mercator[] == 0 | u[] == 0, 0,
        if(show.warnings) {
          a * sqrt( (beta_mercator * cos(m_phirad)) / u)
        } else {
          suppressWarnings(a * sqrt( (beta_mercator*cos(m_phirad)) / u))
        }
      )))

  ks_mercator[] <- ifelse(ks_mercator[] >= 16 &
                            ks_mercator[] != 30, 20, ks_mercator[])

  if(plots) graphics::hist(ks_mercator, main = "Ks")
  if(plots) graphics::filled.contour(ks_mercator, main = "Ks")

  # Adding sf objects ####
  # grid
  bb <- c(xmin = min(lon - 180),
          ymin = min(lat),
          xmax = max(lon - 180),
          ymax = max(lat))

  bb <- sf::st_bbox(bb)

  g <- sf::st_make_grid(sf::st_as_sfc(bb), n = c(nlon, nlat))

  # Coordinates in x (longitude) are organized from 0 to 360
  # To transform to a spatial object with lat/long coordinates
  # we need to shift the x coordinates to -180 to 180.
  # The matrix for `u`, `beta_mercator`, and `ks_mercator` are R objects with more
  # points on y (rows), i.e., those matrix have more rows than columns, or
  # longitudes are in y and latitudes on x.
  # The functions `filled.contour` and `image` read those transposed matrix,
  # adjusting them automatically on the plot, so we see longitude on x and
  # latitude on y.
  # To shift the spatial coordinates of the matrix we need to deal with the
  # rows of our matrix in R in the following way.

  u_sf <- rbind(u[(nrow(u)/2 + 1):nrow(u), ],
                u[1:(nrow(u)/2), ])

  betam_sf <- rbind(beta_mercator[(nrow(beta_mercator)/2 + 1):nrow(beta_mercator), ],
                    beta_mercator[1:(nrow(beta_mercator)/2), ])

  ksm_sf <- rbind(ks_mercator[(nrow(ks_mercator)/2 + 1):nrow(ks_mercator), ],
                  ks_mercator[1:(nrow(ks_mercator)/2), ])

  sfpoly <- sf::st_sf(data.frame(u = as.vector(u_sf)[1:length(g)],
                                 betam = as.vector(betam_sf)[1:length(g)],
                                 ksm = as.vector(ksm_sf)[1:length(g)]),
              geometry = g,
              crs = 4326)

  if(missing(ofile)){
    return(list(lat = phi,
                u = u,
                betam = beta_mercator,
                ksm = ks_mercator,
                sfpoly = sfpoly))
  } else {
    cat("Writting the netcdf in ",ofile,"...\n")

    # Write the netcdf ####
    # definition of dimensions
    XLONG <- ncdf4::ncdim_def(name = "longitude",
                              units = "degrees_east",
                              vals = lon)
    XLAT <- ncdf4::ncdim_def(name = "latitude",
                             units = "degrees_north",
                             vals = -lat)

    # definition of variables
    UWND <- ncdf4::ncvar_def(name = "uwnd",
                             units = "m s^-1",
                             dim = list(XLONG, XLAT),
                             longname="Zonal wind")
    BETAM <- ncdf4::ncvar_def(name = "betam",
                              units = "s^-1 m^-1",
                              dim = list(XLONG, XLAT),
                              longname="Meridional Gradient of the Absolute Vorticity in mercator coordinates")
    KSM <- ncdf4::ncvar_def(name = "Ks",
                            units = "",
                            dim = list(XLONG, XLAT),
                            longname = "Stationary Wavenumber (Ks)")
    vars_file <- ncdf4::nc_create(filename = ofile,
                                  vars = list(UWND, BETAM, KSM))

    cat(paste("The file has", vars_file$nvars, "variables\n"))
    cat(paste("The file has", vars_file$ndim, "dimensions\n"))

    # Global attribute to the file when varid = 0
    # otherwise write the variable name
    ncdf4::ncatt_put(nc = vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "title",
                     attval = "Basic state for calculate ray tracing")
    ncdf4::ncatt_put(nc = vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "Author",
                     attval = "")
    ncdf4::ncatt_put(nc = vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "institution",
                     attval = "Climate Group of Studies (GrEC)/University of Sao Paulo (USP)")
    ncdf4::ncatt_put(nc = vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "title",
                     attval = "Basic state for calculate ray tracing")
    ncdf4::ncatt_put(nc = vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "history",
                     attval = paste0("Created on ", Sys.time()))
    ncdf4::ncatt_put(nc = vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "references",
                     attval = "See: Hoskins and Ambrizzi (1993), Yang and Hoskins (1996), and Rehbein et al. (2020)")

    # Add variables to the ofile
    ncdf4::ncvar_put(nc = vars_file,
                     varid = UWND,
                     vals = u,
                     start = c(1,1),
                     count = c(nlon,nlat))
    ncdf4::ncvar_put(nc = vars_file,
                     varid = BETAM,
                     vals = beta_mercator,
                     start = c(1,1),
                     count = c(nlon,nlat))
    ncdf4::ncvar_put(nc = vars_file,
                     varid = KSM,
                     vals = ks_mercator,
                     start = c(1,1),
                     count = c(nlon,nlat))

    ncdf4::nc_close(nc = vars_file)

    # Returning a list with the calculated variables
    return(list(lat = phi,
                u = u,
                betam = beta_mercator,
                ksm = ks_mercator,
                sfpoly = sfpoly))

  }
}
