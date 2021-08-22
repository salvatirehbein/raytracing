#' Calculates Total Wavenumber for Stationary Rossby Waves (Ks)
#'
#' \code{Ks} ingests the time-mean zonal wind (u) and calculates the Total
#' Wavenumber for Stationary Rossby waves (Ks) in mercator coordinates
#' (see: Hoskins and Ambrizzi, 1993). Stationary Rossby waves are found when
#' zonal wave number (k) is constant along the trajectory, which leads to wave
#' frequency (omega) zero.
#' In this code Ks is used to distinguish the total wavenumber for Stationary
#' Rossby Waves (Ks) from the total wavenumber for Rossby waves (K), and
#' zonal wave number (k).
#' \code{Ks} returns a list with Ks in mercator coordinates (ksm).
#'
#' @param u String indicating the input data filename. The file to be
#' passed consists in a netCDF file with only time-mean zonal wind at one
#' pressure level, latitude in ascending order (not a requisite), and longitude
#' from 0 to 360. It is required that the read dimensions express
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
#' @return list with one vector (lat) and 1 matrix (Ksm)
#' @importFrom ncdf4 nc_open ncvar_get nc_close ncdim_def ncvar_def
#'  nc_create ncatt_put
#' @importFrom graphics filled.contour
#' @export
#' @examples {
#' # u is NetCDF and lat and lon characters
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' Ks <- Ks(u = input, plots = TRUE)
#' Ks$ksm[] <- ifelse(Ks$ksm[] >= 16 |
#'                    Ks$ksm[] <= 0, NA, Ks$ksm[])
#' cores <- c("#ff0000","#ff5a00","#ff9a00","#ffce00","#f0ff00")
#' graphics::filled.contour(Ks$ksm[, -c(1:5, 69:73)] ,
#'                          col = rev(colorRampPalette(cores, bias = 0.5)(20)),
#'                          main = "Ks")
#' }
Ks <- function(u,
               lat = "lat",      # lon,
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
  # The matrix for `ks_mercator` is an R object with more
  # points on y (rows), i.e., this matrix have more rows than columns, or
  # longitudes are in y and latitudes on x.
  # The functions `filled.contour` and `image` read this `ks_mercator` transposed
  # matrix, adjusting them automatically on the plot, so we see longitude on
  # x and latitude on y.
  # To shift the spatial coordinates of the matrix we need to deal with the
  # rows of our matrix in R in the following way.

  ksm_sf <- rbind(ks_mercator[(nrow(ks_mercator)/2 + 1):nrow(ks_mercator), ],
                  ks_mercator[1:(nrow(ks_mercator)/2), ])

  sfpoly <- sf::st_sf(data.frame(ksm = as.vector(ksm_sf)[1:length(g)]),
                      geometry = g,
                      crs = 4326)

  if(missing(ofile)){
    return(list(lat = phi,
                ksm = ks_mercator,
                sfpoly = sfpoly))
  } else {
    cat("Writting the netcdf in ",ofile,"...\n")

    # Write the netcdf ####
    # definition of dimensions
    west_east <- ncdf4::ncdim_def("west_east",
                                  units = "",
                                  vals = 1:nlon)
    south_north <- ncdf4::ncdim_def("south_north",
                                    units = "",
                                    vals = 1:nlat)

    # definition of variables
    XLONG <- ncdf4::ncvar_def(name = "XLONG",
                              units = "",
                              dim = list(west_east,south_north),
                              prec = "float")
    XLAT <- ncdf4::ncvar_def(name = "XLAT" ,
                             units = "",
                             dim = list(west_east, south_north),
                             prec = "float")
    KSM <- ncdf4::ncvar_def(name = "ks_mercator" ,
                            units = "",
                            dim = list(west_east, south_north),
                            prec = "float")

    vars_file <- ncdf4::nc_create(filename = ofile,
                                  vars = c(list('XLAT' = XLAT,
                                                'XLONG' = XLONG,
                                                'ks_mercator' = KSM)),
                                  force_v4 = FALSE)
    cat(paste("The file has", vars_file$nvars,"variables\n"))
    cat(paste("The file has", vars_file$ndim,"dimensions\n"))

    # Global attribute to the file when varid = 0
    # otherwise write the variable name
    ncdf4::ncatt_put(vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "title",
                     attval = "Basic state for calculate ray tracing")
    ncdf4::ncatt_put(vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "Author",
                     attval = "")
    ncdf4::ncatt_put(vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "institution",
                     attval = "Climate Group of Studies (GrEC)/University of Sao Paulo (USP)")
    ncdf4::ncatt_put(vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "title",
                     attval = "Basic state for calculate ray tracing")
    ncdf4::ncatt_put(vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "history",
                     attval = paste0("Created on ", Sys.time()))
    ncdf4::ncatt_put(vars_file,
                     varid = 0, # 0 para o arquivo
                     attname = "references",
                     attval = "See: Hoskins and Ambrizzi (1993), Yang and Hoskins (1996), and Rehbein et al. (2020)")

    # values for the basic variables ####
    ncdf4::ncvar_put(nc = vars_file,
                     varid = "XLONG",
                     vals = matrix(data = lon,
                                   ncol = nlat,
                                   nrow = nlon, byrow = FALSE))
    ncdf4::ncatt_put(vars_file,
                     varid = "XLONG",
                     attname = "MemoryOrder",
                     attval = "XY")
    ncdf4::ncatt_put(vars_file,
                     varid = "XLONG",
                     attname = "description",
                     attval = "LONGITUDE, 0 to 360")
    ncdf4::ncatt_put(vars_file,
                     varid = "XLONG",
                     attname = "units",
                     attval = "degree WMO")
    ncdf4::ncatt_put(vars_file,
                     varid = "XLONG",
                     attname = "stagger",
                     attval = "")
    ncdf4::ncatt_put(vars_file,
                     varid = "XLONG",
                     attname = "FieldType",
                     attval = 104)
    ncdf4::ncvar_put(nc = vars_file,
                     varid = "XLAT",
                     vals = matrix(data = phi,
                                   ncol = nlat,
                                   nrow = nlon, byrow = TRUE))
    ncdf4::ncatt_put(vars_file,
                     varid = "XLAT",
                     attname = "MemoryOrder",
                     attval = "XY")
    ncdf4::ncatt_put(vars_file,
                     varid = "XLAT",
                     attname = "description",
                     attval = "LATITUDE, SOUTH IS NEGATIVE")
    ncdf4::ncatt_put(vars_file,
                     varid = "XLAT",
                     attname = "units",
                     attval = "degree north")
    ncdf4::ncatt_put(vars_file,
                     varid = "XLAT",
                     attname = "stagger",
                     attval = "")
    ncdf4::ncatt_put(vars_file,
                     varid = "XLAT",
                     attname = "FieldType",
                     attval = 104)
    # Ks mercator
    ncdf4::ncvar_put(nc = vars_file,
                     varid = "ks_mercator",
                     vals = ks_mercator)
    ncdf4::ncatt_put(vars_file,
                     varid = "ks_mercator",
                     attname = "MemoryOrder",
                     attval = "XYZ")
    ncdf4::ncatt_put(vars_file,
                     varid = "ks_mercator",
                     attname = "description",
                     attval = "Basic state")
    ncdf4::ncatt_put(vars_file,
                     varid = "ks_mercator",
                     attname = "units",
                     attval = "")
    ncdf4::ncatt_put(vars_file,
                     varid = "ks_mercator",
                     attname = "stagger",
                     attval = "betam")
    ncdf4::ncatt_put(vars_file,
                     varid = "ks_mercator",
                     attname = "FieldType",
                     attval = 104)

    ncdf4::nc_close(vars_file)

    # Returning a list with the calculated variables
    return(list(lat = phi,
                ksm = ks_mercator,
                sfpoly = sfpoly))

  }
}
