#' Calculate the meridional gradient of the absolute vorticity (beta)
#' and the stationary Rossby wavenumber (Ks)
#'
#' This function calculates the meridional gradient of the absolute
#'  vorticity (beta), zonal wind (u), and stationary Rossby wave in
#'  mercator coordinates as in Hoskins & Ambrizzi (1993).
#'  The input data is zonal wind with latitude in ascending order.
#' @param ifile Input file name consisting in zonal wind at one pressure
#' level and one time
#' @param ofile Output file name
#' @param varname Zonal wind variable name at a specified level and time
#' @param latname Latitude variable name
#' @param lonname Longitude variable name
#' @param a Earth's radio (m)
#' @param plots Logical, if TRUE will produce filled.countour plots
#' @param show.warnings Logical, if TRUE will warns about NaNs in sqrt(<0)
#' @return The meridional group velocity
#' @importFrom ncdf4 nc_open ncvar_get nc_close ncdim_def ncvar_def
#'  nc_create ncatt_put
#' @importFrom graphics filled.contour
#' @export
#' @examples {
#' input <- system.file("extdata",
#'                      "uwnd.mon.mean_200hPa_2014JFM.nc",
#'                       package = "raytracing")
#' b <- betaks(ifile = input)
#' output <- paste0(tempfile(), ".nc")
#' b <- betaks(ifile = input, ofile = output)
#' }
betaks <- function(ifile,
                   varname = "uwnd",   # lat
                   latname = "lat",   # lon,
                   lonname = "lon",
                   ofile,
                   a = 6371000,
                   plots = FALSE,
                   show.warnings = FALSE
) {

  ncu <- ncdf4::nc_open(filename = ifile)
  uwnd <- ncdf4::ncvar_get(nc = ncu, varid = varname)
  lat <- ncdf4::ncvar_get(nc = ncu, varid = latname)
  lon <- ncdf4::ncvar_get(nc = ncu, varid = lonname)
  ncdf4::nc_close(ncu)

  dx <- abs(lon[2] - lon[1])
  dy <- abs(lat[2] - lat[1])
  nlat <- dim(lat)
  nlon <- dim(lon)

  # Define parameters and constants
  if(is.unsorted(lat)) {
    phi <- -lat
  } else {
    phi <- lat # nocov
  }

  omega <- 2*pi/(24*60*60)
  dphi <- dy*pi/180
  dphia <- dphi*a
  phirad <- phi*(pi/180)
  m_phirad <- matrix(rep(phirad, nlon),
                     ncol = nlat,
                     nrow = nlon,
                     byrow = TRUE)

  if(plots) graphics::filled.contour(m_phirad, main = "phirad")

  #  U Mercator ####
  if(is.unsorted(lat)) {
    u <- uwnd[, nlat:1]
  } else {
    u <- uwnd[,]  # nocov
  }
  u_mercator <- u / cos(m_phirad)

  if(plots) graphics::filled.contour(uwnd, main = "uwnd")
  if(plots) graphics::filled.contour(u, main = "u")
  if(plots) graphics::filled.contour(u_mercator[,-c(1:5,(nlat-5):nlat)],
                                     main = "u mercator")

  # Calculate beta terms
  # 1st term: df/dy ####
  fdfdy <- function(j){
    (2*omega*sin(phirad[j + 1]) - 2*omega*sin(phirad[j - 1]))/(2*dphia)
  }

  dfdy <- fdfdy(j = 2:(length(phi)-1))
  length(dfdy)

  m_dfdy <- matrix(rep(dfdy, nlon),
                   ncol = nlat - 2 ,
                   nrow = nlon,
                   byrow = TRUE)

  # unique(m_dfdy[1,] == m_dfdy[12,]) # TRUE

  #  message("df/dy must be close to 0 at the poles and > 0 at the Equator\n")
  if(plots) graphics::filled.contour(m_dfdy,main= "df/dy")

  # 2nd term: d2u/dy2 ####
  fd2u <- function(m, i, j){
    (m[i, j + 1] - 2*m[i, j] + m[i, j-1]) / (dphia^2)
  }
  d2udy2 <- fd2u(m = u,
                 i = 1:nlon,
                 j = 2:(nlat-1))
  if(plots) graphics::filled.contour(d2udy2, main = expression("d2u/dy2"))

  # Calculate Beta ####
  beta <- m_dfdy - d2udy2

  beta_f <- cbind(matrix(NA, ncol = 1, nrow = nlon),
                  beta,
                  matrix(NA, ncol = 1, nrow = nlon))

  if(plots) graphics::filled.contour(beta_f[,], main = "Beta")

  # Calculate Beta Mercartor --> beta * cos(phi) #####
  mercator <- function(m, m2, i, j){
    m[i,j]*cos(m2[i,j])
  }
  beta_mercator <- mercator(m = beta_f,
                            m2 = m_phirad,
                            i = 1:nlon,
                            j = 1:nlat)
  if(plots) graphics::filled.contour(beta_mercator,
                                     main = "Beta Mercator")

  # Ks mercartor #######
  if(show.warnings) {
    fks <- function(beta, phirad, u, i, j){
      a * sqrt(( beta * cos(phirad[i,j])) / u[i,j])
    }

  } else {
    fks <- function(beta, phirad, u, i, j){
      a * suppressWarnings(sqrt(( beta * cos(phirad[i,j]))) / u[i,j])
    }

  }
  ks_mercator <- fks(beta = beta_mercator,
                     phirad = m_phirad,
                     u = u,
                     i = 1:nlon,
                     j = 1:nlat)

  if(plots) graphics::filled.contour(ks_mercator, main = "Ks")

  if(missing(ofile)){
    return(list(lat = phi,
                um = u_mercator,
                betam = beta_mercator,
                ksm = ks_mercator))
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
    UM <- ncdf4::ncvar_def(name = "u_mercator" ,
                           units = "",
                           dim = list(west_east, south_north),
                           prec = "float")
    BETAM <- ncdf4::ncvar_def(name = "beta_mercator" ,
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
                                                'u_mercator' = UM,
                                                'beta_mercator' = BETAM,
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
                     attval = "Hoskins and Ambrizzi (1993)")

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
    # u_mercator
    ncdf4::ncvar_put(nc = vars_file,
                     varid = "u_mercator",
                     vals = u_mercator)
    ncdf4::ncatt_put(vars_file,
                     varid = "u_mercator",
                     attname = "MemoryOrder",
                     attval = "XYZ")
    ncdf4::ncatt_put(vars_file,
                     varid = "u_mercator",
                     attname = "description",
                     attval = "Basic state")
    ncdf4::ncatt_put(vars_file,
                     varid = "u_mercator",
                     attname = "units",
                     attval = "m/s")
    ncdf4::ncatt_put(vars_file,
                     varid = "u_mercator",
                     attname = "stagger",
                     attval = "um")
    ncdf4::ncatt_put(vars_file,
                     varid = "u_mercator",
                     attname = "FieldType",
                     attval = 104)
    # beta_mercator
    ncdf4::ncvar_put(nc = vars_file,
                     varid = "beta_mercator",
                     vals = beta_mercator)
    ncdf4::ncatt_put(vars_file,
                     varid = "beta_mercator",
                     attname = "MemoryOrder",
                     attval = "XYZ")
    ncdf4::ncatt_put(vars_file,
                     varid = "beta_mercator",
                     attname = "description",
                     attval = "Basic state")
    ncdf4::ncatt_put(vars_file,
                     varid = "beta_mercator",
                     attname = "units",
                     attval = "s^-1 m^-1")
    ncdf4::ncatt_put(vars_file,
                     varid = "beta_mercator",
                     attname = "stagger",
                     attval = "betam")
    ncdf4::ncatt_put(vars_file,
                     varid = "beta_mercator",
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
                um = u_mercator,
                betam = beta_mercator,
                ksm = ks_mercator))

  }
}
