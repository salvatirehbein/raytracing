% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betam.R
\name{betam}
\alias{betam}
\title{Calculates Meridional Gradient of the Absolute Vorticity (beta)
in mercator coordinates (betam)}
\usage{
betam(
  u,
  lat = "lat",
  lon = "lon",
  uname = "uwnd",
  ofile,
  a = 6371000,
  plots = FALSE,
  show.warnings = FALSE
)
}
\arguments{
\item{u}{String indicating the input data filename. The file to be
passed consists in a netCDF file with only time-mean zonal wind at one
pressure level, latitude in ascending order (not a requisite), and longitude
from 0 to 360.
It is required that the read dimensions express
longitude (in rows) x latitude (in columns).
\strong{u} also can be a numerical matrix with time-mean zonal wind at one
pressure level, latitude in ascending order (not a requisite), and longitude
from 0 to 360. It is required that the read dimensions express longitude
(in rows) x latitude (in columns).}

\item{lat}{String indicating the name of the latitude field. If
\strong{u} is a matrix, \strong{lat} must be numeric.}

\item{lon}{String indicating the name of the longitude field.If
\strong{u} is a matrix, \strong{lon} must be numeric from 0 to 360.}

\item{uname}{String indicating the variable name field}

\item{ofile}{String indicating the file name for store output data.
If missing, it will not return a netCDF file}

\item{a}{Numeric indicating the Earth's radio (m)}

\item{plots}{Logical, if TRUE will produce filled.countour plots}

\item{show.warnings}{Logical, if TRUE will warns about NaNs in sqrt(<0)}
}
\value{
list with one vector (lat) and 2 matrices (u and betam)
}
\description{
\code{betam} ingests the time-mean zonal wind (u), transform it in
mercator coordinates (um) and then calculates the meridional gradient of
the absolute vorticity (beta) in mercator coordinates (betam) using
equation Karoly (1983). \code{betam} returns a list with the \strong{u},
\strong{betam}, and \strong{lat} for being ingested in \code{\link{Ktotal}},
\code{\link{Ks}}, \code{\link{ray}} or \code{\link{ray_source}}.
}
\examples{
{
# u is NetCDF and lat and lon characters
input <- system.file("extdata",
                     "uwnd.mon.mean_200hPa_2014JFM.nc",
                      package = "raytracing")
b <- betam(u = input, plots = TRUE)
cores <- c("#ff0000","#ff5a00","#ff9a00","#ffce00","#f0ff00")
graphics::filled.contour(b$betam/10e-12,
                         zlim = c(0, 11),
                         col = rev(colorRampPalette(cores)(24)),
                         main = "Beta Mercator (*10e-11)")
# u, lat and lon as numeric
input <- system.file("extdata",
                     "uwnd.mon.mean_200hPa_2014JFM.bin",
                      package = "raytracing")
u <- readBin(input,
             what = numeric(),
             size = 4,
             n = 144*73*4)
lat <- seq(-90, 90, 2.5)
lon <- seq(-180, 180 - 1, 2.5)
u <- matrix(u,
            nrow = length(lon),
            ncol = length(lat))
graphics::filled.contour(u, main = "Zonal Wind Speed [m/s]")
}
}
