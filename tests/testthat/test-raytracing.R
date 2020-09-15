# ypos function test ####
test_that("ypos works", {
  expect_equal(ypos(y = -30,
                    lat = seq(-90, 90, 2.5)),
               25)
})

# betaks function test ####
input <- system.file("extdata",
                     "uwnd.mon.mean_200hPa_2014JFM.nc",
                     package = "ray")
# Test the plots
b <- betaks(ifile = input, plots = TRUE)
# Test the if's

# Test the netcdf output
b <- betaks(ifile = input, ofile = tempfile(), show.warnings = TRUE)

# ray function test ####
input <- system.file("extdata",
                     "uwnd.mon.mean_200hPa_2014JFM.nc",
                     package = "raytracing")
b <- betaks(ifile = input)
a <- ray(betamz = colMeans(b$betam, na.rm = TRUE),
                umz = colMeans(b$um, na.rm = TRUE),
                lat = b$lat,
                K = 3,
                itime = 6,
                x0 = -135 + 360,
                y0 = -30,
                dt = 6 * 60 * 60,
                direction = -1,
                ofile = tempfile())

test_that("ray works", {
  expect_equal(round(a$umz_y1[1]), 22)
  expect_message(ray(betamz = colMeans(b$betam, na.rm = TRUE),
                            umz = colMeans(b$um, na.rm = TRUE),
                            lat = b$lat,
                            K = 3,
                            itime = 30,
                            x0 = -135 + 360,
                            y0 = -30,
                            dt = 6 * 60 * 60,
                            direction = -1,
                            verbose = TRUE),
                 "x & y differences almost 0")
})
