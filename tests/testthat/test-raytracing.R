# betaks function test ####
input <- system.file("extdata",
                     "uwnd.mon.mean_200hPa_2014JFM.nc",
                     package = "raytracing")
# Test the plots
b <- betaks(u = input, plots = TRUE)

# ypos function test ####
test_that("ypos works", {
  expect_equal(ypos(y = -30,
                    yk = rev(colMeans(b$u)),
                    lat = seq(90, -90, -2.5)),
               16.842706048)
})

# Test the netcdf output
b <- suppressWarnings(betaks(u = input,
                             ofile = tempfile(),
                             show.warnings = TRUE))

# ray function test ####
input <- system.file("extdata",
                     "uwnd.mon.mean_200hPa_2014JFM.nc",
                     package = "raytracing")
b <- betaks(u = input)
test_that("ray stops", {
  expect_error(ray(betam = b$betam,
                        u = b$u,
                        lat = b$lat,
                        K = 3,
                        itime = 6,
                        x0 = -135 + 360,
                        y0 = -30,
                        dt = 6 * 60 * 60),
  ".?")
})

a <- ray(betam = b$betam,
                u = b$u,
                lat = b$lat,
                K = 3,
                itime = 6,
                x0 = -135,
                y0 = -30,
                dt = 6 * 60 * 60,
                direction = -1,
                ofile = tempfile())

test_that("ray works", {
  expect_equal(round(a$iday[1]), 0)
  expect_message(ray(betam = b$betam,
                            u = b$u,
                            lat = b$lat,
                            K = 3,
                            itime = 30,
                            x0 = -135 ,
                            y0 = -30,
                            dt = 6 * 60 * 60,
                            direction = -1,
                            verbose = TRUE),
                 ".?")
})

