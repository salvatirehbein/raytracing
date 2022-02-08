NEWS
===========
# TODO

- Review descriptions in order to substitute "source" by "departure point" or
  other
- Implement RWS function
- Implement WFA function
- Depreciate function ray
- Include function for plotting
- Export uwnd in beta and betaks with lat/lon in ascending order

### raytracing 0.5.0 (Date: 08 feb 2022)
- Implement function for filter by area of interest (wave_arrival)
- ray and ray_source functions now exports in csv instead of RDS.

### raytracing 0.4.0 (Date: 06 feb 2022)
- Substitute "expand.grid" by "data.frame" in ray_source function
- fix bugs in documentation

### raytracing 0.3.0 (Date: 22 aug 2021)
- Create betam, Ktotal, and Ks function from betaks
- Export netcdf files to be read by GrADS 

### raytracing 0.2.0 (Date: 23 dez 2020)
- Add correction for longitudes bigger than 180
- Add DOI with zenodo in README

### raytracing 0.1.1 (Date: 17 dez 2020)
- Add tryCatch to the ray_path function in order to avoid stopping
  when the wave does not propagate: [geo(i + 1)] problem
- Add sep = ";" in write.csv function at ray and ray_source functions
- Fix typo error in documentation of betaks function, replacing lat by lon

### raytracing 0.1.0 (Date: 22 out 2020)
- ray_source returns great circles
- betaks now returns spatial polygons
- trigonometric interpolation was added
- Preparing submission to CRAN

### raytracing 0.0.0.9000 (Date: 16 sep 2020)
- First news
- Preparing submission to CRAN and JOSS
