
# quick notes about how placement is calculated:
#   - get the upper capacity from WUA/territory size (note: threshold WUA calculation of cells with value >= 0.4)
#   - create a filtered suitability raster by calculating mean value over window sized to territory
#   - place redds at highest filtered values until upper capacity limit is reached

library(tidyverse)
library(raster)

rm(list = ls(all = TRUE))
options(digits = 10)

# set parent folder
pf = "C:/etal/Shared/Projects/USA/CHaMP/HabitatSuitability/wrk_Data"

# for multiplot function
source('C:/etal/Shared/Projects/USA/CHaMP/HabitatSuitability/wrk_Code/07_ReddPlacement/newNreiFunctionsDraft_01.R')
source('C:/etal/Shared/Projects/USA/CHaMP/HabitatSuitability/wrk_Code/07_ReddPlacement/reddPlaceAlg.R')


# read in visit csv
visits.df = read.csv("C:/etal/Shared/Projects/USA/CHaMP/HabitatSuitability/wrk_Data/00_Projectwide/RunLists/visitList_Placement_Rerun.csv") %>%
  mutate(visit.dir = file.path(pf, WatershedName, SiteName, VisitYear, VisitDir))

# create tibble of visit directories
visit.dirs = visits.df %>% 
  # filter(VisitID == 4883) %>%
  pull(visit.dir)


# redd and territory dimensions

# sthd
sthd.redd.area_m2 <- 4.8
sthd.terr.area_m2 <- 4 * sthd.redd.area_m2
sthd.redd.rad_m <- sqrt(sthd.redd.area_m2 / pi)
# sthd.terr.rad_m <- sqrt(sthd.terr.area_m2 / pi)

# chk
chk.redd.area_m2 <- 6.7
chk.terr.area_m2 <- 4 * chk.redd.area_m2
chk.redd.rad_m <- sqrt(chk.redd.area_m2 / pi)
# chk.terr.rad_m <- sqrt(chk.terr.area_m2 / pi)

# chk juvenile
# Grant and Kramer, 1990
juv.length.cm <- 10.0
juv.terr.area_m2 = (10 ** (2.61 * log10(juv.length.cm) - 2.83))
juv.terr.rad_m <- sqrt(juv.terr.area_m2 / pi)

# hsi threshold
hsi.threshold = 0.4



place_fish = function(in.ras, terr.radius, terr.area, out.folder.name, out.file.name){
  # if spawner out.folder.name = 'ReddPlacement' 
  # if juvnile out.folder.name = 'FishPlacement'
  # if spawner out.file.name = 'SteelheadSpawner_PredReddLocations'
  
  # read in the raster data
  hsi.ras <- raster(in.ras)
  # str(hsi.ras)
  
  # convert raster to spdf points
  hsi.spdf <- rasterToPoints(hsi.ras, spatial = TRUE)
  names(hsi.spdf@data) <- 'fuz.spawn.hsi'
  
  
  # library(ggplot2)
  # ggplot(data = as.data.frame(hsi.spdf), aes(x = x, y = y, color = fuz.spawn.hsi)) + 
  #   geom_point() + 
  #   coord_fixed(ratio = 1) +
  #   scale_color_distiller(type = 'div', 'palette' = 'Spectral', direction = 1) +
  #   labs(title = NULL, x = 'UTM X', y = 'UTM Y', color = 'Fuzzy HQ') +
  #   theme_bw()
  
  mean.ra.hsi.list <- CalculateMeanReddAreaHsi(
    hsi.ras = hsi.ras,
    hsi.spdf = hsi.spdf,
    # search.radius = sthd.redd.rad_m, 
    search.radius = terr.radius, 
    col.to.summarize = 'fuz.spawn.hsi',
    new.col.name = 'mean.redd.area.hsi')
  
  # names(mean.ra.hsi.list)
  
  # extract raster object from list
  mean.redd.area.hsi.ras <- mean.ra.hsi.list$mean.redd.area.hsi.ras
  
  # extract spdf object from list
  hsi.spdf2 <- mean.ra.hsi.list$hsi.spdf
  
  rm(mean.ra.hsi.list)
  
  # write mean hsi raster and hsi.spdf points to file
  output.dir <- file.path(dirname(in.ras), out.folder.name)
  if (!dir.exists(output.dir)) {dir.create(output.dir, recursive = TRUE)}
  out.ras.fp <- file.path(output.dir, 'MeanSuitability.tif')
  
  library(rgdal)
  writeRaster(mean.redd.area.hsi.ras, filename = out.ras.fp,
              format = 'GTiff', overwrite = TRUE)
  
  # library(readr)
  # out.pts.fp <- file.path(output.dir, 'sthdHsiSpdfWithMeanRaHsi.csv')
  # write_csv(cbind(hsi.spdf@coords, hsi.spdf@data), path = out.pts.fp)
  
  redd.locs <- PredictReddLocations(
    mean.ra.hsi.spdf = hsi.spdf2,
    hsi.col.name = 'fuz.spawn.hsi',
    mean.ra.hsi.col.name = 'mean.redd.area.hsi',
    hsi.thresh = hsi.threshold,
    terr.area_m2 = terr.area
  )
  
  head(redd.locs, 25)
  dim(redd.locs)
  
  # write redd locs to file
  if(nrow(redd.locs) > 0){
    print(file.path(output.dir, paste(out.file.name, '.csv', sep = '')))
    redd.loc.fp <- file.path(output.dir, paste(out.file.name, '.csv', sep = ''))
    write_csv(redd.locs, path = redd.loc.fp)}
  
  
  #-----------------------------------------------------------
  # plots
  #-----------------------------------------------------------
  
  # head(as.data.frame(hsi.spdf2))
  
  library(ggplot2)
  hsi <- ggplot(data = as.data.frame(hsi.spdf2), aes(x = x, y = y, color = fuz.spawn.hsi)) +
    geom_point() +
    coord_fixed(ratio = 1) +
    scale_color_distiller(type = 'div', 'palette' = 'Spectral', direction = 1) +
    labs(title = 'fuzzyHS', x = 'UTM X', y = 'UTM Y', color = 'Fuz_Spawn_DVSC') +
    theme_bw()
  
  rp.hsi <- ggplot() +
    geom_point(data = as.data.frame(hsi.spdf2), aes(x = x, y = y, color = fuz.spawn.hsi)) +
    coord_fixed(ratio = 1) +
    scale_color_distiller(type = 'div', 'palette' = 'Spectral', direction = 1) +
    geom_point(data = redd.locs, aes(x = x, y = y), color = 'black') +
    labs(
      title = 'fuzzyHS modeled locations',
      x = 'UTM X',
      y = 'UTM Y',
      color = 'Fuz_Spawn_DVSC'
    ) +
    theme_bw()
  
  mean.redd.area.hsi <- ggplot(
    data = as.data.frame(hsi.spdf2),
    aes(x = x, y = y, color = mean.redd.area.hsi)
  ) +
    geom_point() +
    coord_fixed(ratio = 1) +
    scale_color_distiller(type = 'div', 'palette' = 'Spectral', direction = 1) +
    labs(title = 'Territory area mean fuzzyHS', x = 'UTM X', y = 'UTM Y', color = 'Mean_RA_HSI') +
    theme_bw()
  
  rp.mean <- ggplot() +
    geom_point(
      data = as.data.frame(hsi.spdf2),
      aes(x = x, y = y, color = mean.redd.area.hsi)) +
    coord_fixed(ratio = 1) +
    scale_color_distiller(type = 'div', 'palette' = 'Spectral', direction = 1) +
    geom_point(data = redd.locs, aes(x = x, y = y), color = 'black') +
    labs(
      title = 'Territory area mean fuzzyHS modeled locations',
      x = 'UTM X',
      y = 'UTM Y',
      color = 'Mean_RA_HSI'
    ) +
    theme_bw()
  
  plot.fn <- file.path(output.dir, paste(out.file.name, '.png', sep = ''))
  png(filename = plot.fn, width = 20 * 300, height = 10 * 300, units = 'px', res = 300)
  CreateMultiplot(hsi, rp.hsi, mean.redd.area.hsi, rp.mean, cols = 2)
  dev.off()
  
}

#-------------------------------------------------------------------------------
# begin looping through visits
#-------------------------------------------------------------------------------
# vis.fold.no = 11  # vis.fold.no = vis.fold.no + 1
# note: stopped here

for(visit.dir in visit.dirs){

  print(visit.dir)

  st.sp.ras = unlist(list.files(path = visit.dir, pattern = 'FuzzySteelheadSpawner_DVSC.tif$', full.names = TRUE, recursive = TRUE, include.dirs = FALSE))
  st.sp.locs = unlist(list.files(path = visit.dir, pattern = 'SteelheadSpawner_PredReddLocations.csv$', full.names = TRUE, recursive = TRUE, include.dirs = FALSE))
  
  ch.sp.ras = unlist(list.files(path = visit.dir, pattern = 'FuzzyChinookSpawner_DVSC.tif$', full.names = TRUE, recursive = TRUE, include.dirs = FALSE))
  ch.sp.locs = unlist(list.files(path = visit.dir, pattern = 'ChinookSpawner_PredReddLocations.csv$', full.names = TRUE, recursive = TRUE, include.dirs = FALSE))
  
  ch.juv.ras = unlist(list.files(path = visit.dir, pattern = 'FuzzyChinookJuvenile_DVS.tif$', full.names = TRUE, recursive = TRUE, include.dirs = FALSE))
  ch.juv.locs = unlist(list.files(path = visit.dir, pattern = 'ChinookJuvenile_PredFishLocations.csv$', full.names = TRUE, recursive = TRUE, include.dirs = FALSE))

  if(all(length(st.sp.ras) > 0, length(st.sp.locs) == 0)){
    place_fish(in.ras = st.sp.ras, terr.radius = sthd.redd.rad_m, terr.area = sthd.terr.area_m2, out.folder.name = 'ReddPlacement', out.file.name = 'SteelheadSpawner_PredReddLocations')
  }

  if(all(length(ch.sp.ras) > 0, length(ch.sp.locs) == 0)){
    place_fish(in.ras = ch.sp.ras, chk.redd.rad_m, chk.terr.area_m2, 'ReddPlacement', 'ChinookSpawner_PredReddLocations')
  }
  
  if(all(length(ch.juv.ras) > 0, length(ch.juv.locs) == 0)){
    place_fish(ch.juv.ras, juv.terr.rad_m, juv.terr.area_m2, 'FishPlacement', 'ChinookJuvenile_PredFishLocations')
  }


}  # end of for (vis.fold.no
