# load required packages
library(rgdal)
library(rgeos)
library(raster)
library(tidyverse)
library(gridExtra)
library(grid)

# remove any objects in workspace
rm(list = ls())

# read in data
visits = read.csv("C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FishData/ReddData/UGR/UGR_PotentialSitesWithRedds_Subsample.csv", header = TRUE, stringsAsFactors = FALSE)
bfw = read.csv("C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/AllBasins/BankfullWidth.csv", header = TRUE, stringsAsFactors = FALSE)
redds = "C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FishData/ReddData/UGR/UGR_ChinookRedds_20142015_ChampSites_Snapped.shp"

# set output paths
out.reddCount = "C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FishData/ReddData/UGR/UGR_ChampSites_ReddCount.csv"
out.runVisits = "C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FishData/ReddData/UGR/UGR_ValidationSites.csv"

redds.pt = readOGR(dirname(redds), strsplit(basename(redds), '[.]')[[1]][1])

visitRedd.fn = function(x){

  redds.pt@data = mutate(redds.pt@data, redd.no = 1)
  visit.id = strsplit(strsplit(x, '/')[[1]][13], '_')[[1]][2]
  site.name = strsplit(x, '/')[[1]][12]
  visit.year = strsplit(x, '/')[[1]][10]
  basin = strsplit(x, '/')[[1]][11]
  
  we.path = file.path(x, 'Sims/FIS/Inputs') 
  
  we.poly = readOGR(we.path, layer = 'WaterExtent')
  
  visit.redds = as.data.frame(redds.pt[we.poly,]) %>%
    summarise(n.redds = sum(redd.no)) %>%
    mutate(WatershedName = basin, SiteName = site.name, VisitYear = visit.year, VisitID = visit.id)
  
  return(visit.redds)
}


tot.redds = visits %>% rowwise() %>% do(visitRedd.fn(.$visit.dir)) %>% bind_rows

write.csv(tot.redds, out.reddCount, row.names = FALSE)

visits.wredds = tot.redds %>% filter(n.redds > 0) %>% select(VisitID, n.redds) %>% mutate(VisitID = as.integer(VisitID))

bfw = bfw %>% select(VisitID, AverageBFWidth) %>% rename(AveBFW = AverageBFWidth)

val.sites = visits %>% filter(VisitID %in% visits.wredds$VisitID) %>% left_join(visits.wredds, by = 'VisitID') %>% left_join(bfw, by = 'VisitID') %>% select(-visit.url)

write.csv(val.sites, out.runVisits, row.names = FALSE)