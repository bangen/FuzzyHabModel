library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)
library(reshape2)

# remove any objects in workspace
rm(list = ls())

# # read in data
# visits = read.csv('C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FISValidation/ChinookSpawner/UGR_ValidationSites.csv', header = TRUE, stringsAsFactors = FALSE)
# redds = 'C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FishData/ReddData/UGR/UGR_ChinookRedds_20142015_ChampSites_Snapped.shp'
# 
# redd.pts = readOGR(dirname(redds), strsplit(basename(redds), '[.]')[[1]][1])
# 
# # redd point 
# reddVal.fn = function(x){
#   
#   visit.id = strsplit(strsplit(x, '/')[[1]][13], '_')[[1]][2]
#   site.name = strsplit(x, '/')[[1]][12]
#   visit.year = strsplit(x, '/')[[1]][10]
#   basin = strsplit(x, '/')[[1]][11]
#   
#   rpath = file.path(x, 'Sims/FIS/Output')
#   
#   fis.r = raster(file.path(rpath, 'FuzzyChinookSpawner_DVSC.tif'))
#   fis.df = as.data.frame(fis.r)
#   names(fis.df) = 'fis.val'
#   
#   #quantile(fis.df$fis.val, 0.08, na.rm = T)
#   #ecdf(fis.df$fis.val)(0.75)
#   
#   redd.fis.mean = extract(fis.r, redd.pts, buffer = 1.5, fun = mean, na.rm = TRUE, df = TRUE)
#   colnames(redd.fis.mean) = c('redd.id', 'redd.fisval') 
#   redd.fis.mean = redd.fis.mean %>% filter(!is.na(redd.fisval)) %>% dplyr::mutate(buffer.diam = 3.0, extract.stat = 'mean', fis.quantile = round(ecdf(fis.df$fis.val)(redd.fisval), 2))
#   
#   redd.fis.max = extract(fis.r, redd.pts, buffer = 1.5, fun = max, na.rm = TRUE, df = TRUE)
#   colnames(redd.fis.max) = c('redd.id', 'redd.fisval') 
#   redd.fis.max = redd.fis.max %>% filter(!is.na(redd.fisval)) %>% dplyr::mutate(buffer.diam = 3.0, extract.stat = 'max', fis.quantile = round(ecdf(fis.df$fis.val)(redd.fisval), 2))
#   
#   redd.results = rbind(redd.fis.max,redd.fis.mean) %>%
#     mutate(WatershedName = basin, SiteName = site.name, VisitYear = visit.year, VisitID = visit.id)
#   
#   return(redd.results)
# }
# 
# 
# redd.vals = visits %>% rowwise() %>% do(reddVal.fn(.$visit.dir)) %>% bind_rows
# write.csv(redd.vals, 'C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FISValidation/ChinookSpawner/UGR_Validation_ReddFISValues.csv', row.names = FALSE)

redd.vals = read.csv('C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FISValidation/ChinookSpawner/YankeeFork_Validation_ReddFISValues.csv', header = TRUE, stringsAsFactors = FALSE)
redd.vals %>% dplyr::group_by(extract.stat) %>% dplyr::summarise(ave.fis = mean(redd.fisval))
redd.vals %>% dplyr::group_by(extract.stat) %>% dplyr::summarise(ave.fis = mean(fis.quantile))

# Contigency Tables
quantile.tbl.mean = redd.vals %>% filter(extract.stat == 'mean') %>% mutate(quantile.group = ifelse(fis.quantile <= 0.5, '< q50', '> q50')) %>% group_by(quantile.group) %>% summarise(n = n())
chisq.test(as.matrix(quantile.tbl.mean[,-1]))

quantile.tbl.max = redd.vals %>% filter(extract.stat == 'max') %>% mutate(quantile.group = ifelse(fis.quantile <= 0.5, '< q50', '> q50')) %>% group_by(quantile.group) %>% summarise(n = n())
chisq.test(as.matrix(quantile.tbl.max[,-1]))

# Plots

#p1 = ggplot(data = redd.vals, aes(x = redd.fisval, y = fis.quantile, color = extract.stat)) +
p1 = ggplot(data = redd.vals, aes(x = redd.fisval, y = fis.quantile)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ extract.stat, ncol = 1) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggtitle('Redd FIS Value x FIS Value Percentile') +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = 'bottom')

png(filename = 'C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FISValidation/ChinookSpawner/UGR_ReddValuexPercentile.png', width = 8, height = 6, units = 'in', res = 500)
print(p1)
dev.off()

redd.vals.m = redd.vals %>% dplyr::select(SiteName, extract.stat, redd.fisval, fis.quantile) %>% melt(id = c('SiteName', 'extract.stat')) %>% rowwise() %>% mutate(SiteName = strsplit(SiteName, '-')[[1]][2])

p2 = ggplot(data = redd.vals.m, aes(x = SiteName, y = value, color = extract.stat)) +
  geom_point(alpha = 0.6) +
  facet_grid(variable ~ extract.stat) +
  theme_minimal() +
  ggtitle('Redd FIS Value and Percentile x Statistic Extracted in Redd Buffer') + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 90, hjust = 1))

png(filename = 'C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FISValidation/ChinookSpawner/UGR_ReddValuePercentilexStatExtracted.png', width = 12, height = 6, units = 'in', res = 500)
print(p2)
dev.off()

# summarise by fis membership function group
redd.mf = redd.vals %>% mutate(redd.mf = ifelse(redd.fisval < 0.1, 'Poor',
                                                ifelse(redd.fisval >= 0.9, 'High',
                                                       ifelse(redd.fisval > 0.1 & redd.fisval <= 0.2, 'Poor-Low',
                                                              ifelse(redd.fisval > 0.2 & redd.fisval <= 0.4, 'Low',
                                                                     ifelse(redd.fisval > 0.4 & redd.fisval <= 0.5, 'Low-Moderate',
                                                                            ifelse(redd.fisval > 0.5 & redd.fisval >= 0.8, 'Moderate', 'Moderate-High')))))))

p3 = ggplot(data = redd.mf, aes(x = redd.mf, color = extract.stat)) +
  geom_bar(alpha = 0.6) +
  facet_grid(~ extract.stat) +
  theme_minimal() +
  ggtitle('Redd FIS Value x MF Group') + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 90, hjust = 1))

png(filename = 'C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FISValidation/ChinookSpawner/JohnDay_ReddValuePercentilexStatExtracted.png', width = 12, height = 6, units = 'in', res = 500)
print(p2)
dev.off()
