

CalculateMeanReddAreaHsi <- function(hsi.ras, hsi.spdf, search.radius, 
   col.to.summarize, new.col.name) {  

  require(raster)
  require(nabor)

  # get cell size
  cell.size_m <- res(hsi.ras)[1]

  # calculate number of points in the square that contains the circle
  #  of the radius of interest; increasing by 10% adds padding
  max.nn <- ceiling(((ceiling(2 * search.radius / cell.size_m)) ** 2) * 1.1)

  # find nearest neighbors to include in calculations
  #  umcomment "query = ..." to query a subset of points for testing
  #  or development; leave it commented to query the entire set against itself
  nabor.nn2 <- knn(
    data = hsi.spdf@coords,
    # query = hsi.pts[1:30000, c('x', 'y')],  # for querying subsets
    k = max.nn
  )

  # add mean redd area hsi as an attribute
  hsi.spdf[[new.col.name]] <- vapply(
    1:nrow(hsi.spdf),
    function(x) {
      inds <- nabor.nn2$nn.idx[x, ][nabor.nn2$nn.dists[x, ] <= search.radius]
      mean(hsi.spdf@data[inds, col.to.summarize])
    },
    FUN.VALUE = numeric(1)
  )

  # create raster with the calculated mean redd area hsi values 
  mean.redd.area.hsi.ras <- rasterize(
    x = hsi.spdf,
    y = hsi.ras,
    field = new.col.name
  )

  return(
    list(
      hsi.spdf = hsi.spdf,
      mean.redd.area.hsi.ras = mean.redd.area.hsi.ras
    )
  )

}


PredictReddLocations <- function(mean.ra.hsi.spdf, hsi.col.name, mean.ra.hsi.col.name,
  hsi.thresh, terr.area_m2) {

  require(dplyr)
  require(fields)

  # red capacity based on WUA and grey-lit-sourced defended redd territory area
  inds <- which(mean.ra.hsi.spdf[[hsi.col.name]] >= hsi.thresh)
  wua.pred.hsi.gte0.8 <- floor(
    sum(mean.ra.hsi.spdf[[hsi.col.name]][inds] * 0.1 * 0.1 ) / terr.area_m2
  )

  # use terr.area_m2 to calculate territory radius in m
  terr.rad_m <- sqrt(terr.area_m2 / pi)

  print(hsi.col.name)
  print(paste('WUA: ', as.character(sum(mean.ra.hsi.spdf[[hsi.col.name]][inds] * 0.1 * 0.1 ))))
  print(paste('upper capacity limit: ', as.character(sum(mean.ra.hsi.spdf[[hsi.col.name]][inds] * 0.1 * 0.1 / terr.area_m2)))) 
  
  #-----------------------------------------------------------------------------
  # work with locations above the hsi threshold
  #-----------------------------------------------------------------------------

  pdf.gte0.8 <- as.data.frame(mean.ra.hsi.spdf) %>%
    filter(fuz.spawn.hsi >= hsi.thresh) %>%
    dplyr::select(x, y, fuz.spawn.hsi, mean.redd.area.hsi) %>%
    arrange(desc(mean.redd.area.hsi)) %>%
    mutate(hsi.threshold = hsi.thresh)
  # head(pdf.gte0.8)

  qualifying.locs.exist <- (nrow(pdf.gte0.8) >= 1) & 
    (wua.pred.hsi.gte0.8 >= 1)

  if (qualifying.locs.exist) {

    # record first predicted location
    redd.loc.df.2cat <- pdf.gte0.8[1, ]

    cand.row <- 2
    
    while (nrow(redd.loc.df.2cat) < wua.pred.hsi.gte0.8 & cand.row <= nrow(pdf.gte0.8)) {
    
      # print(paste('working on', cand.row, 'of', nrow(pdf.gte0.8), sep = ' '))
      
      candidate.point <- pdf.gte0.8[cand.row, c('x', 'y')]
       
      dists <- rdist(redd.loc.df.2cat[, c('x', 'y')], candidate.point)
      
      if (all(dists > terr.rad_m)) {
      
        redd.loc.df.2cat <- rbind(redd.loc.df.2cat, pdf.gte0.8[cand.row, ])
      
      }
      
      cand.row <- cand.row + 1    
    
    }

    # check
    all(dist(redd.loc.df.2cat[, c('x', 'y')]) <= terr.rad_m)

  } else {

    print('no redds met the hsi threshold')
    redd.loc.df.2cat <- pdf.gte0.8[0, ]

  }

  redd.loc.df.2cat

}

