library(tidyverse)
library(data.table)
options(scipen = 999)

visits = read.csv("C:/et_al/Shared/Projects/USA/CHaMP/ResearchProjects/HabitatSuitability/wrk_Data/FISValidation/ChinookSpawner/UGR_ValidationSites.csv", header = TRUE, stringsAsFactors = FALSE)

hydro2inputs.fn = function(x){
  
  print(x)
  visit.id = as.integer(strsplit(strsplit(x, '/')[[1]][13], '_')[[1]][2])
  out.path = file.path(x, 'Sims/FIS/Inputs')

  hydro.path = list.files(file.path(x, 'Hydro'), pattern = 'dem_grid_results.csv', full.names = TRUE, recursive = FALSE)

  if(length(hydro.path) > 0){
    hydro.path = hydro.path[[1]][1]
    print(hydro.path)
    hydro = fread(input = hydro.path, nrows = -1, header = TRUE, data.table = FALSE, select = c('x', 'y', 'X', 'Y', 'Velocity.Magnitude', 'Depth'))
    hydro = if('X' %in% names(hydro)){rename(hydro, x = X, y = Y, Vel = Velocity.Magnitude)}else{rename(hydro, Vel = Velocity.Magnitude)}
    #hydro = mutate(hydro, Vel = replace(Vel, Vel == 0, 0.00001))
    hydro = mutate(hydro, Vel = round(as.numeric(Vel), 4), Depth = round(as.numeric(Depth), 4))
    hydro = mutate(hydro, Vel = replace(Vel, Vel == 0, 0.0001))
    write.csv(hydro, file.path(out.path, 'FuzzyHSI_Inputs.csv'), row.names = FALSE)
  }else{
    print(paste('WARNING: No hydo results for visit', visit.id))
  }
}

visits %>% rowwise() %>% do(hydro2inputs.fn(.$visit.dir))

