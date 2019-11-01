
#-------------------------------------------------------------------------------
# Written for Eco Logical Research, Inc. and South Fork Research, Inc.
# December 2014
#
# For information, questions, or suggestions, please contact:
# Eric Wall
# c.eric.wall@gmail.com
#
# Description:  A set of functions to implement a net rate of energy intake
#  model for drift-feeding salmonids. Functions in this script take 
#  rectilinear hydraulic output, create an on-the-fly 'wall of water' for each
#  raster cell normal to the flow direction of that cell and calculate
#  NREI.  The foraging model is Addley's (1993) derivation (adapted from
#  Hughes and Dill (1990)) which accounts for velocity shear near the foraging
#  point.  The method for calculating GREI is that of Hughes et al. (2003).
#  Swimming costs calculations are those of Hayes et al. (2007).
#
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Section I: Inputs
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Section II: Hydraulic Functions
#-------------------------------------------------------------------------------


#---------------------------------------------------------------------------
# Function to read and interpret rectilinear output from Delft3D and build
#  a 3D grid of the site's wetted points; parallel adaptation of MN's original
#  Build3DGrid function
#
# Args:
#  DEMfilename: name of the rectilinear Delft3D results file
#  DZ: spacing for 3D grid in Z direction (may need to be tighter than X 
#    and Y direction)
#  grd.reduction.factor: factor by which to increase X and Y grid spacing
#    relative to the 0.1m standard grid spacing; My be necessary for big
#    sites (as in the Entiat)
#  k.roughness: roughness height in m
#
# Returns:
#  a list object consisting of:
#    1. Grid3D, the 3D grid of all wetted points with their logarithmic
#      velocity values calculated
#    2. the Data.Use df, which is the original wetted df from Matt's output
#    3. a float number giving the spacing in X and Y of the resulting Grid3D
#      point cloud
#   
#---------------------------------------------------------------------------

Build3DGrid_par <- function(DEMfilename, DZ, grd.reduction.factor, k.roughness,
  sim.cores) {

  print("Reading Hydraulic Model Data File (This may take a minute!)")

  # Read the data.  Screen to only wetted areas.
  require(data.table)
  DEM.Results = fread(DEMfilename, header = TRUE, data.table = FALSE)
  
  # calculate the site area
  area.calc.temp <- DEM.Results[DEM.Results$Depth > 0,]
  # head(area.calc.temp)
  # dim(area.calc.temp)
  # summary(area.calc.temp)
  area.calc.temp <- area.calc.temp[area.calc.temp$BedLevel != -9999,]
  # plot(area.calc.temp$X, area.calc.temp$Y)
  # summary(area.calc.temp)
  my.area <- nrow(area.calc.temp) * 0.1 * 0.1
  my.volume <- sum(area.calc.temp$Depth * 0.1 * 0.1)
  rm(area.calc.temp)
  
  Xlev = levels(as.factor(DEM.Results$X))
  xidx = seq(1, length(Xlev), by = grd.reduction.factor)

  Ylev = levels(as.factor(DEM.Results$Y))
  yidx = seq(1, length(Ylev),by = grd.reduction.factor)

  DEM.Results = DEM.Results[(DEM.Results$X %in% Xlev[xidx]),]
  DEM.Results = DEM.Results[DEM.Results$Y %in% Ylev[yidx],]

  # dim(DEM.Results)
  Data.Use = DEM.Results[DEM.Results$Depth > 0,]
  Data.Use = Data.Use[Data.Use$BedLevel != -9999,]
  Data.Use$Z = rep(0, nrow(Data.Use))
  rm(DEM.Results)

  # # Check - Plot the velocity and depth
  # par(mfrow = c(2,1))

  # plot(Data.Use$X, Data.Use$Y, 
    # col = colorp[round(50*Data.Use$Depth/max(Data.Use$Depth))+1],main="Depth",
    # asp = 1
  # )

  # plot(Data.Use$X, Data.Use$Y, 
    # col = colorp[round(50*Data.Use$Velocity.Magnitude/max(Data.Use$Velocity.Magnitude))+1],main="Velocity",
    # asp = 1
  # )

  # Convert 2D Grid into 3D Grid.  We'll let z-grid spacing by tighter than X and Y grid spacing of .1 m

  # Set the maximum number of levels in Z. We'll build a rectangular prism encompassing all the wetted volume,
  # then trim out dry cells afterwards
  Zlevels = seq(trunc(min(Data.Use$BedLevel),1),round(max(c(Data.Use$WSE,Data.Use$BedLevel))),, by = DZ)
  #Zlevels = as.numeric(levels(factor(round(Data.Use$BedLevel,digits=1))))
  NZ = length(Zlevels)

  # register parallel backend
  require(foreach)
  require(doParallel)
  cl <- makeCluster(sim.cores)
  registerDoParallel(cl)

  # build 3D point cloud
  Grid3D <- foreach(nz = 1:NZ, .combine = 'rbind') %dopar% {
    
    # add appropriate Z values to Data.Use
    temp = Data.Use
    temp$Z = as.numeric(Zlevels[nz])
    
    # remove points outside of usable area
    temp = temp[temp$Z < temp$WSE,]
    temp = temp[temp$Z > temp$BedLevel,]

  }
  stopCluster(cl)
  
  # just to be sure
  Grid3D = Grid3D[Grid3D$Z < Grid3D$WSE,]
  Grid3D = Grid3D[Grid3D$Z > Grid3D$BedLevel,]

  order.z = order(Grid3D$Z)
  Grid3D = Grid3D[order.z,]
  order.y = order(Grid3D$Y)
  Grid3D = Grid3D[order.y,]
  order.x = order(Grid3D$X)
  Grid3D = Grid3D[order.x,]

  Grid3D$DA_X.Velocity = Grid3D$X.Velocity
  Grid3D$DA_Y.Velocity = Grid3D$Y.Velocity

  #######################################################################
  # logarithmic vertical velocity profile

  Grid3D$DA.Velocity.Magnitude = sqrt(
    (Grid3D$DA_X.Velocity * Grid3D$DA_X.Velocity) + (Grid3D$DA_Y.Velocity * Grid3D$DA_Y.Velocity)
  )

  #k.roughness= 10 # roughness

  vstar = Grid3D$DA.Velocity.Magnitude / (5.75* log(12.3 * (Grid3D$WSE-Grid3D$BedLevel) / k.roughness))
  Grid3D$Velocity.Magnitude =  5.75* log(30*(Grid3D$Z-Grid3D$BedLevel)/k.roughness)*vstar

  # Correct for conservation of flow
  XY = (factor(paste(Grid3D$X, Grid3D$Y)))
  Correction = tapply(Grid3D$Velocity.Magnitude, XY, sum)/(tapply(Grid3D$DA.Velocity.Magnitude,XY,sum)+.000001)


  Grid3D$Velocity.Magnitude = Grid3D$Velocity.Magnitude/(Correction[match(XY, levels(XY))]+.000000001)

  angle = atan2(Grid3D$Y.Velocity, Grid3D$X.Velocity)

  Grid3D$X.Velocity = Grid3D$Velocity.Magnitude *cos(angle)
  Grid3D$Y.Velocity = Grid3D$Velocity.Magnitude *sin(angle)


  Grid3D$Dir = atan2(Grid3D$Y.Velocity,  Grid3D$X.Velocity)+pi/2
  
  # add a colum, 'zrank', to indicate the relative position from the bed of 
  #  3d points for the various x-y combos
  temp.dt <- data.table(Grid3D[, c('X', 'Y', 'Z')])
  # names(temp.dt)
  temp.dt[, zrank := rank(Z), by = list(X, Y)]
  Grid3D$zrank <- temp.dt$zrank
  
  # min(Grid3D$Y.Velocity)
  # min(Grid3D$X.Velocity)
  # min(Grid3D$Velocity.Magnitude)
  return(
    list(
      "Grid3D"=Grid3D,
      "Data.Use"=Data.Use,
      "DX_DY" = .1*grd.reduction.factor,
      "my.area" = my.area,
      "my.volume" = my.volume
    )
  )

}


#---------------------------------------------------------------------------
# Function to initialize a radial grid we'll use over and over for every
#  XYZ point in our 3D grid
#  
# Args:
#  pts: number of angles in the radial grid (sorry about the bad name!)
#  max.Dist: maximum radius on any radial in the radial grid.
#    Note: this should be at least as far as fish reaction distance
#  Dist: distance between points along a radial in the radial grid.  Again,
#    lousy name!
#
# Returns:
#  a list object consisting of:
#    1. rad.grid.r: an array containing the radii lengths for each point in the 
#      radial grid; one row for each radial, one column for each point
#    2. rad.grid.theta: an array containing the thetas for each point in the 
#      radial grid; one row for each radial, one column for each point
#    3. Rmax: total number of points on each of the radials
#    4. VMult: a fixed matrix to help with mean normal velocity calculations
#   
#---------------------------------------------------------------------------

Initialize.Radial.Grid <- function(pts, max.Dist, Dist) {

  # Initialize some stuff

  angles = seq(-1 * pi, pi, by = (2 * pi) / pts)[1:pts]
  Dists = seq(Dist, max.Dist, by = Dist)

  rad.grid.r = (t(array(rep(c(Dists), pts), c(max.Dist/Dist,pts))))
  rad.grid.theta = array(angles[rep(rep((1:pts),max.Dist/Dist))], c(pts,max.Dist/Dist))

  # Initializing this once - outside the loops - for later... used to calculate mean
  # Velocities from NREI location every point in radial grid
  Rmax = max.Dist/Dist
  VMult = array(0, c(Rmax,Rmax+1))
  
  for (k in 1:Rmax) {
  
    VMult[k,] = c(.5, rep(1, (k-1)), .5,rep(0,(Rmax-k)))
  
  }
    VMult = t(VMult)

  
  #-------------------------------------------------------------
  # initialize an array to hold the 'slice' areas
  #-------------------------------------------------------------
  rad.grid.areas <- rad.grid.r * 0

  # first slices
  rad.grid.areas[, 1] <- pi * (1.5 * Dist) * (1.5 * Dist) / pts

  # second to second-to-last slices
  for (my.rad in 2:(Rmax - 1)) {
    
    rad.grid.areas[, my.rad] <- 2 * my.rad * pi * Dist * Dist / pts
  
  }

  # last slices
  rad.grid.areas[, Rmax] <- pi * Dist * (max.Dist - 0.25 * Dist) / pts

  # check on areas; threshold is 1 square mm
  if (sum(rad.grid.areas) - (pi * (max.Dist * max.Dist)) >= 0.000001) {
  
    print('!!!!!!!!  AREAS DID NOT ADD UP CORRECTLY  !!!!!!!!')
  
  }

  return(
    list(
     "rad.grid.r" = rad.grid.r,
     "rad.grid.theta" = rad.grid.theta,
     "Rmax" = Rmax,
     "VMult" = VMult,
     "rad.grid.areas" = rad.grid.areas
    )
  )
  
}


#---------------------------------------------------------------------------
# Helper function for mean normal velocity calculations; used in the
#  radial.grid function below.
#  
# Args:
#  rootVel: velocity at XYZ point being analyzed
#  Rmax: total number of points on each of the radials
#  VelMult: a matrix returned by the Initialize.Radial.Grid function above;
#    using this helps avoid for loops
#  VelVec: vector of velocities normal to flow along radius.
#
# Returns:
#   
#---------------------------------------------------------------------------

Mean.Zero.To.r <- function(rootVel, Rmax, VelVec, VelMult) {

  meanVel =array(c(rootVel,VelVec),c(1,Rmax+1)) %*% VelMult / (1:Rmax)

  return(meanVel)
  
}


#---------------------------------------------------------------------------
# Function to build a radial grid at the "wall of water" perpendicular to the
#  flow at an X-Y-Z grid point.
#  
# Args:
#  idx: index of XYZ grid; used to find location and velocity from the master
#    XYZ array
#  pts: number of angles in the radial grid
#  max.Dist: maximum radius on any radial in the radial grid
#    Note: this should be at least as far as fish reaction distance
#  Dist: distance between points along a radial in the radial grid
#  plots: boolean to turn QA plots on or off
#  rad.grid.r: an array containing the radii lengths for each point in the 
#    radial grid; one row for each radial, one column for each point;
#    of dimensions [pts, Rmax]
#  rad.grid.theta: an array containing the thetas for each point in the 
#    radial grid; one row for each radial, one column for each point;
#    of dimensions [pts, Rmax]
#  Rmax: number of possible points along radials (max.Dist/Dist)
#  VMult: See above function
#  DX_DY: grid spacing in X and Y in the 3D grid.
#
# Returns:
#  a list object containing many things.
#
#------------------
# Matt's notes:
#
# Function radial.grid
# Build a radial grid at the "wall of water" perpendicular to the flow at an X-Y-Z grid point.
# Return the grid, the normal velocities at each point in the grid, and the average
# velocity along each radius from the origin to each point at the radius.
#
# idx = index of X-Y-Z grid used to find location and velocity from master X-Y-Z array
# pts = number of angles
# max.Dist = maximum distance to look along radius (should be >= fish reaction distance)
# Dist = distance between points along radials
# plots: boolean to turn QA plots on or off.  if TRUE, plots will be generated at every 500th idx
# rad.grid.r: radial distances in radial grid of dimensions [pts, Rmax]
# rad.grid.theta: radial angles in radial grids of dimsions [pts, Rmax]
# Rmax: number of possible points along radials (max.Dist/Dist)
# VMult: See above function
# DX_DY: gird spacing in X and Y in the 3D grid.
#
#---------------------------------------------------------------------------

radial.grid <- function(idx, pts, max.Dist, Dist, plots,
  rad.grid.r, rad.grid.theta, Rmax, VMult,DX_DY) {

  require(RANN)

  #Data = Data.Use[
  #     Data.Use$X > (Grid3D$X[idx] - max.Dist) & 
  #     Data.Use$X < (Grid3D$X[idx] + max.Dist) & 
  #     Data.Use$Y > (Grid3D$Y[idx] - max.Dist) & 
  #     Data.Use$Y < (Grid3D$Y[idx] + max.Dist),]

  Data = Grid3D[
    Grid3D$X > (Grid3D$X[idx] - max.Dist) & 
    Grid3D$X < (Grid3D$X[idx] + max.Dist) & 
    Grid3D$Y > (Grid3D$Y[idx] - max.Dist) & 
    Grid3D$Y < (Grid3D$Y[idx] + max.Dist),]

  #rad.grid.r=as.vector(t(array(rep(c(Dists), pts), c(max.Dist/Dist,pts))))
  #rad.grid.theta=angles[rep(rep((1:pts),max.Dist/Dist))]

  # nrow(Grid3D)
  # Grid3D$X[idx]
  # Grid3D$Dir[idx]
  rad.grid.X =  Grid3D$X[idx] + cos(Grid3D$Dir[idx])*rad.grid.r*cos(rad.grid.theta)
  rad.grid.Y =  Grid3D$Y[idx] + sin(Grid3D$Dir[idx])*rad.grid.r*cos(rad.grid.theta)
  rad.grid.Z =  Grid3D$Z[idx] + rad.grid.r*sin(rad.grid.theta)

  #rad.grid.nn = nn2(Data[,1:2], data.frame(as.vector(rad.grid.X), as.vector(rad.grid.Y)),1)

  nn.col.idx = match( c("X","Y","Z"),colnames(Data))
  rad.grid.nn = nn2(Data[,nn.col.idx], 
    data.frame(as.vector(rad.grid.X), as.vector(rad.grid.Y),as.vector(rad.grid.Z)),1
  )

  # rad.grid.Z
  # rad.grid.X
  # rad.grid.Y
  rad.grid.nn.dists = array(rad.grid.nn$nn.dists, c(pts, max.Dist/Dist))

  rad.grid.Xvel = array(Data$X.Velocity[rad.grid.nn$nn.idx], c(pts, max.Dist/Dist))
  rad.grid.Yvel = array(Data$Y.Velocity[rad.grid.nn$nn.idx], c(pts, max.Dist/Dist))
  rad.grid.BedLevel = array(Data$BedLevel[rad.grid.nn$nn.idx], c(pts, max.Dist/Dist))
  rad.grid.WSE = array(Data$WSE[rad.grid.nn$nn.idx], c(pts, max.Dist/Dist))

  rad.grid.use = 
    (rad.grid.nn.dists < (DX_DY/1.2)) &
    (rad.grid.BedLevel < rad.grid.Z) &
    (rad.grid.WSE > rad.grid.Z)

  for (th in 2: (max.Dist/Dist)) {
  
    rad.grid.use[,th]= apply(1*rad.grid.use[,1:th], 1, prod)==1
  
  }

  rad.grid.Xvel[rad.grid.use==F] = 0
  rad.grid.Yvel[rad.grid.use==F] = 0
  rad.grid.Vmag = sqrt(rad.grid.Xvel^2 + rad.grid.Yvel^2)

  ##############################################################
  # Find the normal component of velocity for all grid points
  Vel.theta = atan2(rad.grid.Yvel, rad.grid.Xvel)
  Vel.Norm = rad.grid.Vmag * sin(-1*Vel.theta+Grid3D$Dir[idx])

  #################################################
  # For every radial grid point, find the mean velocity for all points along a radius
  # from the center to the given point.  This is a key NREI requirement.

  #initialize..
  rad.grid.mean.normV = rad.grid.X

  #Use the function defined outside the loop
  for (p in 1:pts) {
  
    rad.grid.mean.normV[p,]=Mean.Zero.To.r(Grid3D$Vmag[idx], Rmax, Vel.Norm[p,],VMult)
    
  }

  rad.grid.mean.normV = rad.grid.mean.normV * rad.grid.use

  ##################################################
  z.lim=c(min(Grid3D$Z), max(Grid3D$Z))

  if (plots) {
  
    if ((idx/plot.every) == round(idx/plot.every)) {

      par(mfrow = c(3,1))
      layout(matrix(c(1,2,3),3,1), heights=c(2,2,3.8))  # 
      par(mar=c(5,4,2,1))

      #par(mar=c(0,2,5,2))

      # Fun plots to watch as program runs.. but slows code down a LOT

      x.lim=c(min(rad.grid.X[rad.grid.use]),
        max(rad.grid.X[rad.grid.use])+
        .25*(max(rad.grid.X[rad.grid.use])-
        min(rad.grid.X[rad.grid.use]))
      )

      y.lim=c(min(rad.grid.Y[rad.grid.use]),
        max(rad.grid.Y[rad.grid.use])+
        .25*(max(rad.grid.Y[rad.grid.use])-
        min(rad.grid.Y[rad.grid.use]))
      )

      #plot(rad.grid.X[rad.grid.use], rad.grid.Z[rad.grid.use], col=
      #colorp[1+round(50*((rad.grid.Vmag[rad.grid.use]-min(rad.grid.Vmag[rad.grid.use]))/(max(rad.grid.Vmag[rad.grid.use])-
      #     min(rad.grid.Vmag[rad.grid.use]))))]#
      #, xlim=x.lim, xlab= "X-Coord (m)", ylab="Y-Coord (m)"
      #)
      #points(rad.grid.X, rad.grid.BedLevel,pch=19, cex=.5)
      #points(rad.grid.X, rad.grid.WSE,pch=19, cex=.5)
      #legend.text =as.character(round(seq(min(rad.grid.Vmag),   max(rad.grid.Vmag),by=(max(rad.grid.Vmag)-min(rad.grid.Vmag))/5),2))
      #legend.col = colorp[seq(1,51, by=10)]
      #legend("topright",legend.text, pch=19, col=legend.col, title=" Vel Mag(m/s)", bg="white")

      plot(rad.grid.X[rad.grid.use], rad.grid.Z[rad.grid.use], col=
        colorp[1+round(50*((rad.grid.mean.normV[rad.grid.use]-min(rad.grid.mean.normV[rad.grid.use]))/(max(rad.grid.mean.normV[rad.grid.use])-
        min(rad.grid.mean.normV[rad.grid.use]))))]#
        , xlim=x.lim, xlab= "X-Coord (m)", ylab="Z-Coord (m)", main = paste('idx = ', idx, sep = '')
      )
      points(rad.grid.X, rad.grid.BedLevel,pch=19, cex=.5)
      points(rad.grid.X, rad.grid.WSE,pch=19, cex=.5)

      legend.text =as.character(round(seq(min(rad.grid.mean.normV),   max(rad.grid.mean.normV),
        by=(max(rad.grid.mean.normV)-min(rad.grid.mean.normV))/5),2)
      )
      legend.col = colorp[seq(1,51, by=10)]
      legend("topright",legend.text, pch=19, col=legend.col, title=" Ave Velocity(m/s)", bg="white")

      #plot(rad.grid.Y[rad.grid.use], rad.grid.Z[rad.grid.use], col=
      #colorp[1+round(50*((rad.grid.Xvel[rad.grid.use]-min(rad.grid.Xvel[rad.grid.use]))/(max(rad.grid.Xvel[rad.grid.use])-
      #     min(rad.grid.Xvel[rad.grid.use]))))], xlim=y.lim,
      #     xlab= "Y-Coord (m)", ylab = "X-Coord (m)")
      #points(rad.grid.Y, rad.grid.BedLevel, pch=19, col="black", cex=.5)
      #points(rad.grid.Y, rad.grid.WSE, pch=19, cex=.5)
      #
      #legend.text =as.character(round(seq(min(rad.grid.Xvel),   max(rad.grid.Xvel),by=(max(rad.grid.Xvel)-min(rad.grid.Xvel))/5),2))
      #legend.col = colorp[seq(1,51, by=10)]
      #legend("topright",legend.text, pch=19, col=legend.col, title="X Velocity(m/s)")
      #

      plot(rad.grid.Y[rad.grid.use], rad.grid.Z[rad.grid.use], col=
        colorp[1+round(50*((rad.grid.mean.normV[rad.grid.use]-min(rad.grid.mean.normV[rad.grid.use]))/(max(rad.grid.mean.normV[rad.grid.use])-
        min(rad.grid.mean.normV[rad.grid.use]))))]#
        , xlim=y.lim, xlab= "Y-Coord (m)", ylab="Z-Coord (m)"
      )
      points(rad.grid.Y, rad.grid.BedLevel,pch=19, cex=.5)
      points(rad.grid.Y, rad.grid.WSE,pch=19, cex=.5)

      legend.text =as.character(round(seq(min(rad.grid.mean.normV),   max(rad.grid.mean.normV),
        by=(max(rad.grid.mean.normV)-min(rad.grid.mean.normV))/5),2)
      )
      legend.col = colorp[seq(1,51, by=10)]
      legend("topright",legend.text, pch=19, col=legend.col, title="Ave Velocity (m/s)", bg="white")

      plot(Grid3D$X, Grid3D$Y, col="gray",#, xlim=x.lim, ylim=y.lim,
        main=paste("X-Y-Z index", idx,": Radial Grid Location and Normal Velocity Component"),
        xlab="X-Coord (m)", ylab="Y-Coord(m)"
      )
      points(rad.grid.X[rad.grid.use], rad.grid.Y[rad.grid.use], col= 4*(rad.grid.use==T))

      scalar = 10

      for (i in 1: pts) {
      
        for (j in 1:Rmax) {
        
          if(rad.grid.use[i,j]) {
           
            lines(
              c(rad.grid.X[i,j], rad.grid.X[i,j] +scalar*rad.grid.Xvel[i,j]),
              c(rad.grid.Y[i, j], rad.grid.Y[i, j] + scalar*rad.grid.Yvel[i, j])
            )
            
            lines(
              c(rad.grid.X[i,j], rad.grid.X[i,j]+scalar*cos(Grid3D$Dir[idx]-pi/2)*Vel.Norm[i,j]),
              c(rad.grid.Y[i,j], rad.grid.Y[i,j]+scalar*sin(Grid3D$Dir[idx]-pi/2)*Vel.Norm[i,j]), col=4*(rad.grid.use==T)
            )
      
          }  # end if rad.grid.use
        }  # end of for j
      }  # end of for i
    }  # end of if ((idx/plot.every)...
  }  # end of if (plots)

  return(list(
    "rad.grid.origin.velocity" = Grid3D$Velocity.Magnitude[idx],
    "xy.water.column.depth" = Grid3D$WSE[idx] - Grid3D$BedLevel[idx],
    "rad.grid.origin.wse" = Grid3D$WSE[idx],
    "rad.grid.origin.Z" = Grid3D$Z[idx],
    "rad.grid.origin.bed.level" = Grid3D$BedLevel[idx],
    "rad.grid.X"= rad.grid.X, 
    "rad.grid.Y"=  rad.grid.Y,
    "rad.grid.WSE"= rad.grid.WSE, 
    "rad.grid.BedLevel" =rad.grid.BedLevel,
    "rad.grid.Xvel"=   rad.grid.Xvel,
    "rad.grid.Yvel"= rad.grid.Yvel,
    "rad.grid.Vel.Norm"= Vel.Norm,
    "rad.grid.mean.normV" = rad.grid.mean.normV, 
    "rad.grid.use" = rad.grid.use,
    "rad.grid.Z" = rad.grid.Z
  ))

}

CalculateRd <- function(prey.length.vec_m, fish.length.vec_m) {
  
  n.fls <- length(fish.length.vec_m)

  rd <- matrix(
    vapply(
      prey.length.vec_m,
      function(x) {
        rd.row <- 0.12 * (x * 1000) * (1 - exp(-20 * fish.length.vec_m))
      },
      FUN.VALUE = numeric(n.fls)
    ),
    nrow = n.fls
  )
  rownames(rd) <- fish.length.vec_m
  colnames(rd) <- prey.length.vec_m

  return(rd)

}

CalculateVmax <- function(fish.length.vec_m) {

  # max swim speed for the fish lengths in the simulation 
  vmax <- 0.87 * (fish.length.vec_m ** 0.19)  # in m/s
  names(vmax) <- fish.length.vec_m

  return(vmax)

}


GetMcdMat_01 <- function(vmax, rd, ca.prop, radial.grid.mean.normal.vels,
  radial.grid.norm.vels) {

  # calculate the max capture distance matrix
  mcd.mat <- ca.prop * sqrt(
    rd * rd * (vmax * vmax - radial.grid.mean.normal.vels * radial.grid.mean.normal.vels) / 
      (vmax * vmax + radial.grid.norm.vels * radial.grid.norm.vels - 
        radial.grid.mean.normal.vels * radial.grid.mean.normal.vels)
  )

  return(mcd.mat)

}

# GetFishSuccessMat_09 <- function(x, radial.step.dist.mat, radial.wet.node.mat) {
GetFishSuccessMat_09 <- function(mcd.mat, radial.step.dist.mat, radial.wet.node.mat) {

  # fish could be successful at any contiguous radial points where
  #  the mcd values were greater than the actual radial distances;
  #  we multiply by rad.grid$rad.grid.use to eliminate radial points that are
  #  out of the channel or water
  fish.success.mat <- (mcd.mat > radial.step.dist.mat) * radial.wet.node.mat
  # fish.success.mat
  
  # NAs could be generated if dist-weighted avg vel becomes too high somewhere
  #  on the radial; remove these by setting them equal to zero (because the
  #  fish can't make it to these spots)
  if (any(is.na(fish.success.mat))) {
    fish.success.mat[which(is.na(fish.success.mat))] <- 0L
  }

  # sometimes a fish crossing from slow to fast to slow waters will generate
  #  a dist-weighted avg vel that goes up then down again; sometimes these
  #  up/down patterns result in a dist-weighted avg vel that makes it look
  #  like the fish can make it, then not make it, then make it again; this is
  #  an artifact of the dist-weighted velocity approach; remove these false
  #  positives by setting them to zero too
  fsm.dims <- dim(fish.success.mat)
  for(my.row in 1:fsm.dims[1]) {
    
    # if there are both ones and zeros
    if (any(fish.success.mat[my.row, ] == 0)) {
    
      first.zero.ind <- min(which(fish.success.mat[my.row, ] == 0))
      fish.success.mat[my.row, first.zero.ind:fsm.dims[2]] <- 0L
            
    }

  }

  return(fish.success.mat)

}

GetCapAreaQ <- function(cap.area.success.mat, radial.grid.area.mat,
  radial.grid.velocities) {

  rad.grid.discharges <- cap.area.success.mat * radial.grid.area.mat * 
    radial.grid.velocities
  cap.area.discharge <- sum(rad.grid.discharges)

  return(cap.area.discharge)

}

#---------------------------------------------------------------------------
# Function to 
#  
# Args:
#  prey.energy.density_Jpg: J/g wet body mass of prey items
#  ...
#
# Returns:
#   
#---------------------------------------------------------------------------

CalculateCmax_FishBioE3_ConsEq3_RandParams <- function(water.temp_C, 
  fish.mass_g, prey.energy.density_Jpg, feeding.period_h) {
  
  # consumption equation 3 from fishBioE 3.0 for cool/cold water species
  
  # paremeters from Rand 1993
  kCA <- 0.628      # intercept of mass dependence funtion for a 1 gram fish at the optimum h20 temp
  kCB <- -0.3       # slope of allometric mass function (aka coefficient of mass dependence)
  kCQ <- 5    
  kCTO <- 20        # h20 temp correspoinding to 0.98 of the max consump. rate
  kCTM <- 20
  kCTL <- 24
  kCK1 <- 0.33
  kCK4 <- 0.2

  G1 <- (1 / (kCTO - kCQ)) * (log((0.98 * (1 - kCK1)) / (kCK1 * 0.02)))
  L1 <- exp(G1 * (water.temp_C - kCQ))
  KA <- (kCK1 * L1) / (1 + kCK1 * (L1 - 1))
  G2 <- (1 / (kCTL - kCTM)) * (log((0.98 * (1 - kCK4)) / (kCK4 * 0.02)))
  L2 <- exp(G2 * (kCTL - water.temp_C))
  KB <- (kCK4 * L2) / (1 + kCK4 * (L2 - 1))

  # maximal consumption, unconstrained by temperature effects
  c.max <- kCA * (fish.mass_g ** kCB)   #max specific feeding rate (g_prey/g_pred/d)
  
  # consumption, constrained by temperature and pval
  #  units of (g_prey / g_pred / day); p-val = 1
  spec.cons.rate <- (c.max * 1.0 * KA * KB)
  
  # converting spec.cons.rate to J/h to be compatible with units of grei and sc;
  #  note: might want to change this from 24 hours to 8 or 10 in the future as
  #    24 likely underestimates the maximal feeding rate; check the FishBioE 3.0
  #    literature for confirmation.
  cons.jh <- spec.cons.rate * prey.energy.density_Jpg * fish.mass_g / feeding.period_h  # in J/h
  
  return(cons.jh)
  
}


#---------------------------------------------------------------------------
# Function to 
#  
# Args:
#  prey.energy.density_Jpg: J/g wet body mass of prey items
#  ...
#
# Returns:
#   
#---------------------------------------------------------------------------
CalculateCmax_FishBioE3_ConsEq3_RailsRoseParams <- function(water.temp_C, 
  fish.mass_g, prey.energy.density_Jpg, feeding.period_h) {
  
  # consumption equation 3 from fishBioE 3.0 for cool/cold water species
  
  # parameters from fishBioE 3.0 and Railsback and Rose 1999
  kCA <- 0.628    # intercept of mass dependence funtion for a 1 gram fish at the optimum h20 temp
  kCB <- -0.3     # slope of allometric mass function (aka coefficient of mass dependence)
  kCQ <- 3.5    
  kCTO <- 25      # h20 temp correspoinding to 0.98 of the max consump. rate
  kCTM <- 22.5
  kCTL <- 24.3
  kCK1 <- 0.2
  kCK4 <- 0.2

  G1 <- (1 / (kCTO - kCQ)) * (log((0.98 * (1 - kCK1)) / (kCK1 * 0.02)))
  L1 <- exp(G1 * (water.temp_C - kCQ))
  KA <- (kCK1 * L1) / (1 + kCK1 * (L1 - 1))
  G2 <- (1 / (kCTL - kCTM)) * (log((0.98 * (1 - kCK4)) / (kCK4 * 0.02)))
  L2 <- exp(G2 * (kCTL - water.temp_C))
  KB <- (kCK4 * L2) / (1 + kCK4 * (L2 - 1))

  # maximal consumption, unconstrained by temperature effects
  c.max <- kCA * (fish.mass_g ** kCB)   #max specific feeding rate (g_prey/g_pred/d)
  
  # consumption, constrained by temperature and pval
  #  units of (g_prey / g_pred / day); p-val = 1
  spec.cons.rate <- (c.max * 1.0 * KA * KB)
  
  # converting spec.cons.rate to J/h to be compatible with units of grei and sc;
  #  note: might want to change this from 24 hours to 8 or 10 in the future as
  #    24 likely underestimates the maximal feeding rate; check the FishBioE 3.0
  #    literature for confirmation.
  cons.jh <- spec.cons.rate * prey.energy.density_Jpg * fish.mass_g / feeding.period_h  # in J/h
  
  return(cons.jh)
  
}

# returns cmax-limited if you input cmax; returns specific idxs if you input 
#  them
CalculateGrei <- function(cap.area.discharge.fn, c.max_Jph = NULL,  
  pt.ids = NULL, prey.length_m, fish.length_m, prey.conc_noPerM3,
  prey.energy.content_J, losses.to.waste_prop) {

  require(readr)
  require(dplyr)

  cad.df <- read_csv(cap.area.discharge.fn, col_types = 'idddid')

  if (is.null(c.max_Jph)) {

    if (is.null(pt.ids)) {

      cad.dfs <- cad.df %>%
        filter(
          near(pl_m, prey.length_m),
          near(fl_m, fish.length_m)
        ) %>%
        mutate(
          unlim.grei_Jph = cap.area.q_m3ps * prey.conc_noPerM3 * 
            prey.energy.content_J * (1 - losses.to.waste_prop) * 3600
        ) %>%
        select(
          idx, ca.prop, max.radial.step, unlim.grei_Jph
        )

    } else {

      if(!is.integer(pt.ids)) {

        pt.ids <- as.integer(pt.ids)

      }

      cad.dfs <- cad.df %>%
        filter(
          idx %in% pt.ids,
          near(pl_m, prey.length_m),
          near(fl_m, fish.length_m)
        ) %>%
        mutate(
          unlim.grei_Jph = cap.area.q_m3ps * prey.conc_noPerM3 * 
            prey.energy.content_J * (1 - losses.to.waste_prop) * 3600
        ) %>% 
        select(
          idx, ca.prop, max.radial.step, unlim.grei_Jph
        )

    }

  } else {

    if (is.null(pt.ids)) {

      cad.dfs <- cad.df %>%
        filter(
          near(pl_m, prey.length_m),
          near(fl_m, fish.length_m)
        ) %>%
        mutate(
          unlim.grei_Jph = cap.area.q_m3ps * prey.conc_noPerM3 * 
            prey.energy.content_J * (1 - losses.to.waste_prop) * 3600,
          lim.grei_Jph = pmin(
            unlim.grei_Jph,
            (1 - losses.to.waste_prop) * c.max_Jph
          )
        ) %>%
        select(
          idx, ca.prop, max.radial.step, unlim.grei_Jph, lim.grei_Jph
        )

    } else {

      if(!is.integer(pt.ids)) {

        pt.ids <- as.integer(pt.ids)

      }
      
      cad.dfs <- cad.df %>%
        filter(
          idx %in% pt.ids,
          near(pl_m, prey.length_m),
          near(fl_m, fish.length_m)
        ) %>%
        mutate(
          unlim.grei_Jph = cap.area.q_m3ps * prey.conc_noPerM3 * 
            prey.energy.content_J * (1 - losses.to.waste_prop) * 3600,
          lim.grei_Jph = pmin(
            unlim.grei_Jph,
            (1 - losses.to.waste_prop) * c.max_Jph
          )
        ) %>% 
        select(
          idx, ca.prop, max.radial.step, unlim.grei_Jph, lim.grei_Jph
        )

    }

  }

  return(cad.dfs)

}


CalculateMinRadialStepToAchievePval <- function(prey.energy.density_Jpg,
  water.temp_C, fish.mass_g, feeding.period_h, cap.area.discharge.fn, 
  c.max_Jph = NULL, pt.ids = NULL, prey.length_m, fish.length_m,
  prey.conc_noPerM3, prey.energy.content_J, losses.to.waste_prop, user.pval) {

  require(readr)
  require(dplyr)

  c.max_Jph <- CalculateCmax_FishBioE3_ConsEq3_RandParams(
    prey.energy.density_Jpg =prey.energy.density_Jpg , 
    water.temp_C = water.temp_C, 
    fish.mass_g = fish.mass_g, 
    feeding.period_h = feeding.period_h
  )

  grei.df <- CalculateGrei(
    cap.area.discharge.fn = cap.area.discharge.fn, 
    prey.length_m = prey.length_m,
    fish.length_m = fish.length_m, 
    prey.conc_noPerM3 = prey.conc_noPerM3, 
    prey.energy.content_J = prey.energy.content_J,
    losses.to.waste_prop = losses.to.waste_prop
  )

  rad.step.gte.user.pval.df <- grei.df %>%
    # filter(idx %in% 10:20) %>%
    mutate(
      unlim.grei.gte.user.pval = unlim.grei_Jph >= (1 - losses.to.waste_prop) * c.max_Jph * user.pval
    ) %>%
    group_by(idx) %>%
    summarize(
      rad.step.gte.user.pval = max.radial.step[min(which(unlim.grei.gte.user.pval))]
    )

  return(rad.step.gte.user.pval.df)

}

#---------------------------------------------------------------------------
# Function to calculate energy expenditure due to holding position at a 
#  foraging location (based on constant velocity at the focal point).
#  Implemented as in Hayes/Hughes 2007 and the Cawthron user's manual
#  appendix for steelhead swim costs; essentially, this is a modification
#  of Fish BioE 3.0 respiration equation 2; instead of a constant
#  activity multiplier, they have substituted the exponential
#  activity multiplier from respiration equation 1  
#
# Args:
#  
# Returns:
#   
#---------------------------------------------------------------------------

CalculateSwimCosts_Hayes07_sthd <- function(velocity.mag.fn, pt.ids = NULL,
  water.temp_C, fish.mass_g) {

  require(readr)
  require(dplyr)

  vel.mag.df <- read_csv(velocity.mag.fn, guess_max = 10000)

  if (!is.null(pt.ids)) {

    if (!is.integer(pt.ids)) {

      pt.ids <- as.integer(pt.ids)

    }

    vel.mag.df <- vel.mag.df %>%
      filter(idx %in% pt.ids)

  }
  
  # parameters
  kRA <- 0.013
  kRB <- -0.217
  kRQ2 <- 2.2
  kRT <- 0.0234
  kRTO <- 22
  kRTM <- 26
  
  V <- (kRTM - water.temp_C)/(kRTM - kRTO)  # this is not velocity
  Z <- (log(kRQ2)) * (kRTM - kRTO)
  Y <- (log(kRQ2)) * (kRTM - kRTO + 2)
  X <- ((Z ** 2) * (1 + (1 + 40 / Y) ** 0.5) ** 2) / 400

  temp.func <- (V ** X) * (exp(X * (1 - V)))
  act <- exp(kRT * (vel.mag.df$Velocity.Magnitude * 100))  # vels in cm/s

  # swim costs in g O2/g fish/day
  swim.costs.hh.ggd <- kRA * (fish.mass_g ** kRB) * temp.func * act

  # convert g O2/g fish/day to J/h
  vel.mag.df$swim.costs_Jph <- swim.costs.hh.ggd * fish.mass_g * 13565 / 24
  
  return(
    vel.mag.df %>%
      select(idx, Velocity.Magnitude, swim.costs_Jph)
  )

}


CalculateNrei <- function(prey.energy.density_Jpg, water.temp_C, fish.mass_g, 
  feeding.period_h, cap.area.discharge.fn, prey.length_m, fish.length_m, 
  prey.conc_noPerM3, prey.energy.content_J, losses.to.waste_prop, 
  velocity.mag.fn) {

  require(readr)
  require(dplyr)

  c.max_Jph <- CalculateCmax_FishBioE3_ConsEq3_RandParams(
    prey.energy.density_Jpg = prey.energy.density_Jpg, 
    water.temp_C = water.temp_C,
    fish.mass_g = fish.mass_g, 
    feeding.period_h = feeding.period_h
  )

  grei.df <- CalculateGrei(
    cap.area.discharge.fn = cap.area.discharge.fn,
    prey.length_m = prey.length_m,
    fish.length_m = fish.length_m,
    prey.conc_noPerM3 = prey.conc_noPerM3,
    prey.energy.content_J = prey.energy.content_J,
    losses.to.waste_prop = losses.to.waste_prop,
    c.max_Jph = c.max_Jph
  )

  sc.df <- CalculateSwimCosts_Hayes07_sthd(
    velocity.mag.fn = velocity.mag.fn,
    water.temp_C = water.temp_C, 
    fish.mass_g = fish.mass_g
  )

  nrei.df <- merge(
    grei.df,
    sc.df,
    all.x = TRUE,
    sort = FALSE
  )

  nrei.df <- nrei.df %>%
    mutate(nrei_Jph = lim.grei_Jph - swim.costs_Jph)

  return(nrei.df)

}

CalculateMinRadialStepToAchieveNreiThresh <- function(prey.energy.density_Jpg, 
  water.temp_C, fish.mass_g, feeding.period_h, cap.area.discharge.fn,
  prey.length_m, fish.length_m, prey.conc_noPerM3, prey.energy.content_J,
  losses.to.waste_prop, velocity.mag.fn, user.nrei.thresh_Jph) {

  require(readr)
  require(dplyr)

  c.max_Jph <- CalculateCmax_FishBioE3_ConsEq3_RandParams(
    prey.energy.density_Jpg = prey.energy.density_Jpg, 
    water.temp_C = water.temp_C,
    fish.mass_g = fish.mass_g, 
    feeding.period_h = feeding.period_h
  )

  grei.df <- CalculateGrei(
    cap.area.discharge.fn = cap.area.discharge.fn,
    c.max_Jph = c.max_Jph,
    prey.length_m = prey.length_m,
    fish.length_m = fish.length_m,
    prey.conc_noPerM3 = prey.conc_noPerM3,
    prey.energy.content_J = prey.energy.content_J,
    losses.to.waste_prop = losses.to.waste_prop
  )

  sc.df <- CalculateSwimCosts_Hayes07_sthd(
    velocity.mag.fn = velocity.mag.fn,
    water.temp_C = water.temp_C, 
    fish.mass_g = fish.mass_g
  )

  nrei.df <- merge(
    grei.df,
    sc.df,
    all.x = TRUE,
    sort = FALSE
  )

  nrei.df <- nrei.df %>%
    mutate(nrei_Jph = lim.grei_Jph - swim.costs_Jph)

  rad.step.gte.nrei.thresh.df <- nrei.df %>%
    # filter(idx %in% 10:20) %>%
    mutate(
      nrei.gte.user.nrei.thresh = nrei_Jph >= user.nrei.thresh_Jph
    ) %>%
    group_by(idx) %>%
    summarize(
      rad.step.gte.user.nrei.thresh = max.radial.step[min(which(nrei.gte.user.nrei.thresh))]
    )

  return(rad.step.gte.nrei.thresh.df)

}


#---------------------------------------------------------------------------
# Function to calculate...
#
# based on Keeley and McPhail '98, adjusted based on Imre, Grant, Keeley '04    
#
# Args:
#  
# Returns:
#   
#---------------------------------------------------------------------------

CalculateFishTerritory <- function(fish.length_m) {

  # fish length must be in cm; territory area is in m2;
  fish.territory.area <- 10 ** (1.56 * log10(fish.length_m * 100) - 1.81 - 0.07)

  # radius is in m
  fish.territory.radius <- sqrt(fish.territory.area / pi)

  return(
    list(
      area = fish.territory.area,
      radius = fish.territory.radius
    )
  )
    
}



CreateFishPlacementDf <- function(prey.energy.density_Jpg,
  water.temp_C, fish.mass_g, feeding.period_h, cap.area.discharge.fn, 
  c.max_Jph = NULL, pt.ids = NULL, prey.length_m, fish.length_m,
  prey.conc_noPerM3, prey.energy.content_J, losses.to.waste_prop, user.pval,
  velocity.mag.fn, user.nrei.thresh_Jph, exp.steepness.coeff, max.terr.length_m) {

  require(readr)
  require(dplyr)

  c.max_Jph <- CalculateCmax_FishBioE3_ConsEq3_RandParams(
    prey.energy.density_Jpg =prey.energy.density_Jpg , 
    water.temp_C = water.temp_C, 
    fish.mass_g = fish.mass_g, 
    feeding.period_h = feeding.period_h
  )

  grei.df <- CalculateGrei(
    cap.area.discharge.fn = cap.area.discharge.fn, 
    prey.length_m = prey.length_m,
    fish.length_m = fish.length_m, 
    prey.conc_noPerM3 = prey.conc_noPerM3, 
    prey.energy.content_J = prey.energy.content_J,
    losses.to.waste_prop = losses.to.waste_prop
  )

  rad.step.gte.user.pval.df <- grei.df %>%
    # filter(idx %in% 10:20) %>%
    mutate(
      unlim.grei.gte.user.pval = unlim.grei_Jph >= (1 - losses.to.waste_prop) * c.max_Jph * user.pval
    ) %>%
    group_by(idx) %>%
    summarize(
      rad.step.gte.user.pval = max.radial.step[min(which(unlim.grei.gte.user.pval))]
    )

  nrei.df <- CalculateNrei(
    prey.energy.density_Jpg = prey.energy.density_Jpg,
    water.temp_C = water.temp_C,
    fish.mass_g = fish.mass_g,
    feeding.period_h = feeding.period_h,
    cap.area.discharge.fn = cap.area.discharge.fn,
    prey.length_m = prey.length_m,
    fish.length_m = fish.length_m,
    prey.conc_noPerM3 = prey.conc_noPerM3,
    prey.energy.content_J = prey.energy.content_J,
    losses.to.waste_prop = losses.to.waste_prop,
    velocity.mag.fn = velocity.mag.fn 
  )

  rad.step.gte.nrei.thresh.df <- CalculateMinRadialStepToAchieveNreiThresh(
    prey.energy.density_Jpg = prey.energy.density_Jpg,
    water.temp_C = water.temp_C,
    fish.mass_g = fish.mass_g,
    feeding.period_h = feeding.period_h,
    cap.area.discharge.fn = cap.area.discharge.fn,
    prey.length_m = prey.length_m,
    fish.length_m = fish.length_m,
    prey.conc_noPerM3 = prey.conc_noPerM3,
    prey.energy.content_J = prey.energy.content_J,
    losses.to.waste_prop = losses.to.waste_prop,
    velocity.mag.fn = velocity.mag.fn,
    user.nrei.thresh_Jph = user.nrei.thresh_Jph
  )

  plotting.df <- read_csv(velocity.mag.fn) %>%
    select(idx, X, Y, Z, zrank, Depth, X.Velocity, Y.Velocity, Velocity.Magnitude)

  placement.df <- merge(
    plotting.df,
    rad.step.gte.user.pval.df,
    all.y = TRUE,
    sort = FALSE
  )

  nrei.for.placement.df <- nrei.df %>%
    filter(
      near(ca.prop, 1.0)
    ) %>%
    select(
      idx, nrei_Jph
    )

  placement.df2 <- merge(
    placement.df,
    nrei.for.placement.df,
    all.y = TRUE,
    sort = FALSE
  )

  placement.df3 <- merge(
    placement.df2,
    rad.step.gte.nrei.thresh.df,
    all.x = TRUE,
    sort = FALSE
  )

  fish.terr <- CalculateFishTerritory(fish.length_m)
  vmax <- CalculateVmax(fish.length_m)

  placement.df4 <- placement.df3 %>%
    filter(
      nrei_Jph >= 0,
      !is.na(rad.step.gte.user.pval)
    ) %>%
    mutate(
      los.min.radius = pmax(rad.step.gte.user.pval, rad.step.gte.user.nrei.thresh) / 100,
      lateral.b1 = (los.min.radius - fish.terr$radius) / (exp(exp.steepness.coeff * vmax) - 1),
      lateral.b3 = fish.terr$radius - lateral.b1,
      los.radius = lateral.b1 * exp(exp.steepness.coeff * Velocity.Magnitude) + lateral.b3
    ) %>%
    select(
      -c(
        rad.step.gte.user.pval, rad.step.gte.user.nrei.thresh, los.min.radius, 
        lateral.b1, lateral.b3
      )
    ) %>%
    mutate(
      long.b1 = (los.radius - max.terr.length_m) / (exp(exp.steepness.coeff * vmax) - 1),
      long.b3 = max.terr.length_m - long.b1,
      los.cent.to.cent.dist = long.b1 * exp(exp.steepness.coeff * Velocity.Magnitude) + long.b3,
      vel.angle = atan2(Y.Velocity, X.Velocity),
      cyl.end.delta.x = cos(vel.angle) * los.cent.to.cent.dist,
      cyl.end.delta.y = sin(vel.angle) * los.cent.to.cent.dist,
      cyl.end.delta.z = 0,
      cyl.end.x.coord = X + cyl.end.delta.x,
      cyl.end.y.coord = Y + cyl.end.delta.y,
      cyl.end.z.coord = Z + cyl.end.delta.z
    ) %>%
    select(
      -c(long.b1, long.b3, cyl.end.delta.x, cyl.end.delta.y, cyl.end.delta.z)
    ) %>%
    arrange(desc(nrei_Jph)) 

  return(placement.df4)

}

CalculateFishLocations <- function(prey.energy.density_Jpg,
  water.temp_C, fish.mass_g, feeding.period_h, cap.area.discharge.fn, 
  c.max_Jph = NULL, pt.ids = NULL, prey.length_m, fish.length_m,
  prey.conc_noPerM3, prey.energy.content_J, losses.to.waste_prop, user.pval,
  velocity.mag.fn, user.nrei.thresh_Jph, Dist, exp.steepness.coeff,
  max.terr.length_m, los.cent.line.max.step_m, plots = FALSE,
  plot.increment = NULL, placement.messages = FALSE) {

  require(readr)
  require(dplyr)
  require(fields)

  # calculate cmax
  c.max_Jph <- CalculateCmax_FishBioE3_ConsEq3_RandParams(
    prey.energy.density_Jpg = prey.energy.density_Jpg , 
    water.temp_C = water.temp_C, 
    fish.mass_g = fish.mass_g, 
    feeding.period_h = feeding.period_h
  )

  # calculate grei
  grei.df <- CalculateGrei(
    cap.area.discharge.fn = cap.area.discharge.fn, 
    prey.length_m = prey.length_m,
    fish.length_m = fish.length_m, 
    prey.conc_noPerM3 = prey.conc_noPerM3, 
    prey.energy.content_J = prey.energy.content_J,
    losses.to.waste_prop = losses.to.waste_prop
  )
  
  # find the radial step where unlimited grei exceeds a user-inputted 
  #  proportion of maximal consumption
  rad.step.gte.user.pval.df <- grei.df %>%
    # filter(idx %in% 10:20) %>%
    mutate(
      unlim.grei.gte.user.pval = unlim.grei_Jph >= (1 - losses.to.waste_prop) * c.max_Jph * user.pval
    ) %>%
    group_by(idx) %>%
    summarize(
      rad.step.gte.user.pval = max.radial.step[min(which(unlim.grei.gte.user.pval))]
    )

  # calculate nrei
  nrei.df <- CalculateNrei(
    prey.energy.density_Jpg = prey.energy.density_Jpg,
    water.temp_C = water.temp_C,
    fish.mass_g = fish.mass_g,
    feeding.period_h = feeding.period_h,
    cap.area.discharge.fn = cap.area.discharge.fn,
    prey.length_m = prey.length_m,
    fish.length_m = fish.length_m,
    prey.conc_noPerM3 = prey.conc_noPerM3,
    prey.energy.content_J = prey.energy.content_J,
    losses.to.waste_prop = losses.to.waste_prop,
    velocity.mag.fn = velocity.mag.fn 
  )
  # head(nrei.df, 100)

  # find the radial step where nrei exceeds a user-inputted value
  rad.step.gte.nrei.thresh.df <- CalculateMinRadialStepToAchieveNreiThresh(
    prey.energy.density_Jpg = prey.energy.density_Jpg,
    water.temp_C = water.temp_C,
    fish.mass_g = fish.mass_g,
    feeding.period_h = feeding.period_h,
    cap.area.discharge.fn = cap.area.discharge.fn,
    prey.length_m = prey.length_m,
    fish.length_m = fish.length_m,
    prey.conc_noPerM3 = prey.conc_noPerM3,
    prey.energy.content_J = prey.energy.content_J,
    losses.to.waste_prop = losses.to.waste_prop,
    velocity.mag.fn = velocity.mag.fn,
    user.nrei.thresh_Jph = user.nrei.thresh_Jph
  )

  # begin constructing placement dataframe
  plotting.df <- read_csv(velocity.mag.fn, guess_max = 10000) %>%
    select(idx, X, Y, Z, zrank, Depth, X.Velocity, Y.Velocity, Velocity.Magnitude)

  # add the radial step at which the fish exceeds the user-specified pval
  placement.df1 <- merge(
    plotting.df,
    rad.step.gte.user.pval.df,
    all.y = TRUE,
    sort = FALSE
  )

  # fish placement is based on the full-mcd nrei
  nrei.for.placement.df <- nrei.df %>%
    filter(
      near(ca.prop, 1.0)
    ) %>%
    select(
      idx, nrei_Jph
    )

  # add in full-mcd nrei values
  placement.df2 <- merge(
    placement.df1,
    nrei.for.placement.df,
    all.y = TRUE,
    sort = FALSE
  )
  
  # add in the radial step at which the fish exceeds the user-specified nrei
  #  threshold; has idx, x, y, z, zrank, depth, vel, and nrei => will serve
  #  as a plotting df
  placement.df3 <- merge(
    placement.df2,
    rad.step.gte.nrei.thresh.df,
    all.x = TRUE,
    sort = FALSE
  )
  
  # calculate fish-related metrics for this fish and prey size
  fish.terr <- CalculateFishTerritory(fish.length_m)
  vmax <- CalculateVmax(fish.length_m)
  rd <- CalculateRd(prey.length.vec_m = prey.length_m, 
    fish.length.vec_m = fish.length_m)

  # LOS notes
  # min radius for LOSs will be the radial step at which a user-spec pval is achieved
  # b1 + b3 is the y-intercept
  # we will solve for b1 and b3 such that min radius for LOSs is a function of velocity
  #  bracketed between a fish-length-based territory size and the radial step at which
  #  the user-spec pval is achieved; low velocity = wider spheres; high velocity = narrower spheres
  # in the longitudinal direciton, we will bracket between a user-spec max territory
  #  length and the calculated LOS sphere radius.  So, we will either have a line of spheres
  #  or a single sphere  

  # finished placement df
  pdf <- placement.df3 %>%
    filter(
      nrei_Jph >= user.nrei.thresh_Jph,
      !is.na(rad.step.gte.user.pval)
    ) %>%
    mutate(
      los.min.radius = pmax(rad.step.gte.user.pval, rad.step.gte.user.nrei.thresh) * Dist,
      lateral.b1 = (los.min.radius - fish.terr$radius) / (exp(exp.steepness.coeff * vmax) - 1),
      lateral.b3 = fish.terr$radius - lateral.b1,
      los.radius_m = lateral.b1 * exp(exp.steepness.coeff * Velocity.Magnitude) + lateral.b3
    ) %>%
    select(
      -c(
        rad.step.gte.user.pval, rad.step.gte.user.nrei.thresh, los.min.radius, 
        lateral.b1, lateral.b3
      )
    ) %>%
    mutate(
      long.b1 = (los.radius_m - max.terr.length_m) / (exp(exp.steepness.coeff * vmax) - 1),
      long.b3 = max.terr.length_m - long.b1,
      los.cent.to.cent.dist = long.b1 * exp(exp.steepness.coeff * Velocity.Magnitude) + long.b3,
      vel.angle = atan2(Y.Velocity, X.Velocity),
      cyl.end.delta.x = cos(vel.angle) * los.cent.to.cent.dist,
      cyl.end.delta.y = sin(vel.angle) * los.cent.to.cent.dist,
      cyl.end.delta.z = 0,
      cyl.end.x.coord = X + cyl.end.delta.x,
      cyl.end.y.coord = Y + cyl.end.delta.y,
      cyl.end.z.coord = Z + cyl.end.delta.z
    ) %>%
    select(
      -c(long.b1, long.b3, cyl.end.delta.x, cyl.end.delta.y, cyl.end.delta.z)
    ) %>%
    arrange(desc(nrei_Jph))
  # pdf %>%
  #   filter(idx == 24720)

  pdf.rows <- nrow(pdf)
  
  if (nrow(pdf) == 0) {

    print('no fish')
    global.fish.df <- pdf
    global.fish.terr.df <- NULL
    global.buffer.df <- NULL

  } else if (nrow(pdf) == 1) {

    print('one fish')
    global.fish.df <- pdf[1, ]
    global.fish.terr.df <- NULL
    global.buffer.df <- NULL

  } else {

    print('attempting to place fish in the reach')

    # df for candidate fish
    fish.df <- pdf[1, ]

    # determine search radius to identify nearby points
    search.rad <- fish.df$los.cent.to.cent.dist + (2 * fish.df$los.radius_m) + (2 * rd[1, 1])

    # get nearby points
    close.df <- pdf %>%
      filter(
        abs(X - fish.df$X) <= search.rad,
        abs(Y - fish.df$Y) <= search.rad  # adding Z doesn't save much
      )

    # calculate coordinates of points along the los centerline
    num.steps <- ceiling(fish.df$los.cent.to.cent.dist / los.cent.line.max.step_m)
    true.step.length <- fish.df$los.cent.to.cent.dist / num.steps
    cent.line.pts <- data.frame(
      X = fish.df$X + ((0:num.steps * true.step.length) * cos(fish.df$vel.angle)),
      Y = fish.df$Y + ((0:num.steps * true.step.length) * sin(fish.df$vel.angle)),
      Z = fish.df$Z + ((0:num.steps * true.step.length) * 0)
    )

    # calculate distances betweeen close.df and cent.line.pts
    dists <- rdist(
      cent.line.pts,
      close.df[, c('X', 'Y', 'Z')]
    )

    #-----------------------------------------------------------
    # find points to exclude from future fish placement
    #-----------------------------------------------------------

    # future fish must be rd away from any previous fish's territory; therefore,
    #  find and exclude points within los.radius_m + rd
    all.excluded.inds <- unique(which(dists <= (fish.df$los.radius_m + rd[1, 1]), 
      arr.ind = TRUE)[, 'col'])
    all.excluded.idxs <- close.df[all.excluded.inds, 'idx']

    # get idxs of points in the fish's territory 
    fish.terr.inds <- unique(which(dists <= fish.df$los.radius_m, 
      arr.ind = TRUE)[, 'col'])
    fish.terr.idxs <- close.df[fish.terr.inds, 'idx']
    
    # buffer idxs = those too close to an existing territory for future 
    #  consideration of fish placement
    if (any(!all.excluded.idxs %in% fish.terr.idxs)) {

      buffer.idxs <- all.excluded.idxs[!all.excluded.idxs %in% fish.terr.idxs]

    } else {

      buffer.idxs <- NA

    }

    #-----------------------------------------------------------
    # establish global dfs
    #-----------------------------------------------------------

    global.fish.df <- fish.df
    
    fish.terr.df <- data.frame(
      idx = fish.terr.idxs,
      fish.no = nrow(global.fish.df)
    )

    buffer.df <- data.frame(
      idx = buffer.idxs,
      fish.no = nrow(global.fish.df)
    )

    global.fish.terr.df <- fish.terr.df
    global.buffer.df <- buffer.df
    
    working.pdf <- pdf[!pdf$idx %in% all.excluded.idxs, ]

    #-----------------------------------------------------------
    # loop through potential fish 2+
    #-----------------------------------------------------------

    while (nrow(working.pdf) > 0) {

      print(paste('nrow of working.pdf = ', nrow(working.pdf), sep = ''))
      
      # new candidate fish
      fish.df <- working.pdf[1, ]

      # determine search radius to identify nearby points
      search.rad <- fish.df$los.cent.to.cent.dist + (2 * fish.df$los.radius_m) + (2 * rd[1, 1])

      # get nearby points
      close.df <- pdf %>%
        filter(
          abs(X - fish.df$X) <= search.rad,
          abs(Y - fish.df$Y) <= search.rad  # adding Z doesn't save much
        )

      # calculate coordinates of points along the los centerline
      num.steps <- ceiling(fish.df$los.cent.to.cent.dist / los.cent.line.max.step_m)
      true.step.length <- fish.df$los.cent.to.cent.dist / num.steps
      cent.line.pts <- data.frame(
        X = fish.df$X + ((0:num.steps * true.step.length) * cos(fish.df$vel.angle)),
        Y = fish.df$Y + ((0:num.steps * true.step.length) * sin(fish.df$vel.angle)),
        Z = fish.df$Z + ((0:num.steps * true.step.length) * 0)
      )

      # calculate distances betweeen close.df and cent.line.pts
      dists <- rdist(
        cent.line.pts,
        close.df[, c('X', 'Y', 'Z')]
      )

      #-----------------------------------------------------------
      # find points to exclude from future fish placement
      #-----------------------------------------------------------

      # future fish must be rd away from any previous fish's territory; therefore,
      #  find and exclude points within los.radius_m + rd
      all.excluded.inds <- unique(which(dists <= (fish.df$los.radius_m + rd[1, 1]), 
        arr.ind = TRUE)[, 'col'])
      all.excluded.idxs <- close.df[all.excluded.inds, 'idx']

      # get idxs of points in the fish's territory 
      fish.terr.inds <- unique(which(dists <= fish.df$los.radius_m, 
        arr.ind = TRUE)[, 'col'])
      fish.terr.idxs <- close.df[fish.terr.inds, 'idx']

      # candidate fish's territory points must be at least rd from any
      #  previously placed fish
      fish.terr.pts <- close.df[fish.terr.inds, c('X', 'Y', 'Z')]
      dists.2 <- rdist(
        global.fish.df[, c('X', 'Y', 'Z')],
        fish.terr.pts,
      )
      # min(dists.2)
      
      # buffer idxs = those too close to an existing territory for future 
      #  consideration of fish placement
      if (any(!all.excluded.idxs %in% fish.terr.idxs)) {

        buffer.idxs <- all.excluded.idxs[!all.excluded.idxs %in% fish.terr.idxs]

      } else {

        buffer.idxs <- NA

      }

      # if this fish's territory doesn't overlap a previously-placed fish AND
      #  its territory is at least rd from any previously placed fish,
      #  add it; else, remove the row and go on to the next candidate
      if (all(!fish.terr.idxs %in% global.fish.df$idx) && all(dists.2 > rd[1, 1])) {

        global.fish.df <- rbind(global.fish.df, fish.df)
        
        fish.terr.df <- data.frame(
          idx = fish.terr.idxs,
          fish.no = nrow(global.fish.df)
        )

        buffer.df <- data.frame(
          idx = buffer.idxs,
          fish.no = nrow(global.fish.df)
        )

        global.fish.terr.df <- rbind(global.fish.terr.df, fish.terr.df)
        global.buffer.df <- rbind(global.buffer.df, buffer.df)        
        
        working.pdf <- working.pdf[!working.pdf$idx %in% all.excluded.idxs, ]

        #-----------------------------------------------------------
        # optional plotting check
        #-----------------------------------------------------------        
        if (plots && nrow(global.fish.df) %% plot.increment == 0) {

          par(mfrow = c(1, 2))
          
          # nearby points, territory, and buffer for this fish
          # local view
          plot(close.df$X, close.df$Y, asp = 1, col = 'lightgray')
          points(
            close.df$X[close.df$idx %in% buffer.df$idx],
            close.df$Y[close.df$idx %in% buffer.df$idx],
            col = 'brown',
            pch = '*'
          )
          points(
            close.df$X[close.df$idx %in% fish.terr.df$idx],
            close.df$Y[close.df$idx %in% fish.terr.df$idx],
            col = 'yellow'
          )
          points(fish.df$X, fish.df$Y, col = 'red')
          
          # whole-reach view
          plot(pdf$X, pdf$Y, asp = 1, col = 'lightgray')
          points(
            pdf$X[pdf$idx %in% global.buffer.df$idx],
            pdf$Y[pdf$idx %in% global.buffer.df$idx],
            col = 'brown',
            pch = '*'
          )
          points(
            pdf$X[pdf$idx %in% global.fish.terr.df$idx],
            pdf$Y[pdf$idx %in% global.fish.terr.df$idx],
            col = 'yellow'
          )
          points(global.fish.df$X, global.fish.df$Y, col = as.factor(global.fish.df$idx))  

        }

      } else {

        if (placement.messages) {
        
          if (any(fish.terr.idxs %in% global.fish.df$idx) && all(dists.2 > rd[1, 1])) {

            print(
              paste(
                'territory for idx', fish.df$idx,
                'overlaps an existing fish so a fish cannot be placed here',
                sep = ' '
              )
            )

          } else if (all(!fish.terr.idxs %in% global.fish.df$idx) && any(dists.2 < rd[1, 1])) {

            print(
              paste(
                'territory for idx', fish.df$idx,
                'is too close to an existing fish so a fish cannot be placed here',
                sep = ' '
              )
            )

          } else if (any(fish.terr.idxs %in% global.fish.df$idx) && any(dists.2 < rd[1, 1])) {

            print(
              paste(
                'territory for idx', fish.df$idx, 'both overlaps an existing fish',
                'AND is too close to an existing fish so a fish cannot be placed here',
                sep = ' '
              )
            )

          } else {

            print('Something else went wrong')

          }

        }

        print(paste('removing idx', fish.df$idx, 'from consideration', sep = ' '))
        working.pdf <- working.pdf[-1, ]

      }  # end of if the candidate point doesn't satisfy our rules

    }  # end of while loop for fish 2+

  }  # end of if, if else, else statement for nrow(pdf)

  return(
    list(
      nrei.df = nrei.df,
      plotting.df = placement.df3,
      placement.df = pdf,
      fish.df = global.fish.df,
      fish.terr.df = global.fish.terr.df,
      buffer.pts.df = global.buffer.df
    )
  )

}  # end of function def


CreateMultiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# determine grd.reduction.factor for ilt
SetGRF <- function(DEMfilename, grf.start, dx_dy) {

  print("Reading Hydraulic Model Data File (This may take a minute!)")

  # Read the data.  Screen to only wetted areas.
  require(data.table)
  DEM.Results = fread(DEMfilename, header = TRUE, data.table = FALSE)

  dz.start <- dx_dy * grf.start

  approx.grid.size <- sum(floor(DEM.Results$Depth / dz.start)) / (grf.start * grf.start)

  while(approx.grid.size > 500000) {

    grf.start <- grf.start + 1
    dz.start <- dx_dy * grf.start

    approx.grid.size <- floor(sum(floor(DEM.Results$Depth / dz.start)) / (grf.start * grf.start))

  }

  data.frame(gf.end = grf.start, dz.end = dz.start, ags = approx.grid.size)

}