#################################################################################
### To run, this code requires:
  # data_L: nxd matrix of data (transformed to Laplace margins - see separate R file) 
  # loc: dx2 matrix of locations

#################################################################################
### Fix necessary model parameters:
  # Choose number of latent variables in alpha and gamma spline functions
    knot.num <- 14 # corresponds to 15 knots including boundary knots (with spline constrained to 0 at 0)
  # Choose degree in spline functions ###
    spline.deg <- 2 # quadratic B-splines
  # Choose SPDE parameter alpha = nu + D/2
    alpha <- 1.5

#################################################################################
### Select conditioning site and extract extreme observations:
  # Set threshold
    u <- quantile(data_L, 0.95)
    
  # Choose conditioning site and place in first column. Order columns in increasing distance.
    cond <- 3363
    dist2s0 <- as.matrix(dist(loc))[cond,] # Distance of locations from conditioning site
    ord <- sort.int(dist2s0,index.return = T)$ix
    
  # Extract exceedances at conditioning site
    data_L <- data_L[,ord]
    loc <- loc[ord,]
    obs <- data_L[data_L[,1]>u,]

##########################################################################################
### Set up data for use with R-INLA
  # Number of observations and locations
    n.repl <- nrow(obs)
    nsites <- ncol(obs)
    
  # First site fixed as the conditioning location
    dist2s0=as.matrix(dist(loc))[1,]
    y.cond <- obs[,1]
    
  # set INLA data (can set to NA if used for validation) 
    yinla = as.numeric(t(obs[,-1]))
    
##########################################################################################
### Generate the mesh - same for all models ###
    bnd.int = inla.nonconvex.hull(loc, convex =-0.05)
    bnd.ext = inla.nonconvex.hull(loc, convex =-1)
    mesh <- inla.mesh.2d(loc.domain=loc, boundary = list(bnd.int,bnd.ext),  offset=c(-0.05,-0.25), max.edge=c(0.4,5))
    
    tmp = rbind(loc[1,], mesh$loc[,1:2])
    dist2s0.knots = as.matrix(dist(tmp))[1,-1]

##########################################################################################
### Select model form and run code
  # Model number from 0-6
    model.num <- 0
    source("Spatial_models.R")
  # Will save the output of the model as defined at the end of the 'Spatial_models.R' file.