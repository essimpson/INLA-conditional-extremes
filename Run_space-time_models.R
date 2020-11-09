#################################################################################
### First run file 'Space-time_declustering.R' to obtain independent clusters of observations

#################################################################################
### Fix necessary model parameters:
  # Choose number of latent variables in alpha and gamma spline functions
    knot.num <- 14 # corresponds to 15 knots including boundary knots (with spline constrained to 0 at 0)
  # Choose degree in spline functions ###
    spline.deg <- 2 # quadratic B-splines
  # Choose SPDE parameter alpha = nu + D/2
    alpha <- 2

#################################################################################
### Set up data for use with R-INLA
  # Summaries of data used within space-time code:
    n.repl <- sum(data[,6])
    ntimes <- max(data[,4]) 
    nsites <- max(data[,3])
    obs <- data[,1]
    
  # Use first point as conditioning location:
    dist2s0=as.matrix(dist(loc))[1,]
    y.cond <- data[data[,6]==1,1]
    
##########################################################################################
### Generate the mesh - same for all models 
    bnd.int = inla.nonconvex.hull(loc, convex =-0.05)
    bnd.ext = inla.nonconvex.hull(loc, convex =-1)
    mesh = inla.mesh.2d(loc.domain=loc, boundary = list(bnd.int,bnd.ext),  offset=c(-0.05,-0.25), max.edge=c(0.4,5))
    
  # Need to define some elements even if they are not needed (since they are "saved" at the end)
    mesh.dist = mesh.dist2 =  NULL
    
##########################################################################################
### Select model form and run code
  # Model number (1,3,4,5)
    model.num <- 1
    source("Space-time_models.R")
  # Will save the output of the model as defined at the end of the 'Space-time_models.R' file.