#################################################################################
### Load required packages ###
library(INLA)

#################################################################################
### Code for all spatial models ###

### Model 0 ###
# Z(s)-Z(s0) #
if(model.num==0){
  set.seed(1)
  ## Generate the projection matrix to project from latent variables to observation sites:
  #(where we have to take into account the "-W(s0)" term)
  A.obs=inla.spde.make.A(mesh,loc=loc[-1,],index=rep(1:(nsites-1),n.repl),repl=rep(1:n.repl,each=nsites-1),n.spde=mesh$n)
  A.s0 = inla.spde.make.A(mesh, loc=loc[1,,drop=FALSE], index=rep(1,n.repl*(nsites-1)), repl=rep(1:n.repl,each=nsites-1), n.spde=mesh$n)
  dim(A.obs); dim(A.s0)
  A.obs = A.obs - A.s0
  
  #index over latent variables, with replication:
  idx.spatial=inla.spde.make.index("spatial",n.spde=mesh$n,n.repl=n.repl)
  
  #SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh,alpha=alpha,prior.range=c(1,.5),prior.sigma=c(1,.5))
  
  #covariate data frame for INLA (if any):
  covar.inla=data.frame(y.cond=rep(y.cond,each=nsites-1))

  #stack with all the information:
  mystack=inla.stack(data=list(y=yinla),A=list(A.obs,1),effects=list(idx.spatial,covar.inla))
  
  #model formula:
  myform=y~-1+f(spatial,model=myspde,replicate=spatial.repl,nrep=n.repl)
}

##########################################################################################
### Model 1 ###
# x0*alpha(s-s0) + Z(s)-Z(s0) #
if(model.num==1){
  set.seed(1)
  ## Generate the projection matrix from latent variables to observation sites:
  #(again, we have to take into account the "-W(s0)" term)
  A.obs=inla.spde.make.A(mesh,loc=loc[-1,],index=rep(1:(nsites-1),n.repl),repl=rep(1:n.repl,each=nsites-1),n.spde=mesh$n)
  A.s0 = inla.spde.make.A(mesh, loc=loc[1,,drop=FALSE], index=rep(1,n.repl*(nsites-1)), repl=rep(1:n.repl,each=nsites-1), n.spde=mesh$n)
  #dim(A.obs); dim(A.s0)
  A.obs = A.obs - A.s0
  
  #index over latent variables, with replication:
  idx.spatial=inla.spde.make.index("spatial",n.spde=mesh$n,n.repl=n.repl)
  
  #SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh,alpha=alpha,prior.range=c(1,.5),prior.sigma=c(1,.5))
  
  #set up a piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  #at 0, we fix the value to 0 (Dirichlet boundary conditions):
  mesh.dist=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("dirichlet","free")) 
  #define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1,NA),prior.sigma=c(.5,NA))
  idx.spline=inla.spde.make.index("spline1d", n.spde=myspde.dist$n.spde, n.group=1)
  A.dist=inla.spde.make.A(mesh=mesh.dist,loc=rep(dist2s0[-1],n.repl))
  #multiply with conditioning value X(s0):
  A.dist=A.dist * rep(y.cond,each=nsites-1)
  
  #covariate data frame for INLA (if any):
  covar.inla=data.frame(y.cond=rep(y.cond,each=nsites-1))
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=yinla),A=list(A.obs,A.dist,1),effects=list(idx.spatial,idx.spline,covar.inla))
  
  #model formula:
  myform=y~-1+f(spatial,model=myspde,replicate=spatial.repl,nrep=n.repl)+f(spline1d,model=myspde.dist)
}

##########################################################################################
### Model 2 ###
# x0 + gamma(s-s0) + Z(s)-Z(s0) #
if(model.num==2){
  set.seed(1)
  ## Generate the projection matrix from latent variables to observation sites:
  #(again, we have to take into account the "-W(s0)" term)
  A.obs=inla.spde.make.A(mesh,loc=loc[-1,],index=rep(1:(nsites-1),n.repl),repl=rep(1:n.repl,each=nsites-1),n.spde=mesh$n)
  A.s0 = inla.spde.make.A(mesh, loc=loc[1,,drop=FALSE], index=rep(1,n.repl*(nsites-1)), repl=rep(1:n.repl,each=nsites-1), n.spde=mesh$n)
  #dim(A.obs); dim(A.s0)
  A.obs = A.obs - A.s0
  
  #index over latent variables, with replication:
  idx.spatial=inla.spde.make.index("spatial",n.spde=mesh$n,n.repl=n.repl)
  
  #SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh,alpha=alpha, prior.range=c(1,.5), prior.sigma=c(1,.5))
  
  #set up second piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  #at 0, we fix the value to 0 (Dirichlet boundary conditions):
  mesh.dist2 = inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("dirichlet","free")) 
  #define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist2=inla.spde2.pcmatern(mesh=mesh.dist2,alpha=1.5,prior.range=c(1,NA),prior.sigma=c(.5,NA))
  idx.spline2=inla.spde.make.index("spline1d2", n.spde=myspde.dist2$n.spde, n.group=1)
  A.dist2=inla.spde.make.A(mesh=mesh.dist2,loc=rep(dist2s0[-1],n.repl))
  
  #covariate data frame for INLA (if any):
  covar.inla=data.frame(y.cond=rep(y.cond,each=nsites-1))
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=yinla),A=list(A.obs,A.dist2,1),effects=list(idx.spatial,idx.spline2,covar.inla))
  
  #model formula:
  myform=y~-1+f(spatial,model=myspde,replicate=spatial.repl,nrep=n.repl)+f(spline1d2,model=myspde.dist2)
}
  
##########################################################################################
### Model 3 ###
# x0*alpha(s-s0) + gamma(s-s0) + Z(s)-Z(s0) #
if(model.num==3){
  set.seed(1)
  ## Generate the projection matrix from latent variables to observation sites:
  #(again, we have to take into account the "-W(s0)" term)
  A.obs=inla.spde.make.A(mesh,loc=loc[-1,],index=rep(1:(nsites-1),n.repl),repl=rep(1:n.repl,each=nsites-1),n.spde=mesh$n)
  A.s0 = inla.spde.make.A(mesh, loc=loc[1,,drop=FALSE], index=rep(1,n.repl*(nsites-1)), repl=rep(1:n.repl,each=nsites-1), n.spde=mesh$n)
  dim(A.obs); dim(A.s0)
  A.obs = A.obs - A.s0
  
  #index over latent variables, with replication:
  idx.spatial=inla.spde.make.index("spatial",n.spde=mesh$n,n.repl=n.repl)
  
  #SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh,alpha=alpha,prior.range=c(1,.5),prior.sigma=c(1,.5))
  
  #set up a piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  #at 0, we fix the value to 0 (Dirichlet boundary conditions):
  mesh.dist=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("dirichlet","free"))
  #define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1,NA),prior.sigma=c(.5,NA))
  idx.spline=inla.spde.make.index("spline1d", n.spde=myspde.dist$n.spde, n.group=1)
  A.dist=inla.spde.make.A(mesh=mesh.dist,loc=rep(dist2s0[-1],n.repl))
  #multiply with conditioning value X(s0):
  A.dist=A.dist * rep(y.cond,each=nsites-1)
  
  #set up second piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  #at 0, we fix the value to 0 (Dirichlet boundary conditions):
  mesh.dist2=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("dirichlet","free"))
  #define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist2=inla.spde2.pcmatern(mesh=mesh.dist2,alpha=1.5,prior.range=c(1,NA),prior.sigma=c(.5,NA))
  idx.spline2=inla.spde.make.index("spline1d2", n.spde=myspde.dist2$n.spde, n.group=1)
  A.dist2=inla.spde.make.A(mesh=mesh.dist2,loc=rep(dist2s0[-1],n.repl))
  
  #covariate data frame for INLA (if any):
  covar.inla=data.frame(y.cond=rep(y.cond,each=nsites-1))
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=yinla),A=list(A.obs,A.dist,A.dist2,1),effects=list(idx.spatial,idx.spline,idx.spline2,covar.inla))
  
  #model formula:
  myform=y~-1+f(spatial,model=myspde,replicate=spatial.repl,nrep=n.repl)+f(spline1d,model=myspde.dist)+f(spline1d2,model=myspde.dist2)
}  

#################################################################################
### Model 4  ###
# x0*alpha(s-s0) + gamma(s-s0) + (x0^beta)*{Z(s)-Z(s0)} #
if(model.num==4){
  set.seed(1)
  library(Matrix)
  
  # Implement replicated spatial random effect with parametric b-function b(x)=x^beta.
  # tau (SPDE), kappa (SPDE), and beta are hyperparameters to be estimated. 
  # Within this model, we use covariates additional information:
  # xcond, nrep, mymesh, alpha.spde, prior.sigma (u, alpha), prior.range (u, alpha)
  rgeneric.beta = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = parent.env(environment())
    
    interpret.theta = function() {
      return(list(sigma=exp(theta[1L]), rho=exp(theta[2L]), beta=exp((theta[3L]))))
    }
    
    mu = function() {
      #par = interpret.theta() 
      #x0vec = rep(xcond, each = mymesh$n)
      #return(x0vec * exp(-(rep(dist2s0.knots, nrep)/par$range)^par$power))
      return(numeric(0))
    }
    
    # precision matrix
    Q = function() {
      par = interpret.theta()
      Q = inla.spde2.precision(myspde, theta = c(log(par$rho), log(par$sigma)))
      I = Diagonal(n = nrep, x = 1)
      Q = kronecker(I, Q)
      #fact = (1+mu( )^0.5)^{-1} # minus since we work with the precision matrix
      x0vec = rep(xcond, each = mymesh$n)
      fact = x0vec^{-par$beta}
      ##print(length(fact))
      len = length(fact)
      mydiag = spMatrix(len, len, i = 1:len, j = 1:len, x = fact)
      Q =  mydiag %*% Q 
      #print(class(Q))
      Q = Q %*% mydiag 
      #print(class(Q))
      return(Q)
    }
    
    
    # graph structure (binary matrix)
    graph = function() {
      return (Q())
    }
    
    # numeric(0) to tell it to calculate it from Q()-function
    log.norm.const = function() { 
      return(numeric(0))
    }
    
    # Put Gaussian prior with mean -log(2) and standard deviation 1 on log(beta).
    # Implement PC prior for sigma and rho
    log.prior = function() {
      val = 0
      lambda0 = -log(prior.sigma[2])/prior.sigma[1]
      ## Prior for standard deviation
      val = val + log(lambda0) - lambda0*exp(theta[1]) + theta[1]
      lambda1 = -log(prior.range[2])*prior.range[1]
      ## Prior for range
      val = val + log(lambda1) - lambda1*exp(-theta[2]) -theta[2]
      val = val + dnorm(theta[3], mean=-log(2), sd=1, log=TRUE) # beta
      return(val)
    }
    
    # initial values for theta
    initial = function() { 
      return(c(-0.5, 0, -log(2))) 
    }
    
    quit = function() { 
      return(invisible())
    }
    
    if (is.null(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    
    return(val) 
  }
  
  tmp = rbind(loc[1,], mesh$loc[,1:2])
  dist2s0.knots = as.matrix(dist(tmp))[1,-1]
  
  # Generate the projection matrix from latent variables to observation sites:
  A.obs=inla.spde.make.A(mesh,loc=loc[-1,],index=rep(1:(nsites-1),n.repl),repl=rep(1:n.repl,each=nsites-1),n.spde=mesh$n)
  A.s0 = inla.spde.make.A(mesh, loc=loc[1,,drop=FALSE], index=rep(1,n.repl*(nsites-1)), repl=rep(1:n.repl,each=nsites-1), n.spde=mesh$n)
  #dim(A.obs); dim(A.s0)
  A.obs = A.obs - A.s0
  
  #index over latent variables, with replication:
  spatial=1:(mesh$n*n.repl)
  
  #SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh,alpha=alpha, prior.range=c(1,.5), prior.sigma=c(1,.5))
  
  #set up a piecewise linear spline model (=1D SPDE):
  #(at 0, we fix the value to 0 corresponding to Dirichlet boundary conditions)
  mesh.dist=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("dirichlet","free")) 
  #define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1,NA),prior.sigma=c(.5,NA))
  idx.spline=inla.spde.make.index("spline1d", n.spde=myspde.dist$n.spde, n.group=1)
  A.dist=inla.spde.make.A(mesh=mesh.dist,loc=rep(dist2s0[-1],n.repl))
  #multiply with conditioning value X(s0):
  A.dist=A.dist * rep(y.cond,each=nsites-1)
  
  #set up second piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  #at 0, we fix the value to 0 (Dirichlet boundary conditions):
  mesh.dist2=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("dirichlet","free"))
  #define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist2=inla.spde2.pcmatern(mesh=mesh.dist2,alpha=1.5,prior.range=c(1,NA),prior.sigma=c(.5,NA))
  idx.spline2=inla.spde.make.index("spline1d2", n.spde=myspde.dist2$n.spde, n.group=1)
  A.dist2=inla.spde.make.A(mesh=mesh.dist2,loc=rep(dist2s0[-1],n.repl))
  
  #covariate data frame for INLA (if any):
  covar.inla=data.frame(y.cond=rep(y.cond,each=nsites-1))
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=yinla), A = list(A.obs, A.dist, A.dist2, 1), effects=list(list(spatial=spatial), idx.spline,idx.spline2,covar.inla))
  
  #model formula:
  prior.range=c(1,.5)
  prior.sigma=c(1,.5)
  mymodel = inla.rgeneric.define(model = rgeneric.beta, mymesh=mesh, dist2s0.knots=dist2s0.knots, alpha.spde=alpha, nrep=n.repl, xcond = y.cond, myspde = myspde, prior.sigma = prior.sigma, prior.range = prior.range)
  
  myform=y~-1 + f(spatial, model = mymodel) +
                f(spline1d,model=myspde.dist) + 
                f(spline1d2,model=myspde.dist2)
}

#################################################################################
### Model 5 ###
# x0*alpha(s-s0) + (x0^beta)*{Z(s)-Z(s0)} #
if(model.num==5){
  set.seed(1)
  library(Matrix)
  
  # Implement replicated spatial random effect with parametric b-function b(x)=x^beta.
  # tau (SPDE), kappa (SPDE), and beta are hyperparameters to be estimated. 
  # Within this model, we use covariates additional information:
  # xcond, nrep, mymesh, alpha.spde, prior.sigma (u, alpha), prior.range (u, alpha)
  rgeneric.beta = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = parent.env(environment())
    
    interpret.theta = function() {
      return(list(sigma=exp(theta[1L]), rho=exp(theta[2L]), beta=exp((theta[3L]))))
    }
    
    mu = function() {
      #par = interpret.theta() 
      #x0vec = rep(xcond, each = mymesh$n)
      #return(x0vec * exp(-(rep(dist2s0.knots, nrep)/par$range)^par$power))
      return(numeric(0))
    }
    
    # precision matrix
    Q = function() {
      par = interpret.theta()
      Q = inla.spde2.precision(myspde, theta = c(log(par$rho), log(par$sigma)))
      I = Diagonal(n = nrep, x = 1)
      Q = kronecker(I, Q)
      #fact = (1+mu( )^0.5)^{-1} # minus since we work with the precision matrix
      x0vec = rep(xcond, each = mymesh$n)
      fact = x0vec^{-par$beta}
      ##print(length(fact))
      len = length(fact)
      mydiag = spMatrix(len, len, i = 1:len, j = 1:len, x = fact)
      Q =  mydiag %*% Q 
      #print(class(Q))
      Q = Q %*% mydiag 
      #print(class(Q))
      return(Q)
    }
    
    
    # graph structure (binary matrix)
    graph = function() {
      return (Q())
    }
    
    # numeric(0) to tell it to calculate it from Q()-function
    log.norm.const = function() { 
      return(numeric(0))
    }
    
    # Put Gaussian prior with mean -log(2) and standard deviation 1 on log(beta).
    # Implement PC prior for sigma and rho
    log.prior = function() {
      val = 0
      lambda0 = -log(prior.sigma[2])/prior.sigma[1]
      ## Prior for standard deviation
      val = val + log(lambda0) - lambda0*exp(theta[1]) + theta[1]
      lambda1 = -log(prior.range[2])*prior.range[1]
      ## Prior for range
      val = val + log(lambda1) - lambda1*exp(-theta[2]) -theta[2]
      val = val + dnorm(theta[3], mean=-log(2), sd=1, log=TRUE) # beta
      return(val)
    }
    
    # initial values for theta
    initial = function() { 
      return(c(-0.5, 0, -log(2))) 
    }
    
    quit = function() { 
      return(invisible())
    }
    
    if (is.null(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    
    return(val) 
  }
  
  tmp = rbind(loc[1,], mesh$loc[,1:2])
  dist2s0.knots = as.matrix(dist(tmp))[1,-1]
  
  # Generate the projection matrix from latent variables to observation sites:
  A.obs=inla.spde.make.A(mesh,loc=loc[-1,],index=rep(1:(nsites-1),n.repl),repl=rep(1:n.repl,each=nsites-1),n.spde=mesh$n)
  A.s0 = inla.spde.make.A(mesh, loc=loc[1,,drop=FALSE], index=rep(1,n.repl*(nsites-1)), repl=rep(1:n.repl,each=nsites-1), n.spde=mesh$n)
  #dim(A.obs); dim(A.s0)
  A.obs = A.obs - A.s0
  
  #index over latent variables, with replication:
  spatial=1:(mesh$n*n.repl)
  
  #SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh,alpha=alpha, prior.range=c(1,.5), prior.sigma=c(1,.5))
  
  #set up a piecewise linear spline model (=1D SPDE):
  #(at 0, we fix the value to 0 corresponding to Dirichlet boundary conditions)
  mesh.dist=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("dirichlet","free")) 
  #define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1,NA),prior.sigma=c(.5,NA))
  idx.spline=inla.spde.make.index("spline1d", n.spde=myspde.dist$n.spde, n.group=1)
  A.dist=inla.spde.make.A(mesh=mesh.dist,loc=rep(dist2s0[-1],n.repl))
  #multiply with conditioning value X(s0):
  A.dist=A.dist * rep(y.cond,each=nsites-1)
  
  #covariate data frame for INLA (if any):
  covar.inla=data.frame(y.cond=rep(y.cond,each=nsites-1))
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=yinla), A = list(A.obs, A.dist, 1), effects=list(list(spatial=spatial), idx.spline,covar.inla))
  
  #model formula:
  prior.range=c(1,.5)
  prior.sigma=c(1,.5)
  mymodel = inla.rgeneric.define(model = rgeneric.beta, mymesh=mesh, dist2s0.knots=dist2s0.knots, alpha.spde=alpha, nrep=n.repl, xcond = y.cond, myspde = myspde, prior.range=prior.range, prior.sigma=prior.sigma)
  
  myform=y~-1 + f(spatial, model = mymodel) + f(spline1d,model=myspde.dist) 
}

#################################################################################
### Model 6 ###
# x0*alpha(s-s0) + gamma(s-s0) #
if(model.num==6){
  set.seed(1)
  
  #set up a piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  #at 0, we fix the value to 0 (Dirichlet boundary conditions):
  mesh.dist=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("dirichlet","free")) 
  #define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1,NA),prior.sigma=c(.5,NA))
  idx.spline=inla.spde.make.index("spline1d", n.spde=myspde.dist$n.spde, n.group=1)
  A.dist=inla.spde.make.A(mesh=mesh.dist,loc=rep(dist2s0[-1],n.repl))
  #multiply with conditioning value X(s0):
  A.dist=A.dist * rep(y.cond,each=nsites-1)
  
  #set up second piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  #at 0, we fix the value to 0 (Dirichlet boundary conditions):
  mesh.dist2=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("dirichlet","free"))
  #define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist2=inla.spde2.pcmatern(mesh=mesh.dist2,alpha=1.5,prior.range=c(1,NA),prior.sigma=c(.5,NA))
  idx.spline2=inla.spde.make.index("spline1d2", n.spde=myspde.dist2$n.spde, n.group=1)
  A.dist2=inla.spde.make.A(mesh=mesh.dist2,loc=rep(dist2s0[-1],n.repl))
  
  #covariate data frame for INLA (if any):
  covar.inla=data.frame(y.cond=rep(y.cond,each=nsites-1))
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=yinla),A=list(A.dist,A.dist2,1),effects=list(idx.spline,idx.spline2,covar.inla))
  
  #model formula:
  myform=y~-1+f(spline1d,model=myspde.dist)+f(spline1d2,model=myspde.dist2)
}


##########################################################################################
##########################################################################################
### Run INLA - same for all models ###
fit=inla(myform, 
         data=inla.stack.data(mystack), 
         family="gaussian",
         offset=covar.inla$y.cond,
         control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = FALSE),
         control.family = list(
           hyper = list(theta = list(initial = log(100), prior = "pc.prec", param = c(.1, .5)))),
         control.predictor=list(A=inla.stack.A(mystack), compute=FALSE, link=1),
         control.inla=list(strategy="simplified.laplace", int.strategy="eb"),
         control.fixed=list(prec=.1),
         num.thr = 2,
         verbose=FALSE)

##########################################################################################
### Select necessary output ###
  hyper = fit$summary.hyper
  fitted = fit$summary.fitted.values
  summ = summary(fit)
  spatial = fit$summary.random$spatial
  spline = fit$summary.random$spline1d
  spline2 = fit$summary.random$spline1d2
  obs = mystack$data$data$y
  pit = fit$cpo$pit
  cpo = fit$cpo$cpo


### Save results ###
  if(model.num == 0){
    save(model.num, hyper, fitted, summ, spatial, spline, obs, mesh, cpo, pit, file = paste0("AllSpatialModels-", "OUTPUT", model.num, "-alpha", alpha,  ".RData"))
  }else if(model.num==1){
    save(model.num, hyper, fitted, summ, spatial, spline, obs, mesh, mesh.dist, cpo, pit, file = paste0("AllSpatialModels-", "OUTPUT", model.num,  "-alpha", alpha,  ".RData"))
  }else if(model.num==2){
    save(model.num, hyper, fitted, summ, spatial, spline, spline2, obs, mesh, mesh.dist2, cpo, pit, file = paste0("AllSpatialModels-", "OUTPUT", model.num,  "-alpha", alpha,  ".RData"))
  }else if(model.num==3){
    save(model.num, hyper, fitted, summ, spatial, spline, spline2, obs, mesh, mesh.dist, mesh.dist2, cpo, pit, file = paste0("AllSpatialModels-", "OUTPUT", model.num,  "-alpha", alpha, ".RData"))
  }else if(model.num==4){
    save(model.num, hyper, fitted, summ, spatial, spline, spline2, obs, mesh, mesh.dist, mesh.dist2, cpo, pit, file = paste0("AllSpatialModels-", "OUTPUT", model.num, "-alpha", alpha,".RData"))
  }else if(model.num==5){
    save(model.num, hyper, fitted, summ, spatial, spline, obs, mesh, mesh.dist, cpo, pit, file = paste0("AllSpatialModels-", "OUTPUT", model.num, "-alpha", alpha,".RData"))
  }else if(model.num==6){
    save(model.num, hyper, fitted, summ, spline, spline2, obs,  mesh.dist, mesh.dist, cpo, pit, file = paste0("AllSpatialModels-", "OUTPUT", model.num, "-alpha", alpha, ".RData"))
  }
##########################################################################################
