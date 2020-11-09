#################################################################################
### Load library and data ###
library(INLA)

#################################################################################
### Code for the four space-time models ###
### Model 1 ###
# x0*S(dist,t)+(Z(w)-Z(w_0)) with constrained ST quadratic spline (generic+A+extraconstr) 
if(model.num==1){
  
  # Implement quadratic spline (over space), with AR1-time-effect
  # hyperparameters: sigma, range (spline S), rho_ar
  # Within this model, we use additional information:
  # mymesh, myspde, alpha.spde, ntimes
  # so far, no prior.sigma (u, alpha), prior.range (u, alpha), prior.rho
  rgeneric.st_alpha = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = parent.env(environment())
    
    interpret.theta = function() {
      return()
      #return(list(sigma=exp(theta[1L]), range=exp(theta[2L]), rho_ar = exp(theta[3])/(1+exp(theta[3]))))
    }
    
    # precision matrix of AR1-model 
    # (adapted from inla.rgeneric.ar1.model)
    Q_AR1 = function(rho_ar, ntimes) {
      param = list(rho=rho_ar, prec = 1)
      i = c(1, ntimes, 2:(ntimes - 1), 1:(ntimes - 1))
      j = c(1, ntimes, 2:(ntimes - 1), 2:ntimes)
      x = param$prec/(1 - param$rho^2) * c(1, 1, rep(1 + param$rho^2, ntimes - 2), rep(-param$rho, ntimes - 1))
      Q = sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE)
      return(Q)
    }
    
    mu = function() {
      return(numeric(0))
    }
    
    # precision matrix
    Q = function() {
      par = list(sigma = 0.5, range = 1, rho_ar = .75)
      
      # par = interpret.theta()
     
      Q_S = inla.spde2.precision(myspde, theta = c(log(par$range), log(par$sigma)))
      Q_T = Q_AR1(par$rho_ar, ntimes)
      Q = kronecker(Q_T, Q_S)
      #fact = rep(xcond, each = mymesh$n*ntimes)
      #len = length(fact)  
      #mydiag = spMatrix(len, len, i = 1:len, j = 1:len, x = fact)
      #Q =  mydiag %*% Q 
      #Q = Q %*% mydiag 
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
      
      #lambda0 = -log(prior.sigma[2])/prior.sigma[1]
      #### Prior for standard deviation
      #val = val + log(lambda0) - lambda0*exp(theta[1]) + theta[1]
      #lambda1 = -log(prior.range[2])*prior.range[1]
      #### Prior for range
      #val = val + log(lambda1) - lambda1*exp(-theta[2]) - theta[2]
      #val = val + inla.pc.dcor1(exp(theta[3])/(1+exp(theta[3])), prior.ar[1], prior.ar[2], lambda, log = TRUE)
      
      return(val)
    }
    
    # initial values for theta
    initial = function() { 
      return()
      return(c(-0.5, 0, 0.5))#, -log(10))) 
    }
    
    quit = function() { 
      return(invisible())
    }
    
    if (is.null(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    
    return(val) 
  }
  
  set.seed(1)
  #generate the projection matrix from latent variables to observation sites:
  #(again, we have to take into account the "-W(s0)" term)
  A.obs=inla.spde.make.A(mesh, loc=loc, index = data[,3], repl=data[,5], group = data[,4], n.spde=mesh$n) 
  A.s0=inla.spde.make.A(mesh, loc=loc[1, , drop=FALSE], index = 1, repl=1, group=1,  n.spde=mesh$n, n.group = ntimes, n.repl = n.repl) 
  dim(A.obs); dim(A.s0)
  #which(A.s0[1,]>0)
  idx2mod=which(A.s0[1,] != 0)
  idx2mod
  for(i in 1:n.repl){
    print(i)
    for(j in 1:length(idx2mod)){
      A.obs[(i-1)*nsites*ntimes+1:(nsites*ntimes), idx2mod[j]] = A.obs[(i-1)*nsites*ntimes + 1:(nsites*ntimes), idx2mod[j]] - A.s0[1, idx2mod[j]] 
    }
  }
  
  y.cond.repl = rep(y.cond, each = nsites*ntimes)
  
  #index over latent variables, with replication:
  idx.spatial=inla.spde.make.index("spatial", n.spde=mesh$n, n.repl=n.repl, n.group = ntimes) 
  
  ##SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh,alpha=alpha,prior.range=c(1,.5),prior.sigma=c(1,.5))
  
  #set up a piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  mesh.dist=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("free","free")) 
  # define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1, NA),prior.sigma=c(.5, NA))
  #idx.spline=inla.spde.make.index("spline1d", n.spde=myspde.dist$n.spde, n.group = ntimes, n.repl = 1)
  A.dist=inla.spde.make.A(mesh=mesh.dist, loc=rep(dist2s0, n.repl * ntimes), n.repl = 1, group = rep(rep(1:ntimes, each = nsites), n.repl)) 
  #A.dist.cond=inla.spde.make.A(mesh.dist, loc=0, index = 1, repl=1, group=1,  n.spde=mesh.dist$n, n.group = ntimes, n.repl = 1) 
  #dim(A.dist); dim(A.dist.cond)
  #idx2mod=which(A.dist.cond[1,] != 0)
  #idx2mod
  #for(j in 1:length(idx2mod)){
  #  A.dist[, idx2mod[j]] = A.dist[ , idx2mod[j]] - A.dist.cond[1, idx2mod[j]] 
  #}
  # multiply with conditioning value X(s0):
  A.dist = A.dist * y.cond.repl 
  dim(A.dist)
  
  # covariate data frame for INLA (if any) (not used currently):
  covar.inla=data.frame(y.cond = rep(y.cond, each=(nsites) * ntimes)) 
  dim(covar.inla)
  
  # we remove the data for the conditioning point (which may help avoid instabilities):
  idx.cond = which(data[,6]==1)
  covar.inla = covar.inla[-idx.cond, , drop = FALSE]
  A.obs = A.obs[-idx.cond, ]
  A.dist = A.dist[-idx.cond, ]
  tmp = data[,1]
  y.inla = as.numeric(tmp)[-idx.cond]
  #y.inla = as.numeric(data[,1])[-idx.cond]
  
  # generate rgeneric model:
  #model formula:
  
  #prior.range=c(1,.5)
  #prior.sigma=c(1,.5)
  #prior.ar = c(0.8, 0.5)
  
  mymodel = inla.rgeneric.define(model = rgeneric.st_alpha, mymesh=mesh.dist, alpha.spde=1.5, myspde = myspde.dist, ntimes = ntimes) #, prior.sigma = prior.sigma, prior.range = prior.range, prior.ar = prior.ar)
  #index:
  idx.spline.st=1:(mesh.dist$m*ntimes)
  # create spline constraints:
  A.dist.cond=inla.spde.make.A(mesh.dist, loc=0, index = 1, repl=1, group=1,  n.spde=mesh.dist$n, n.group = 1, n.repl = 1) 
  A.extraconstr = t(c(.5,.5,rep(0,mesh.dist$m*ntimes-2)))
  e.extraconstr = 0
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=y.inla), A=list(A.obs, A.dist, 1), effects=list(idx.spatial, list(idx.spline.st=idx.spline.st), covar.inla))
  
  #model formula: ####
  hyper.ar1=list(theta=list(prior="pccor0", param=c(.5, .5))) # a priori, probability of having rho > 0.5 is 0.5
  myform=y~-1+
    f(spatial, model=myspde, replicate=spatial.repl, group = spatial.group, nrep=n.repl, ngroup = ntimes, control.group=list(model="ar1", hyper=hyper.ar1)) + 
    f(idx.spline.st, model = mymodel, extraconstr = list(A = A.extraconstr, e = e.extraconstr)) 
    #f(spline1d, model=myspde.dist, group = spline1d.group, control.group = list(model = "ar1", hyper = list(theta = list(initial = log((1+0.75)/(1-0.75)), fixed = TRUE)))) 
}

#################################################################################
### Model 3 ###
# x0*S(dist,t)+ S2(dist,t) + (Z(w)-Z(w_0)) with constrained ST quadratic spline (generic+A+extraconstr) 
if(model.num==3){
  
  # Implement quadratic spline (over space), with AR1-time-effect
  # hyperparameters: sigma, range (spline S), rho_ar
  # Within this model, we use additional information:
  # mymesh, myspde, alpha.spde, ntimes
  # so far, no prior.sigma (u, alpha), prior.range (u, alpha), prior.rho
  rgeneric.st_alpha = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = parent.env(environment())
    
    interpret.theta = function() {
      return()
      #return(list(sigma=exp(theta[1L]), range=exp(theta[2L]), rho_ar = exp(theta[3])))
    }
    
    # precision matrix of AR1-model 
    # (adapted from inla.rgeneric.ar1.model)
    Q_AR1 = function(rho_ar, ntimes) {
      param = list(rho=rho_ar, prec = 1)
      i = c(1, ntimes, 2:(ntimes - 1), 1:(ntimes - 1))
      j = c(1, ntimes, 2:(ntimes - 1), 2:ntimes)
      x = param$prec/(1 - param$rho^2) * c(1, 1, rep(1 + param$rho^2, ntimes - 2), rep(-param$rho, ntimes - 1))
      Q = sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE)
      return(Q)
    }
    
    mu = function() {
      return(numeric(0))
    }
    
    # precision matrix
    Q = function() {
      #par = interpret.theta()
      par = list(sigma = 0.5, range = 1, rho_ar = .75)
      Q_S = inla.spde2.precision(myspde, theta = c(log(par$range), log(par$sigma)))
      Q_T = Q_AR1(par$rho_ar, ntimes)
      Q = kronecker(Q_T, Q_S)
      #fact = rep(xcond, each = mymesh$n*ntimes)
      #len = length(fact)  
      #mydiag = spMatrix(len, len, i = 1:len, j = 1:len, x = fact)
      #Q =  mydiag %*% Q 
      #Q = Q %*% mydiag 
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
      #lambda0 = -log(prior.sigma[2])/prior.sigma[1]
      ### Prior for standard deviation
      #val = val + log(lambda0) - lambda0*exp(theta[1]) + theta[1]
      #lambda1 = -log(prior.range[2])*prior.range[1]
      ### Prior for range
      #val = val + log(lambda1) - lambda1*exp(-theta[2]) - theta[2]
      return(val)
    }
    
    # initial values for theta
    initial = function() { 
      return()
      #return(c(-0.5, 0, -log(10)))#, -log(10))) 
    }
    
    quit = function() { 
      return(invisible())
    }
    
    if (is.null(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    
    return(val) 
  }
  
  set.seed(1)
  #generate the projection matrix from latent variables to observation sites:
  #(again, we have to take into account the "-W(s0)" term)
  A.obs=inla.spde.make.A(mesh, loc=loc, index = data[,3], repl=data[,5], group = data[,4], n.spde=mesh$n) 
  A.s0=inla.spde.make.A(mesh, loc=loc[1, , drop=FALSE], index = 1, repl=1, group=1,  n.spde=mesh$n, n.group = ntimes, n.repl = n.repl) 
  dim(A.obs); dim(A.s0)
  #which(A.s0[1,]>0)
  idx2mod=which(A.s0[1,] != 0)
  idx2mod
  for(i in 1:n.repl){
    print(i)
    for(j in 1:length(idx2mod)){
      A.obs[(i-1)*nsites*ntimes+1:(nsites*ntimes), idx2mod[j]] = A.obs[(i-1)*nsites*ntimes + 1:(nsites*ntimes), idx2mod[j]] - A.s0[1, idx2mod[j]] 
    }
  }
  
  y.cond.repl = rep(y.cond, each = nsites*ntimes)
  
  #index over latent variables, with replication:
  idx.spatial=inla.spde.make.index("spatial", n.spde=mesh$n, n.repl=n.repl, n.group = ntimes) 
  
  ##SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh,alpha=alpha,prior.range=c(1,.5),prior.sigma=c(1,.5))
  
  #set up a piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  mesh.dist=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("free","free")) 
  mesh.dist2 = inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("free","free")) 
  # define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1, NA),prior.sigma=c(.5, NA))
  myspde.dist2=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1, NA),prior.sigma=c(.5, NA))
  #idx.spline=inla.spde.make.index("spline1d", n.spde=myspde.dist$n.spde, n.group = ntimes, n.repl = 1)
  A.dist=inla.spde.make.A(mesh=mesh.dist, loc=rep(dist2s0, n.repl * ntimes), n.repl = 1, group = rep(rep(1:ntimes, each = nsites), n.repl)) 
  #A.dist.cond=inla.spde.make.A(mesh.dist, loc=0, index = 1, repl=1, group=1,  n.spde=mesh.dist$n, n.group = ntimes, n.repl = 1) 
  #dim(A.dist); dim(A.dist.cond)
  #idx2mod=which(A.dist.cond[1,] != 0)
  #idx2mod
  #for(j in 1:length(idx2mod)){
  #  A.dist[, idx2mod[j]] = A.dist[ , idx2mod[j]] - A.dist.cond[1, idx2mod[j]] 
  #}
  A.dist2 = A.dist
  # multiply with conditioning value X(s0):
  A.dist = A.dist * y.cond.repl 
  dim(A.dist)
  
  # covariate data frame for INLA (if any) (not used currently):
  covar.inla=data.frame(y.cond = rep(y.cond, each=(nsites) * ntimes)) 
  dim(covar.inla)
  
  # we remove the data for the conditioning point (which may help avoid instabilities):
  idx.cond = which(data[,6]==1)
  covar.inla = covar.inla[-idx.cond, , drop = FALSE]
  A.obs = A.obs[-idx.cond, ]
  A.dist = A.dist[-idx.cond, ]
  A.dist2 = A.dist2[-idx.cond, ]
  #y.inla = as.numeric(data[,1])[-idx.cond]
  tmp = data[,1]
  y.inla = as.numeric(tmp)[-idx.cond]
  
  # generate rgeneric model:
  #model formula:
  #prior.range=c(1,.5)
  #prior.sigma=c(1,.5)
  #prior.rho=...+
  mymodel = inla.rgeneric.define(model = rgeneric.st_alpha, mymesh=mesh.dist, alpha.spde=1.5, myspde = myspde.dist, ntimes = ntimes)#, prior.sigma = prior.sigma, prior.range = prior.range)
  mymodel2 = inla.rgeneric.define(model = rgeneric.st_alpha, mymesh=mesh.dist, alpha.spde=1.5, myspde = myspde.dist, ntimes = ntimes)
  #index:
  idx.spline.st=1:(mesh.dist$m*ntimes)
  idx.spline2.st=1:(mesh.dist$m*ntimes)
  # create spline constraints:
  A.dist.cond=inla.spde.make.A(mesh.dist, loc=0, index = 1, repl=1, group=1,  n.spde=mesh.dist$n, n.group = 1, n.repl = 1) 
  A.extraconstr = t(c(.5,.5,rep(0,mesh.dist$m*ntimes-2)))
  e.extraconstr = 0
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=y.inla), A=list(A.obs, A.dist, A.dist2, 1), effects=list(idx.spatial, list(idx.spline.st=idx.spline.st), list(idx.spline2.st=idx.spline.st), covar.inla))
  
  #model formula: ####
  hyper.ar1=list(theta=list(prior="pccor0", param=c(.5, .5))) # a priori, probability of having rho > 0.5 is 0.5
  myform=y~-1+
    f(spatial,model=myspde, replicate=spatial.repl, group = spatial.group, nrep=n.repl, ngroup = ntimes, control.group=list(model="ar1", hyper=hyper.ar1)) + 
    f(idx.spline.st, model = mymodel, extraconstr = list(A = A.extraconstr, e = e.extraconstr)) + 
    f(idx.spline2.st, model = mymodel, extraconstr = list(A = A.extraconstr, e = e.extraconstr)) 
}

#################################################################################
### Model 4 ###
# x0*S(dist,t)+ S2(dist,t) + x^beta*(Z(w)-Z(w_0)) with constrained ST quadratic spline (generic+A+extraconstr) 
if(model.num == 4){
  
  # Implement quadratic spline (over space), with AR1-time-effect
  # hyperparameters: sigma, range (spline S), rho_ar
  # Within this model, we use additional information:
  # mymesh, myspde, alpha.spde, ntimes
  # so far, no prior.sigma (u, alpha), prior.range (u, alpha), prior.rho (u, alpha)
  rgeneric.st_alpha = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = parent.env(environment())
    
    interpret.theta = function() {
      return()
      #return(list(sigma=exp(theta[1L]), range=exp(theta[2L]), rho_ar = exp(theta[3])))
    }
    
    # precision matrix of AR1-model 
    # (adapted from inla.rgeneric.ar1.model)
    Q_AR1 = function(rho_ar, ntimes){
      param = list(rho=rho_ar, prec = 1)
      i = c(1, ntimes, 2:(ntimes - 1), 1:(ntimes - 1))
      j = c(1, ntimes, 2:(ntimes - 1), 2:ntimes)
      x = param$prec/(1 - param$rho^2) * c(1, 1, rep(1 + param$rho^2, ntimes - 2), rep(-param$rho, ntimes - 1))
      Q = sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE)
      return(Q)
    }
    
    mu = function() {
      return(numeric(0))
    }
    
    # precision matrix
    Q = function() {
      #par = interpret.theta()
      par = list(sigma = 0.5, range = 1, rho_ar = .75)
      Q_S = inla.spde2.precision(myspde, theta = c(log(par$range), log(par$sigma)))
      Q_T = Q_AR1(par$rho_ar, ntimes)
      Q = kronecker(Q_T, Q_S)
      #fact = rep(xcond, each = mymesh$n*ntimes)
      #len = length(fact)  
      #mydiag = spMatrix(len, len, i = 1:len, j = 1:len, x = fact)
      #Q =  mydiag %*% Q 
      #Q = Q %*% mydiag 
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
      #lambda0 = -log(prior.sigma[2])/prior.sigma[1]
      ### Prior for standard deviation
      #val = val + log(lambda0) - lambda0*exp(theta[1]) + theta[1]
      #lambda1 = -log(prior.range[2])*prior.range[1]
      ### Prior for range
      #val = val + log(lambda1) - lambda1*exp(-theta[2]) - theta[2]
      return(val)
    }
    
    # initial values for theta
    initial = function() { 
      return()
      #return(c(-0.5, 0, -log(10)))#, -log(10))) 
    }
    
    quit = function() { 
      return(invisible())
    }
    
    if (is.null(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    
    return(val) 
  }
  
  # Implement replicated spatial random effect with parametric b-function b(x)=x^beta.
  # tau (SPDE), kappa (SPDE), and beta are hyperparameters to be estimated. 
  # Within this model, we use covariates additional information:
  # xcond, nrep, ntimes mymesh, alpha.spde, prior.sigma (u, alpha), prior.range (u, alpha)
  rgeneric.beta = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = parent.env(environment())
    
    interpret.theta = function() {
      return(list(sigma=exp(theta[1]), range=exp(theta[2]),  beta=exp(theta[3]), rho = exp(theta[4])/(1+exp(theta[4]))))
      # beta=exp(theta[3])/(1+exp(theta[3]))
    }
    
    # precision matrix of AR1-model 
    # (adapted from inla.rgeneric.ar1.model)
    Q_AR1 = function(rho_ar, ntimes) {
      param = list(rho=rho_ar, prec = 1)
      i = c(1, ntimes, 2:(ntimes - 1), 1:(ntimes - 1))
      j = c(1, ntimes, 2:(ntimes - 1), 2:ntimes)
      x = param$prec/(1 - param$rho^2) * c(1, 1, rep(1 + param$rho^2, ntimes - 2), rep(-param$rho, ntimes - 1))
      Q = sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE)
      return(Q)
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
      Q_S = inla.spde2.precision(myspde, theta = c(log(par$range), log(par$sigma)))
      Q_T = Q_AR1(par$rho, ntimes)
      Q = kronecker(Q_T, Q_S)
      I = Diagonal(n = nrep, x = 1)
      Q = kronecker(I, Q)
      #fact = (1+mu( )^0.5)^{-1} # minus since we work with the precision matrix
      x0vec = rep(xcond, each = mymesh$n*ntimes)
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
      val = val + inla.pc.dcor0(interpret.theta()$rho, prior.rho[1], prior.rho[2], log = TRUE)
      return(val)
    }
    
    # initial values for theta
    initial = function() { 
      return(c(-0.5, 0, -log(2), 0.5)) 
    }
    
    quit = function() { 
      return(invisible())
    }
    
    if (is.null(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    
    return(val) 
  }
  
  set.seed(1)
  A.obs=inla.spde.make.A(mesh, loc=loc, index = data[,3], repl=data[,5], group = data[,4], n.spde=mesh$n) 
  A.s0=inla.spde.make.A(mesh, loc=loc[1, , drop=FALSE], index = rep(1, n.repl), repl=1:n.repl, group=1, n.spde=mesh$n, n.group = ntimes, n.repl = n.repl) 
  dim(A.obs); dim(A.s0)
  for(i in 1:n.repl){
    print(i)
    idx2mod=which(A.s0[i,] != 0)
    print(idx2mod)
    for(j in 1:length(idx2mod)){
      A.obs[(i-1)*nsites*ntimes+1:(nsites*ntimes), idx2mod[j]] = A.obs[(i-1)*nsites*ntimes + 1:(nsites*ntimes), idx2mod[j]] - A.s0[i, idx2mod[j]] 
    }
  }

  y.cond.repl = rep(y.cond, each = nsites*ntimes)
  
  #index over latent variables, with replication:
  spatial=1:(mesh$n*n.repl*ntimes)
  
  ##SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh, alpha=alpha, prior.range=c(1,.5),prior.sigma=c(1,.5))
  
  #set up a piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  mesh.dist=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("free","free")) 
  mesh.dist2 = inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("free","free")) 
  # define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5, prior.range=c(1, NA),prior.sigma=c(.5, NA))
  myspde.dist2=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1, NA),prior.sigma=c(.5, NA))
  #idx.spline=inla.spde.make.index("spline1d", n.spde=myspde.dist$n.spde, n.group = ntimes, n.repl = 1)
  A.dist=inla.spde.make.A(mesh=mesh.dist, loc=rep(dist2s0, n.repl * ntimes), n.repl = 1, group = rep(rep(1:ntimes, each = nsites), n.repl)) 
  #A.dist.cond=inla.spde.make.A(mesh.dist, loc=0, index = 1, repl=1, group=1,  n.spde=mesh.dist$n, n.group = ntimes, n.repl = 1) 
  #dim(A.dist); dim(A.dist.cond)
  #idx2mod=which(A.dist.cond[1,] != 0)
  #idx2mod
  #for(j in 1:length(idx2mod)){
  #  A.dist[, idx2mod[j]] = A.dist[ , idx2mod[j]] - A.dist.cond[1, idx2mod[j]] 
  #}
  A.dist2 = A.dist
  # multiply with conditioning value X(s0):
  A.dist = A.dist * y.cond.repl 
  dim(A.dist)
  
  # covariate data frame for INLA (if any) (not used currently):
  covar.inla=data.frame(y.cond = rep(y.cond, each = nsites * ntimes)) 
  dim(covar.inla)
  
  # we remove the data for the conditioning point (which may help avoid instabilities):
  idx.cond = which(data[,6]==1)
  covar.inla = covar.inla[-idx.cond, , drop = FALSE]
  A.obs = A.obs[-idx.cond, ]
  A.dist = A.dist[-idx.cond, ]
  A.dist2 = A.dist2[-idx.cond, ]
  #y.inla = as.numeric(data[,1])[-idx.cond]
  tmp = data[,1]
  y.inla = as.numeric(tmp)[-idx.cond]
  
  # generate rgeneric model:
  #model formula:
  #prior.range=c(1,.5)
  #prior.sigma=c(1,.5)
  #prior.rho=...+
  mymodel = inla.rgeneric.define(model = rgeneric.st_alpha, mymesh=mesh.dist, alpha.spde=1.5, myspde = myspde.dist, ntimes = ntimes)#, prior.sigma = prior.sigma, prior.range = prior.range)
  mymodel2 = inla.rgeneric.define(model = rgeneric.st_alpha, mymesh=mesh.dist, alpha.spde=1.5, myspde = myspde.dist, ntimes = ntimes)
  #model formula:
  prior.range=c(1,.5)
  prior.sigma=c(1,.5)
  prior.rho=c(.5,.5)
  mymodel3 = inla.rgeneric.define(model = rgeneric.beta, mymesh=mesh, alpha.spde=alpha, nrep=n.repl, ntimes = ntimes, xcond = y.cond, myspde = myspde, prior.sigma = prior.sigma, prior.range = prior.range, prior.rho = prior.rho)
  #index:
  idx.spline.st=1:(mesh.dist$m*ntimes)
  idx.spline2.st=1:(mesh.dist2$m*ntimes)
  # create spline constraints:
  #A.dist.cond=inla.spde.make.A(mesh.dist, loc=0, index = 1, repl=1, group=1,  n.spde=mesh.dist$n, n.group = 1, n.repl = 1) 
  A.extraconstr = t(c(.5,.5,rep(0,mesh.dist$m*ntimes-2)))
  e.extraconstr = 0
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=y.inla), A=list(A.obs, A.dist, A.dist2, 1), effects=list(list(spatial=spatial), list(idx.spline.st=idx.spline.st), list(idx.spline2.st=idx.spline.st), covar.inla))
  
  #model formula: ####
  #hyper.ar1=list(theta=list(prior="pccor0", param=c(.5, .5))) # a priori, probability of having rho > 0.5 is 0.5
  myform=y~-1+
    f(idx.spline.st, model = mymodel, extraconstr = list(A = A.extraconstr, e = e.extraconstr)) + 
    f(idx.spline2.st, model = mymodel2, extraconstr = list(A = A.extraconstr, e = e.extraconstr)) +
    f(spatial, model = mymodel3) 
}

#################################################################################
### Model 5 ###
# x0*S(dist,t)+ x^beta*(Z(w)-Z(w_0)) with constrained ST quadratic spline (generic+A+extraconstr) 
if(model.num == 5){
  # Implement quadratic spline (over space), with AR1-time-effect
  # hyperparameters: sigma, range (spline S), rho_ar
  # Within this model, we use additional information:
  # mymesh, myspde, alpha.spde, ntimes
  # so far, no prior.sigma (u, alpha), prior.range (u, alpha), prior.rho (u, alpha)
  rgeneric.st_alpha = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = parent.env(environment())
    
    interpret.theta = function() {
      return()
      #return(list(sigma=exp(theta[1L]), range=exp(theta[2L]), rho_ar = exp(theta[3])))
    }
    
    # precision matrix of AR1-model 
    # (adapted from inla.rgeneric.ar1.model)
    Q_AR1 = function(rho_ar, ntimes){
      param = list(rho=rho_ar, prec = 1)
      i = c(1, ntimes, 2:(ntimes - 1), 1:(ntimes - 1))
      j = c(1, ntimes, 2:(ntimes - 1), 2:ntimes)
      x = param$prec/(1 - param$rho^2) * c(1, 1, rep(1 + param$rho^2, ntimes - 2), rep(-param$rho, ntimes - 1))
      Q = sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE)
      return(Q)
    }
    
    mu = function() {
      return(numeric(0))
    }
    
    # precision matrix
    Q = function() {
      #par = interpret.theta()
      par = list(sigma = 0.5, range = 1, rho_ar = .75)
      Q_S = inla.spde2.precision(myspde, theta = c(log(par$range), log(par$sigma)))
      Q_T = Q_AR1(par$rho_ar, ntimes)
      Q = kronecker(Q_T, Q_S)
      #fact = rep(xcond, each = mymesh$n*ntimes)
      #len = length(fact)  
      #mydiag = spMatrix(len, len, i = 1:len, j = 1:len, x = fact)
      #Q =  mydiag %*% Q 
      #Q = Q %*% mydiag 
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
      #lambda0 = -log(prior.sigma[2])/prior.sigma[1]
      ### Prior for standard deviation
      #val = val + log(lambda0) - lambda0*exp(theta[1]) + theta[1]
      #lambda1 = -log(prior.range[2])*prior.range[1]
      ### Prior for range
      #val = val + log(lambda1) - lambda1*exp(-theta[2]) - theta[2]
      return(val)
    }
    
    # initial values for theta
    initial = function() { 
      return()
      #return(c(-0.5, 0, -log(10)))#, -log(10))) 
    }
    
    quit = function() { 
      return(invisible())
    }
    
    if (is.null(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    
    return(val) 
  }
  
  # Implement replicated spatial random effect with parametric b-function b(x)=x^beta.
  # tau (SPDE), kappa (SPDE), and beta are hyperparameters to be estimated. 
  # Within this model, we use covariates additional information:
  # xcond, nrep, ntimes mymesh, alpha.spde, prior.sigma (u, alpha), prior.range (u, alpha)
  rgeneric.beta = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = parent.env(environment())
    
    interpret.theta = function() {
      return(list(sigma=exp(theta[1]), range=exp(theta[2]),  beta=exp(theta[3]), rho = exp(theta[4])/(1+exp(theta[4]))))
      # beta=exp(theta[3])/(1+exp(theta[3]))
    }
    
    # precision matrix of AR1-model 
    # (adapted from inla.rgeneric.ar1.model)
    Q_AR1 = function(rho_ar, ntimes) {
      param = list(rho=rho_ar, prec = 1)
      i = c(1, ntimes, 2:(ntimes - 1), 1:(ntimes - 1))
      j = c(1, ntimes, 2:(ntimes - 1), 2:ntimes)
      x = param$prec/(1 - param$rho^2) * c(1, 1, rep(1 + param$rho^2, ntimes - 2), rep(-param$rho, ntimes - 1))
      Q = sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE)
      return(Q)
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
      Q_S = inla.spde2.precision(myspde, theta = c(log(par$range), log(par$sigma)))
      Q_T = Q_AR1(par$rho, ntimes)
      Q = kronecker(Q_T, Q_S)
      I = Diagonal(n = nrep, x = 1)
      Q = kronecker(I, Q)
      #fact = (1+mu( )^0.5)^{-1} # minus since we work with the precision matrix
      x0vec = rep(xcond, each = mymesh$n*ntimes)
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
      val = val + inla.pc.dcor0(interpret.theta()$rho, prior.rho[1], prior.rho[2], log = TRUE)
      return(val)
    }
    
    # initial values for theta
    initial = function() { 
      return(c(-0.5, 0, -log(2), 0.5)) 
    }
    
    quit = function() { 
      return(invisible())
    }
    
    if (is.null(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    
    return(val) 
  }
  
  set.seed(1)
  A.obs=inla.spde.make.A(mesh, loc=loc, index = data[,3], repl=data[,5], group = data[,4], n.spde=mesh$n) 
  A.s0=inla.spde.make.A(mesh, loc=loc[1, , drop=FALSE], index = rep(1, n.repl), repl=1:n.repl, group=1, n.spde=mesh$n, n.group = ntimes, n.repl = n.repl) 
  dim(A.obs); dim(A.s0)
  for(i in 1:n.repl){
    print(i)
    idx2mod=which(A.s0[i,] != 0)
    print(idx2mod)
    for(j in 1:length(idx2mod)){
      A.obs[(i-1)*nsites*ntimes+1:(nsites*ntimes), idx2mod[j]] = A.obs[(i-1)*nsites*ntimes + 1:(nsites*ntimes), idx2mod[j]] - A.s0[i, idx2mod[j]] 
    }
  }
  
  y.cond.repl = rep(y.cond, each = nsites*ntimes)
  
  #index over latent variables, with replication:
  spatial=1:(mesh$n*n.repl*ntimes)
  
  ##SPDE model: Matern with PC prior
  myspde=inla.spde2.pcmatern(mesh, alpha=alpha, prior.range=c(1,.5),prior.sigma=c(1,.5))
  
  #set up a piecewise linear spline model (=1D SPDE) using the same mechanism as for Z_0(s):
  mesh.dist=inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("free","free")) 
  mesh.dist2 = inla.mesh.1d(loc=0:knot.num/(knot.num-3)*max(dist2s0),interval=c(0,knot.num/(knot.num-3)*max(dist2s0)),degree=spline.deg,boundary=c("free","free")) 
  # define the SPDE prior model (here with fixed range and variance parameters):
  myspde.dist=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5, prior.range=c(1, NA),prior.sigma=c(.5, NA))
  myspde.dist2=inla.spde2.pcmatern(mesh=mesh.dist,alpha=1.5,prior.range=c(1, NA),prior.sigma=c(.5, NA))
  #idx.spline=inla.spde.make.index("spline1d", n.spde=myspde.dist$n.spde, n.group = ntimes, n.repl = 1)
  A.dist=inla.spde.make.A(mesh=mesh.dist, loc=rep(dist2s0, n.repl * ntimes), n.repl = 1, group = rep(rep(1:ntimes, each = nsites), n.repl)) 
  #A.dist.cond=inla.spde.make.A(mesh.dist, loc=0, index = 1, repl=1, group=1,  n.spde=mesh.dist$n, n.group = ntimes, n.repl = 1) 
  #dim(A.dist); dim(A.dist.cond)
  #idx2mod=which(A.dist.cond[1,] != 0)
  #idx2mod
  #for(j in 1:length(idx2mod)){
  #  A.dist[, idx2mod[j]] = A.dist[ , idx2mod[j]] - A.dist.cond[1, idx2mod[j]] 
  #}
  A.dist2 = A.dist
  # multiply with conditioning value X(s0):
  A.dist = A.dist * y.cond.repl 
  dim(A.dist)
  
  # covariate data frame for INLA (if any) (not used currently):
  covar.inla=data.frame(y.cond = rep(y.cond, each = nsites * ntimes)) 
  dim(covar.inla)
  
  # we remove the data for the conditioning point (which may help avoid instabilities):
  idx.cond = which(data[,6]==1)
  covar.inla = covar.inla[-idx.cond, , drop = FALSE]
  A.obs = A.obs[-idx.cond, ]
  A.dist = A.dist[-idx.cond, ]
  A.dist2 = A.dist2[-idx.cond, ]
  #y.inla = as.numeric(data[,1])[-idx.cond]
  tmp = data[,1]
  y.inla = as.numeric(tmp)[-idx.cond]
  
  # generate rgeneric model:
  #model formula:
  #prior.range=c(1,.5)
  #prior.sigma=c(1,.5)
  #prior.rho=...+
  mymodel = inla.rgeneric.define(model = rgeneric.st_alpha, mymesh=mesh.dist, alpha.spde=1.5, myspde = myspde.dist, ntimes = ntimes)#, prior.sigma = prior.sigma, prior.range = prior.range)
  mymodel2 = inla.rgeneric.define(model = rgeneric.st_alpha, mymesh=mesh.dist, alpha.spde=1.5, myspde = myspde.dist, ntimes = ntimes)
  #model formula:
  prior.range=c(1,.5)
  prior.sigma=c(1,.5)
  prior.rho=c(.5,.5)
  mymodel3 = inla.rgeneric.define(model = rgeneric.beta, mymesh=mesh, alpha.spde=alpha, nrep=n.repl, ntimes = ntimes, xcond = y.cond, myspde = myspde, prior.sigma = prior.sigma, prior.range = prior.range, prior.rho = prior.rho)
  #index:
  idx.spline.st=1:(mesh.dist$m*ntimes)
  idx.spline2.st=1:(mesh.dist2$m*ntimes)
  # create spline constraints:
  #A.dist.cond=inla.spde.make.A(mesh.dist, loc=0, index = 1, repl=1, group=1,  n.spde=mesh.dist$n, n.group = 1, n.repl = 1) 
  A.extraconstr = t(c(.5,.5,rep(0,mesh.dist$m*ntimes-2)))
  e.extraconstr = 0
  
  #stack with all the information:
  mystack=inla.stack(data=list(y=y.inla), A=list(A.obs, A.dist, 1), effects=list(list(spatial=spatial), list(idx.spline.st=idx.spline.st), covar.inla))
  
  #model formula: ####
  #hyper.ar1=list(theta=list(prior="pccor0", param=c(.5, .5))) # a priori, probability of having rho > 0.5 is 0.5
  myform=y~-1+
    f(idx.spline.st, model = mymodel, extraconstr = list(A = A.extraconstr, e = e.extraconstr)) + 
    #f(idx.spline2.st, model = mymodel2, extraconstr = list(A = A.extraconstr, e = e.extraconstr)) +
    f(spatial, model = mymodel3) 
}

#################################################################################
#################################################################################
### Run INLA - same for all models ####
fit = inla(myform, 
           data=inla.stack.data(mystack), 
           family="gaussian",
           offset=covar.inla$y.cond,
           control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE), 
           control.predictor=list(A=inla.stack.A(mystack), compute=FALSE, link=1),
           control.inla=list(strategy="simplified.laplace", int.strategy="eb"),#,int.strategy="ccd"),
           control.fixed=list(prec=.1),
           num.thr = 2,
           verbose=TRUE)

#################################################################################
#################################################################################

### Save necessary output ####
hyper = fit$summary.hyper
fitted = fit$summary.fitted.values
summ = summary(fit)
spatial = fit$summary.random$spatial
spline = fit$summary.random$spline1d
random = fit$summary.random
obs = mystack$data$data$y
dic = fit$dic
waic = fit$waic
cpo = fit$cpo

save(fit, file = paste0("AllSpaceTimeModels-", model.num,"-fit.RData"))
save(model.num, hyper, random, fitted, summ, obs, mesh, mesh.dist, mesh.dist2, dic, waic, cpo, file=paste0("AllSpaceTimeModels-", model.num,"-OUTPUT.RData"))

#################################################################################
#################################################################################

