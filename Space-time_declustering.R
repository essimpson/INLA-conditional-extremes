##########################################################################################
### To run, this code requires:
  # data_L: nxd matrix of data (transformed to Laplace margins - see separate R file) 
  # loc: dx2 matrix of locations

##########################################################################################
### Choose threshold and conditioning site:
  u <- quantile(data_L, 0.95)

  cond <- 334
  nsites <- ncol(data_L)
  data_L <- data_L[,c(cond, c(1:nsites)[-cond])]
  loc <- loc[c(cond, c(1:nsites)[-cond]),]

##########################################################################################
### Apply runs declustering
  # Choose runs parameter by checking for stability:
    len <- NULL
    iter <- 1
    for(declust.param in c(1:30)){
      exceed <- which(data_L[,1]>u)
      clusters <- split(exceed, cumsum(c(1, diff(exceed) > declust.param)))
      len[iter] <- length(clusters); iter <- iter+1
    }
    plot(len)
    declust.param <- 12 # This is the value chosen in the paper

##########################################################################################
### Decluster observations and choose random observation in each cluster
  set.seed(1234)
  exceed <- which(data_L[,1]>u)
  clusters <- split(exceed, cumsum(c(1, diff(exceed) > declust.param)))
  obs.random <- sapply(clusters, function(x){x[1]})

  # Check that episodes don't cross years:
    check <- NULL
    for(i in 1:length(obs.random)){
      check <- c(check, year[obs.random[i]] == year[obs.random[i]+6])
    }
    table(check)
    obs.random <- obs.random[-which(check==FALSE)]

  # Indices of episodes #
    obs.ind <- NULL
    for(i in 1:length(obs.random)){
      obs.ind <- c(obs.ind, c((obs.random[i]):(obs.random[i]+6)))
    }
    
    obs <- data_L[obs.ind,]

### Create matrix of data ###
  vec1 <- c(t(obs)) # Observations (on Laplace scale)
    v2 <- matrix(vec1, nrow=length(obs.random), byrow=T)
  vec2 <- c(apply(v2, 1, function(x){x-x[1]}))
  vec3 <- rep(c(1:ncol(data)), length(obs.ind)) # Location
  vec4 <- rep(rep(c(1:7),each=nsites), length(obs.random)) # Day of week
  vec5 <- rep(c(1:length(obs.random)), each=(nsites*7)) # Episode number
  vec6 <- rep(0, length(vec5))
  vec6[1+7*nsites*c(0:(length(obs.random)-1))] <- 1 # Indicator variable for whether this is the conditioning site
  
  data <- cbind(vec1, vec2, vec3, vec4, vec5, vec6)
  colnames(data) <- c("Data","minusX0","Location","Day","Episode","Conditioning?")
  
