#################################################################################
### Function to transform data to Laplace margins:
  library(ismev)

  laplaceMargins <- function(univariateData, thrQ=0.95){
      
    size <- length(univariateData)
    empiricalMargins <- rank(univariateData)/(size+1)
      
    empiricalMargins[univariateData>quantile(univariateData, thrQ)] <- NA
      
    fit <- gpd.fit(univariateData, threshold=quantile(univariateData,thrQ))$mle
    scale <- fit[1]
    shape <- fit[2]
      
    empiricalMargins[is.na(empiricalMargins)] <- 1 - pgpd(univariateData[univariateData>quantile(univariateData, thrQ)], u=quantile(univariateData, thrQ), scale, shape, lower.tail=FALSE)*(1-thrQ)
      
    laplaceMargins <- rep(NA, length(empiricalMargins))
    laplaceMargins[empiricalMargins<=0.5] <- log(2*empiricalMargins[empiricalMargins<=0.5])   
    laplaceMargins[empiricalMargins>0.5]  <- -log(2*(1-empiricalMargins[empiricalMargins>0.5]))   
      
    return(laplaceMargins)
  }

#################################################################################
### Apply transformation separately at each spatial location:
  data_L <- apply(data, 2, laplaceMargins)
  
