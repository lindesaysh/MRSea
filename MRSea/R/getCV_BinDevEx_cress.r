
getCV_BinDevEx_cress<-function(data, baseModel, splineParams){

#baseModel <- salsa2dOutput_k2$bestModel        # test only
#splineParams <- salsa2dOutput_k2$splineParams  # test only
#data <- fh_dat# test only
 
  d2k<-splineParams[[1]]$dist
  radiusIndices <-splineParams[[1]]$radiusIndices
  radii <- splineParams[[1]]$radii
  aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
  
  attributes(baseModel$formula)$.Environment<-environment()
  # calculate cross-validation
  nfolds<-length(unique(data$foldid))
  dists<-d2k
  store<- matrix(0, nrow=nfolds, ncol=1)
  
  for(f in 1:nfolds){
    tot  <- length(data$response[data$foldid==f])   # djfr make vector of a certain size
    Ntot <- length(data$response[data$foldid==f])    # djfr make vector of a certain size
    predstrue  <- length(data$response[data$foldid==f])  # djfr make vector of a certain size
    dists<-d2k[data$foldid!=f,]
    foldedFit<- update(baseModel, .~., data=data[data$foldid!=f,])
    if(length(coef(foldedFit))==1){
      dists<-d2k[data$foldid==f,]
      if(is.null(baseModel$offset)){
        predscv<- baseModel$family$linkinv(as.matrix(model.matrix(baseModel)[data$foldid==f])%*%coef(foldedFit))
      }else{
        predscv<- baseModel$family$linkinv(as.matrix(model.matrix(baseModel)[data$foldid==f])%*%coef(foldedFit)) * baseModel$family$linkinv(baseModel$offset)[data$foldid==f]  
      }
    }else{  
      dists<-d2k[data$foldid==f,]
      if(is.null(baseModel$offset)){
        predscv<- baseModel$family$linkinv(model.matrix(baseModel)[data$foldid==f,]%*%coef(foldedFit))
      }else{
        predscv<- baseModel$family$linkinv(model.matrix(baseModel)[data$foldid==f,]%*%coef(foldedFit)) * baseModel$family$linkinv(baseModel$offset)[data$foldid==f]
      }
    } 
    
        
   predstrue <- data$response[data$foldid==f]        # djfr
   predscv <- ifelse(predscv<1,predscv,0.9999)       # djfr
   predscv <- ifelse(predscv>0,predscv,0.0001)       # djfr
   Nu <- rep(sum(predstrue)/length(predstrue),length(predstrue))  # djfr what is proportion of 1s in true data? 

for(z in 1:length(predstrue)){                                         # djfr
tot[z] <- predstrue[z]*log(predscv[z])+(1-predstrue[z])*log(1-predscv[z])  # djfr
}                                                                         # djfr
Rdev <- -2*sum(tot)                                                       # djfr

for(x in 1:length(predstrue)){                                             # djfr
Ntot[x] <- predstrue[x]*log(Nu[x])+(1-predstrue[x])*log(1-Nu[x])           # djfr
}                                                                          # djfr
Ndev <- -2*sum(Ntot)                                                        # djfr
   store[f]<- (Ndev-Rdev)/Ndev * 100                                       # djfr
  }                                                                        # djfr
  return(store)                                                    # djfr
}


