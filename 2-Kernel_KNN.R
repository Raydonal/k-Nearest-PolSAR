#<>--------------------------------------------------<>
#         Kernel K-Nearest Neighbors Algorithm
# 
# date: 24.11.2017
# version:0.1
# author: Jodavid Ferreira; Anny Kerollayny; Raydonal Ospina
# e-mails: jdaf1@de.ufpe.com; akgr1@de.ufpe.br; raydonal@de.ufpe.br
#<>--------------------------------------------------<>

# Prerequisites
library(parallel)

#Cores
cores <- detectCores()-1
cores

#######################################################
#######################################################
#######################################################

kernel_kNN <- function(train, test, cl, viz = 1, kernel="radial", sigma=0, p=0) {
  
  #train: training sample (data.frame);
  #test: test sample (data.frame);
  #cl: TRUE class vector;
  #viz: number Neighbors;
  
  #Condition that the number of training lines is equal to the size of the true response vector
  if( dim( train )[1] != length( cl ) ) stop( "Training dimension is different from the vector of known outcomes!")
  #<> -------------------------------------- <>
  
  #<> -------------------------------------- <>
  #### Euclidian distance
  euclideanDist <- function(x1, x2) {
    d <- (x1-x2)^2;
    return(d);
  }
  #<> -------------------------------------- <>
  
  #<> -------------------------------------- <>
  
  #<> -------------------------------------- <>
  #### Distance between test observation and training sample to find nearest neighbors  
  distNeighbors <- function(x1,vetx) {
    saida2 <- 0
    #Condition below for the training sample vector between with data.frame type
    if( is.null( dim(x1) ) & length(x1)!=1) stop("Incorrect dimension of data. Input  the training sample as data.frame type!")
    for( i in  1:dim(vetx)[2] ){
      saida <- euclideanDist( as.numeric(x1[i]) , vetx[,i])
      saida2 <- saida2+saida;
    }
    return( sqrt(saida2) )
  }
  #<> -------------------------------------- <>
  
  ##################################################################################
  
  #<> ---------------------------------------- <>
  ##              Radial kernel                ##
  #<> ---------------------------------------- <> 
  distNeighbors_radial<-function(x1,vetx,sigma) {
    saida2 <- 2 - 2*exp( - (((distNeighbors(x1 , vetx))^2) / (2*(sigma^2)) ))
    return(saida2)
  }
  #<> -------------------------------------- <>
  
  #<> ---------------------------------------- <>
  ##            Kernel Polinomial               #
  #<> ---------------------------------------- <>
  distNeighbors_polinomial<-function(x1,vetx,p) {
    saida2 <- 0
    res12 <- 0;res22 <- 0;res32 <- 0;
    if( is.null( dim(x1) ) & length(x1)!=1) stop("Incorrect dimension of data. Input  the training sample as data.frame type!")
    
    res32 <- diag(as.matrix(vetx)%*%t(vetx))
    res22 <- as.matrix(x1)%*%as.matrix(t(vetx))
    
    res12 <- as.matrix(x1)%*%as.matrix(t(x1));
    
    saida2 <-  ((1 + res12)^p)  - 2*c((1 + res22)^p) + c((1 + res32)^p) 
    return(saida2)
  }
  #<> -------------------------------------- <>
  
  ##################################################################################  
  
  #<> -------------------------------------- <>
  
  ##Loop to set the groups of the test observations
  grupofin <- c() #Vector with return the group for each observation
  saida_Neighbors <- c()
  grupofin<-simplify2array(mclapply(1:dim(test)[1],function(j) { #Loop for sample test
    #print(paste("Iteração da observação:",j))
    
    t=test[j,] #Observation o be classified
    
    if(kernel=="radial") {
      
      if(sigma==0) { stop("O valor de sigma deve ser diferente de 0") }
      
      saida_Neighbors <- data.frame(distNeighbors_radial(t,train,sigma),cl) #Calculating distance
      saida_Neighbors <- saida_Neighbors[order(saida_Neighbors[,1]),] #Sorted by distance column
      
    }else{
      
      if( kernel=="polynomial"){
        
        if(p==0){ stop("The value of 'p' must be different from 0") }
        
        saida_Neighbors <- data.frame(distNeighbors_polinomial(t,train,p),cl) #Calculating distance
        saida_Neighbors <- saida_Neighbors[order(saida_Neighbors[,1]),] #Sorted by distance column
        
        
      }else{
        
        print( "Kernel not found!" )
        
      }
    }
    
    grupos <- saida_Neighbors[1+(1:viz),2]  #Selecting nearest neighbors
    grup <- table(grupos)                   #Frequency table
    gi <- grup
    
    #grupofin[j]<-
    return( names(which.max(gi)) )## The observation is directed towards the marjoritário group.
  },mc.cores=cores))
  
  return( factor(grupofin) ) #Returns the group of each observation of the test sample
}#finished function
##################################################################################################
##################################################################################################

#Byte Code Compilation 
kernel_kNN <- compiler::cmpfun( kernel_kNN )
