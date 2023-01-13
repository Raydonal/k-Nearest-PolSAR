#<>--------------------------------------------------<>
#     fuzzy Kernel K-Nearest Neighbors Algorithm
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

kernel_f_kNN <- function(train, test, cl, viz = 1, kernel="radial", sigma=0.1) {

  #m = 2
  m=2
  
  #train: training sample (data.frame);
  #test: test sample (data.frame);
  #cl: TRUE class vector;
  #viz: number Neighbors;
  
  #Condition that the number of training lines is equal to the size of the true response vector
  if(dim(train)[1]!=length(cl)) stop("Training dimension is different from the vector of known outcomes!")
  #<> ----------------------- <>
  #<> ----------------------- <>

  #<> -------------------------------------- <>
  #### Euclidian distance
  euclideanDist <- function(x1, x2){
    d<-(x1-x2)^2
    return(d)
  }
  #<> -------------------------------------- <>
  
    
    #<> -------------------------------------- <>
    
  #<> -------------------------------------- <>
  #### Distance between test observation and training sample to find nearest neighbors  
  distNeighbors<-function(x1,vetx){
    saida2<-0
    #Condition below for the training sample vector between with data.frame type
    if( is.null( dim(x1) ) & length(x1)!=1) stop("Incorrect dimension of data. Input  the training sample as data.frame type!")
    for( i in  1:dim(vetx)[2] ){
      saida <- euclideanDist( as.numeric(x1[i]) , vetx[,i])
      saida2=saida2+saida;
    }
    return(sqrt(saida2))
  }
  #<> -------------------------------------- <>
  
    
  #<> ---------------------------------------- <>
  ##              Radial kernel                ##
  #<> ---------------------------------------- <> 
  distNeighbors_radial<-function(x1,vetx,sigma) {
    saida2 <- 2 - 2*exp( - (((distNeighbors(x1 , vetx))^2) / (2*(sigma^2)) ))
    return(saida2)
  }
  #<> -------------------------------------- <>
  
  #### Distance with formula to calculate pertinence
  distPert <- function(x1,x2) {
    if(m<=1) stop( "m must be greater than 1!" )
    if(x1==x2){ saida=0; saida } else{
      saida <- abs(x1 - x2)
      return(saida)
    }
  }
  #<> ----------------------- <>
  
    
  #<> ----------------------- <>
  ##Calculating the degrees of pertinence
  bd_p_pert<-cbind(train,cl)
  #pertinen_fin<-c()
  pertinen_fin<-  simplify2array(mclapply(1:dim(train)[1], function(j){
        t=train[j,] #Observate to be calculated pertinence
        t1<-cl[j] #Selecting Class of the j-observation 
        
        saida_Neighbors<-data.frame(distNeighbors_radial(t,train,sigma),cl) #Calculating the distance between the observation and the other observations
        saida_Neighbors<-saida_Neighbors[order(saida_Neighbors[,1]),] #Sorting through the distance column to pick up the nearest neighbors
        
    
        grupos<-saida_Neighbors[1+(1:viz),2] #Selecting nearest neighbors
        grup<-table(grupos) #viewing in a frequency table
        
        # Calculation of pertinence for training sample
        
        #Generating the first pertinences
        u<-0
        for(i in 1:length(levels(cl))){
          if(names(grup[i])==t1){
            u[i]<-0.51+(grup[i]/viz)*.49
          }else{
            u[i]<-(grup[i]/viz)*.49
          }
        }
          
            #pertinen_fin<-
            return(u)
      },mc.cores = cores))
      pertinen_fin <- t(pertinen_fin)
      bd_p_pert<-cbind(bd_p_pert,pertinen_fin)
      #<> ----------------------- <>
      
    
    
    
    
      ##Loop to set the groups of the test observations
      #grupofin<-c()
      grupofin<-simplify2array(mclapply(1:dim(test)[1],function(j) { #Loop for sample test
        
      
      t=test[j,] #Observation to be classified
      
      saida_Neighbors<-data.frame(sqrt(distNeighbors_radial(t,train,sigma)),cl,train,pertinen_fin) #Calculating the distance between the observation and the other observations
      saida_Neighbors<-saida_Neighbors[order(saida_Neighbors[,1]),]  #Sorting through the distance column to pick up the nearest neighbors
      
      grupos<-saida_Neighbors[1+(1:viz),-1]  #Selecting the nearest neighbors

      # Calculation of pertinence for test sample
      
      ##Defining the gi's
      U<-c()
      
      for(i in 1:length(levels(cl))){
      coluna_i<-1+dim(train)[2]+i  #The first one is of column cl, this is necessary because the column of pertinence is after the column cl + variables_variables
      
      
      
      sai3=0;sai4=0;
      for(l in 1:viz){      
      for(k in 1+(1:dim(train)[2])){
        sai6<-(distPert(t[(k-1)],grupos[l,k]));
        sai5 <- ( 1 / ( 1 - exp( - ( sai6 / (2* (sigma^2)) ) ) ))
        if(sai5==Inf){sai5=.Machine$integer.max}
        sai<-grupos[l,coluna_i] * sai5
        sai2 <- ( 1 / ( 1 - exp( - ( (sai6 / (2* (sigma^2)) ) ) )))
        if(sai2==Inf){sai2=.Machine$integer.max}
        sai3<-sai3+sai;
        sai4=sai4+sai2;
        }
      }  
      
    
     
      U[i] = ifelse(sai3==0 & sai4==0, 1, sai3/sai4);
      }
      
      #grupofin<-
      return(U)
      },mc.cores=cores))
      
      return(t(grupofin))
}#finished function kernel_f_kNN
##################################################################################################
##################################################################################################


#Byte Code Compilation 
kernel_f_kNN <- compiler::cmpfun( kernel_f_kNN )
