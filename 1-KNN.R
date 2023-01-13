#<>--------------------------------------------------<>
#  K-Nearest Neighbors Algorithm
# 
# date: 2023
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

kNNII<-function(train, test, cl, viz = 1){
  
  #train: training sample (data.frame);
  #test: test sample (data.frame);
  #cl: TRUE class vector;
  #viz: number Neighbors;
  
  #Condition that the number of training lines is equal to the size of the true response vector
  if(dim(train)[1]!=length(cl)) stop("Training dimension is different from the vector of known outcomes!")
    #<> ----------------------- <>
    #<> ----------------------- <>

#    #<>----------------------------------------------------<>
#    #Normalizing as training and test variables
#    for(i in 1:dim(train)[2]){
#      train[,i] <- (train[,i]-min(train[,i]) ) / (max(train[,i])-min(train[,i]) )
#      test[,i]  <- (test[,i]-min(test[,i]) ) / (max(test[,i])-min(test[,i]) )
#    }
#    
#    #<> -------------------------------------- <>
      #### Euclidian distance
      euclideanDist <- function(x1, x2){
            d<-(x1-x2)^2
        return(d)
      }
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
    
    ##Loop to set the groups of the test observations
    grupofin<-c() #Vector with return the group for each observation
    
    grupofin<-simplify2array(mclapply(1:dim(test)[1],function(j) { #Loop for sample test
      
      t=test[j,] #observation to be classified
      
      saida_Neighbors<-data.frame(distNeighbors(t,train),cl) #Calculating distance
      saida_Neighbors<-saida_Neighbors[order(saida_Neighbors[,1]),] #Sorted by distance column
    
        grupos<-saida_Neighbors[1+(1:viz),2] #Selecting nearest neighbors
        grup<-table(grupos) #Frequency table

		gi<-grup
		
		#grupofin[j]<-
		return( names(which.max(gi)) )## The observation is directed towards the marjoritário group.
    },mc.cores=cores))
    
    return(factor(grupofin)) #Returns the group of each observation of the test sample
 
}#finished function
##################################################################################################
##################################################################################################

#Byte Code Compilation 
kNNII <- compiler::cmpfun(kNNII)
