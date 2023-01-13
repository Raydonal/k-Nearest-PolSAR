#<>--------------------------------------------------<>
#  Kullback-Leibler Distance
# 
# date: 30.07.2020
# version:0.1
# author: Jodavid Ferreira; Anny Kerollayny; Raydonal Ospina
# e-mails: jdaf1@de.ufpe.com; akgr1@de.ufpe.br; raydonal@de.ufpe.br
#<>--------------------------------------------------<>
library(MASS)

#<> -------------------------------------------------------------------- <>
# ------------------------------ #
#   Gerando Parâmetros  - SIGMAS #
# ------------------------------ #
##Arquivo Unico gerar parametros

gerarpar<-function(train,cl){
  
  ############### Inicio trabalho com amostra de treinamento
  
  inf<-3#dim(train)[2] #Será sempre 3?
  
  #Quantidades de grupos
  qg<-length(levels(cl))
  
  
  vet<-data.frame(x=qg,y=2) 
  for(i in 1:qg){
    vet[i,]<-c(levels(cl)[i],sum(ifelse(cl==levels(cl)[i],1,0)))
  }
  
  ##Criando matrizes
  
  matrizes<-matrix(rep(NA,(qg*inf)*inf),(qg*inf),inf)
  
  j=1; #variável auxiliar
  for(i in 1:qg){
    matrizmeanf<-0
    mt<-train[which(cl==i),]
    maxl<-as.numeric(vet[i,2])
    
    for(t in 1:maxl){
      matrizmean<-sep.HVPolSar(mt,linhadata=t,type=2)
      matrizmeanf<-matrizmeanf+matrizmean
      if(t==maxl){
        matrizmeanf<-matrizmeanf/maxl
      }
    }
    matrizes[j:(i*inf),]<-matrizmeanf
    j=j+inf; #necessário para colocar matrizes na posição correta
  }
  
  return(list(matriz=matrizes,ngroups=qg))
  
}

#<> -------------------------------------------------------------------- <>

#<> -------------------------------------------------------------------- <>
# -------------------------------------------- #
#       Função para separar como matrizes      #
# -------------------------------------------- #
#Tamanho da matriz vai ser o numero inteiro mais proximo da sua raiz quadrada + 1.
sep.HVPolSar<-function(data,linhadata=1,colunadata=1,type=1){
  
  complexconj<-function(x) Re(x)-Im(x)*(1i) #Calcular o número complexo conjugado para matriz
  
  if(type==1){  
    #Número de linhas e colunas da matriz
    n=(round(sqrt(length(data[linhadata,colunadata,])))+1);m=n
    
    
    p1<-data[linhadata,colunadata,1];p4<-data[linhadata,colunadata,4];p7<-data[linhadata,colunadata,5]
    p2<-complexconj(data[linhadata,colunadata,4]);p5<-data[linhadata,colunadata,2];p8<-data[linhadata,colunadata,6]
    p3<-complexconj(data[linhadata,colunadata,5]);p6<-complexconj(data[linhadata,colunadata,6]);p9<-data[linhadata,colunadata,3]
    vet<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9)
    matrizpixel<-matrix(vet,nrow=n,ncol=m)
  }else{
    
    n=(round(sqrt(length(data[linhadata,])))+1);m=n
    p1<-data[linhadata,1];p4<-data[linhadata,4];p7<-data[linhadata,5]
    p2<-complexconj(data[linhadata,4]);p5<-data[linhadata,2];p8<-data[linhadata,6]
    p3<-complexconj(data[linhadata,5]);p6<-complexconj(data[linhadata,6]);p9<-data[linhadata,3]
    vet<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9)
    matrizpixel<-matrix(vet,nrow=n,ncol=m)
  }
  return(matrizpixel)
}
#<> -------------------------------------------------------------------- <>

#<>-------------------------------------------------------<>
# Encontrando vizinhança do pixel
neighbours <- function(imagem,n=1, pix.l=1,pix.c=1){
  
  l <- dim(imagem)[6];
  test_janela<-c()
  for(l in 1:6){
    matriz.dados<-imagem[,,l]
    
    #selecionando matriz
    #Condicao necessaria para valores negativos
    lin1<-ifelse(pix.l-n<0,0,pix.l-n)
    lin2<-ifelse(pix.l+n>dim(matriz.dados)[1],dim(matriz.dados)[1],pix.l+n)
    col1<-ifelse(pix.c-n<0,0,pix.c-n)
    col2<-ifelse(pix.c+n>dim(matriz.dados)[2],dim(matriz.dados)[2],pix.c+n)
    
    matriz<- matriz.dados[lin1:lin2,col1:col2] 
    
    test_janela<-cbind(test_janela,c(matriz))#n=1 -> 9 pixels
    
  }
  
  return(test_janela)
  
}
#<>-------------------------------------------------------<>

# ------------------------------------ #
# <> Determinant for matrix Complex <> #
# ------------------------------------ #

# Determinant for matrix with order 3
detcomplex<-function(A){
  # parameters: A: Hermitian positive definite matrix 3 x 3
  # return: determinant
  R=Re(A); I=Im(A);
  g1=R[1,1]; g2=R[2,2]; g3=R[3,3]; g4=R[1,2]; g5=I[1,2]; g6=R[1,3];
  g7=I[1,3]; g8=R[2,3]; g9=I[2,3];
  det1=(g1*g2*g3)-(g3*(g4^2))-(g3*(g5^2))-(g2*(g6^2))-(g2*(g7^2))-(g1*(g8^2))-(g1*(g9^2))+(2*g4*g6*g8)+(2*g5*g7*g8)-(2*g5*g6*g9)+(2*g4*g7*g9);
  (det1)}



# ---------------------------------------------------------------- #
#   Função para cálcular média da vizinhança, incluindo o pixel    #
# ---------------------------------------------------------------- #

mean.pixel<-function(matriz.dados,n=1, pix.l=1,pix.c=1){
  
  #Pixel selecionado
  #matriz.dados[pix.i,pix.j]
  
  #selecionando matriz
  #Condicao necessaria para valores negativos
  lin1<-ifelse(pix.l-n<0,0,pix.l-n)
  lin2<-ifelse(pix.l+n>dim(matriz.dados)[1],dim(matriz.dados)[1],pix.l+n)
  col1<-ifelse(pix.c-n<0,0,pix.c-n)
  col2<-ifelse(pix.c+n>dim(matriz.dados)[2],dim(matriz.dados)[2],pix.c+n)
  
  matriz<- matriz.dados[lin1:lin2,col1:col2] 
#  matriz[which(matriz==matriz.dados[pix.l,pix.c])[1]] <- NA #Removendo elemento da matriz
  
  #print(matriz) #Verificando se a matriz selecionada é a correta
  return(mean=mean(matriz,na.rm = T)) #retorna a média 
}



# Pacotes Necessários
library(doSNOW)
library(doParallel)
#<>--------------------------------------------------------------<>
#Determinando os Cores
cl <- 2 #Quantidade de Cores
registerDoParallel(cores=cl)    
#<>---------------------------------------------------------------<>    

# ------------------------------------------------------------------- #
#   Função para classificação distância de Kullback-Leibler Complexa   #
# ------------------------------------------------------------------- #

kl.distance<-function(test,cltest,matrizes_par,qg,nlooks=0,nlinhas,ncolunas,k=1,type=1){
  #  print(c(pixl,pixc))
  
  inf<-ifelse(dim(test)[2]>3,3,dim(test)[2]) #Obtendo tamanho da dimensão, máximo é 3
  
	Classificacaofinal <- list();

	if(type==2){
	varaux <- unique(levels(cltest));
	}else{
	varaux = 1;	
	}

 	for( gr in varaux){

  # --------------------------------- #
  #   Converterndo Colunas em Matrizes #
  # -------------------------------- #
  n=k#tamanho da janela k=1=3x3
  imagem <- list()

 if(type==2){
  for(l in 1:6){
    imagem[[l]]<-matrix(test[which(cltest==gr),l],nrow=nlinhas,ncol=ncolunas)
    
  }
  }else{
  for(l in 1:6){
    imagem[[l]]<-test[,,l]#matrix(test[,l],nrow=nlinhas,ncol=ncolunas)
    
  }
	}


  #Matrizes para indexação do i e j
  matrizlinha <- matrix(rep(1:nlinhas,ncolunas),nrow = nlinhas)
  matrizcoluna <- matrix(rep(1:ncolunas,nlinhas),ncol = ncolunas,byrow = T)
 
	



 	Classificacao <- foreach(h=1:(nlinhas*ncolunas),.combine=c) %dopar% {
    
		pixl <- matrizlinha[h];
		pixc <- matrizcoluna[h];    
        


	test_janela<-c()
	for(l in 1:6){ 
	test_janela<-cbind(test_janela,mean.pixel(imagem[[l]],n=n,pix.l=pixl,pix.c=pixc))#n=1 -> 9 pixels
	}

  sigma1<-sep.HVPolSar(data = test_janela,linhadata=1,colunadata=1,type=2)#linhadata e colunadata sempre 1
  

  #Inicio do For dessa função
  Drd1<-matrix(rep(NA,(2*qg)),qg)
  
  j=1; #variável auxiliar
  h=inf #segunda variável auxiliar
  L=nlooks#Valor de n
  for(i in 1:qg){
   
    
    sigma2= matrizes_par[j:h,];
    
    
    Drd1[i,] = c(i,
					L * ( 
						( ( sum( diag( ginv(sigma1)%*%sigma2 + ginv(sigma2)%*%sigma1 ) ) ) / 2 )
						- inf)				
    )#fechamento do vetor
    j=j+inf; #inicio da posicao da linha da matriz de covariancia
    h=h+inf; # fim da matriz de covariância
    
    
  }#fim do for
  Drd1 <- Re(Drd1);
  #Classificação pelas distâncias de Kullback-Leibler
  Classific <- Drd1[which.min(Drd1[,2])]
  
	return(Classific)

	}

	Classificacaofinal[[gr]] <- Classificacao
	
	} # Fim do for dos grupos	
  
	Classificacao = c(unlist(Classificacaofinal));
	names(Classificacao) <-c()

  return(Classificacao = Classificacao)
  
}

