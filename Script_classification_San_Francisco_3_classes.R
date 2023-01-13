#<>--------------------------------------------------<>
#  Classification script
# 
# date: 30.01.2019
# update: 20.07.2020
# version:0.1
# author: Jodavid Ferreira; Anny Kerollayny; Raydonal Ospina
# e-mails: jdaf1@de.ufpe.com; akgr1@de.ufpe.br; raydonal@de.ufpe.br
#
#
# Output: 'lresults' 
# it is a list with the results of the accuracy for the studied parameters of each method,
# each output is a method, following an order of execution.
# order to list:
# 1 - kNN
# 2 - Kernel kNN
# 3 - Fuzzy kNN
# 4 - Kernel Fuzzy kNN
# 5 - Naive Bayes
# 6 - SVM
# 7 - Xboost
# 8 - KL distance
#<>--------------------------------------------------<>
# Packages
library(lattice)
library(factoextra)
library(rgl)
library(viridis)
library(FactoMineR)
library(tidyverse)
library(reshape2)
library(e1071)
library(caret)
#--------------------------------------------------------------------------- #
#--------------------------------------------------------------------------- #
# Kappa function to measure accuracy
#--------------------------------------------------------------------------- #
kappaclas<-function(matriz_classificacao){

  mc = matriz_classificacao;
  mu=mc;
  total =0;
  total_linha=matrix(rep(0,(dim(mu)[1]^2)),ncol=dim(mu)[1])
  total_coluna=matrix(rep(0,(dim(mu)[1]^2)),ncol=dim(mu)[1])
  
  for( i in 1:dim(mu)[1]){
    for(j in 1:dim(mu)[1]){
      total_linha[i]=total_linha[i]+mc[i,j];
      total_coluna[i]=total_coluna[i]+mc[j,i];
    }
    total=total+total_linha[i];
  }
  
  teta1 = 0;
  teta2 = 0;
  teta3 = 0;
  teta4 = 0;
  
  for(i in 1:dim(mu)[1]){
    teta1=teta1+mc[i,i];
    teta2=teta2+total_linha[i]*total_coluna[i];
    teta3=teta3+mc[i,i]*(total_linha[i]+total_coluna[i]);
    for(j in 1:dim(mu)[1]){
      teta4=teta4+mc[i,j]*(total_linha[i]+total_coluna[j])^2;
    }
  }
  
  teta1=teta1/total;
  teta2=teta2/(total^2);
  teta3=teta3/(total^2);
  teta4=teta4/(total^3);
  
  term1 = (teta1 * (1 - teta1)) / ((1 - teta2)^2);
  term2 = ((2 * (1 - teta1)) * (2 * teta1 * teta2 - teta3)) / ((1 - teta2)^3);
  term3 = (((1 - teta1)^2) * (teta4 - 4 * teta2 * teta2)) / ((1 - teta2)^4);
  
  kappa=(teta1 - teta2)/(1-teta2)*100
  VarK = (term1 + term2 + term3) / total

  return(list(kappa=kappa,VarianciaKappa=VarK))
}

####################################################################################
#<>-------------------------------------------------<>
##################################################################################################
#               Reading the function KNN - k-Nearest Neighbour Classification                    #
source("1-KNN.R")
##################################################################################################
##################################################################################################
#               Reading the function Kernel KNN - k-Nearest Neighbour Classification             #
source("2-Kernel_KNN.R")
##################################################################################################
##################################################################################################
#         Reading the function KNN - Fuzzy k-Nearest Neighbour Classification                    #
source("3-FkNN.R")
##################################################################################################
##################################################################################################
#         Reading the function KNN - Fuzzy k-Nearest Neighbour Classification                    #
source("4-Kernel_FkNN.R")
##################################################################################################
##################################################################################################
#         Reading the function Kullback-Leibler Distance Classification                    #
source("5-KL_distance.R")
##################################################################################################
##################################################################################################
#               Reading the function imagematrix							                     #
source("imagematrix_code.R")
##################################################################################################

# Results List
lresults <- list()

#<>-----------------<>
## San Francisco
#<>-----------------<>

# load image
load("../Application-Images/AirSAR_SanFrancisc_Enxu.RData")
image <- San_Francisc_Enxuto; # pattern the code
# - 245 linhas
#<>-------------------------------------------------------------<>
IHH <- ImHH <- Re(image[245:900,,1])
IHV <- Re(image[245:900,,2])
IVV <- Re(image[245:900,,3])
#<>---<>
ImHHeq <- matrix(ecdf(IHH)(IHH), nrow=nrow(IHH), ncol=ncol(IHH))
ImHVeq <- matrix(ecdf(IHV)(IHV), nrow=nrow(IHH), ncol=ncol(IHH))
ImVVeq <- matrix(ecdf(IVV)(IVV), nrow=nrow(IHH), ncol=ncol(IHH))
#<>---<>
# Pauli
#<>---<>
Blue <- matrix(abs(IHH+IVV), nrow=nrow(IHH), ncol=ncol(IHH))
Blue <- Blue/max(Blue)
Red <- matrix(abs(IHH-IVV), nrow=nrow(IHH), ncol=ncol(IHH))
Red <- Red/max(Red)
Green <- matrix(IHV, nrow=nrow(IHV), ncol=ncol(IHV))
Green <- Green/max(Green)
#<>---<>
IntensidadesEqualizadas_r <- array(0, c(nrow=nrow(Blue), ncol=ncol(Blue), 3))
IntensidadesEqualizadas_r[,,1] <- ImHHeq
IntensidadesEqualizadas_r[,,2] <- ImHVeq
IntensidadesEqualizadas_r[,,3] <- ImVVeq
#<>---<>
#Visualizar pedaço da imagem colorida
IntensidadesEqualizadas_r <- imagematrix(IntensidadesEqualizadas_r)
#<>---<>
png("../results/San Francisco/San_Francisco_amostras_new.png",width=600)
plot(IntensidadesEqualizadas_r^.7)
#Selecting the training datasets
#Area 1 - Water
 x.minotr = 20; 
 x.maxotr = 50; 
 y.minotr = 400 - 245; 
 y.maxotr = 430 - 245;
 rect(xleft = x.minotr,ybottom = nrow(ImHH)-y.minotr,xright = x.maxotr, ytop = nrow(ImHH)-y.maxotr,col="NA", border = "red", lwd = 2)
 #meios
 mx <- (x.minotr+x.maxotr)/2
 my <- (2*nrow(ImHH)-(y.minotr+y.maxotr))/2
 am <- 30
 lines(mx:(mx+am),my:(my+am),col="red", lwd=1) 
 text(x = (mx+am+10), (my+am+10) ,"1",cex=1,col="red")
#
#Area 2 - Vegetation
 x.minftr = 450; 
 x.maxftr = 480; 
 y.minftr = 350 - 245;
 y.maxftr = 380 - 245;
 rect(xleft = x.minftr,ybottom = nrow(ImHH)-y.minftr,xright = x.maxftr, ytop = nrow(ImHH)-y.maxftr,col="NA", border = "red", lwd = 2) 
 #meios
 mx <- (x.minftr+x.maxftr)/2
 my <- (2*nrow(ImHH)-(y.minftr+y.maxftr))/2
 am <- 30
 lines(mx:(mx+am),my:(my+am),col="red", lwd=1) 
 text(x = (mx+am+10), (my+am+10) ,"2",cex=1,col="red")
# 
#Area 4 - Low-Density Urban
 x.minldtr = 400;
 x.maxldtr = 430; 
 y.minldtr = 500 - 245; 
 y.maxldtr = 530 - 245;
rect(xleft = x.minldtr,ybottom = nrow(ImHH)-y.minldtr,xright = x.maxldtr, ytop = nrow(ImHH)-y.maxldtr,col="NA", border = "red", lwd = 2)
 #meios
 mx <- (x.minldtr+x.maxldtr)/2
 my <- (2*nrow(ImHH)-(y.minldtr+y.maxldtr))/2
 am <- 30
 lines(mx:(mx+am),my:(my+am),col="red", lwd=1) 
 text(x = (mx+am+10), (my+am+10) ,"3",cex=1,col="red")
#
# #<>------------------------------------------------<>
#____________________________________________________#
# Selecting the test datasets

#Area 1 - Water
 x.minote = 20; 
 x.maxote = 50; 
 y.minote = 700 - 245; 
 y.maxote = 730 - 245;
 rect(xleft = x.minote,ybottom = nrow(ImHH)-y.minote,xright = x.maxote, ytop = nrow(ImHH)-y.maxote,col="NA", border = "blue", lwd = 2)
 #meios
 mx <- (x.minote+x.maxote)/2
 my <- (2*nrow(ImHH)-(y.minote+y.maxote))/2
 am <- 30
 lines(mx:(mx+am),my:(my+am),col="blue", lwd=1) 
 text(x = (mx+am+10), (my+am+10) ,"1",cex=1,col="blue")
#Area 2 - Vegetation
 x.minfte = 450; 
 x.maxfte = 480; 
 y.minfte = 600 - 245;
 y.maxfte = 630 - 245;
 rect(xleft = x.minfte,ybottom = nrow(ImHH)-y.minfte,xright = x.maxfte, ytop = nrow(ImHH)-y.maxfte,col="NA", border = "blue", lwd = 2)
 #meios
 mx <- (x.minfte+x.maxfte)/2
 my <- (2*nrow(ImHH)-(y.minfte+y.maxfte))/2
 am <- 60
 lines(mx:(mx-am),my:(my-am),col="blue",lwd=1) 
 text(x = (mx-am-10), (my-am-10) ,"2",cex=1,col="blue")
#Area 4 - Low-Density Urban
 x.minldte = 800;
 x.maxldte = 830; 
 y.minldte = 330 - 245; 
 y.maxldte = 360 - 245;
rect(xleft = x.minldte,ybottom = nrow(ImHH)-y.minldte,xright = x.maxldte, ytop = nrow(ImHH)-y.maxldte,col="NA", border = "blue", lwd = 2)
 #meios
 mx <- (x.minldte+x.maxldte)/2
 my <- (2*nrow(ImHH)-(y.minldte+y.maxldte))/2
 am <- 30
 lines(mx:(mx+am),my:(my+am),col="blue", lwd=1) 
 text(x = (mx+am+10), (my+am+10) ,"3",cex=1,col="blue")
dev.off()

image <- image[245:900,,]

Waterdata <- sapply(1:3, function(x) c( Re( image[y.minotr:y.maxotr , x.minotr:x.maxotr , x] ) ) )
Vegetationdata <- sapply(1:3, function(x) c( Re( image[y.minftr:y.maxftr , x.minftr:x.maxftr , x] ) ) )
LowUrbandata <- sapply(1:3, function(x) c( Re( image[y.minldtr:y.maxldtr , x.minldtr:x.maxldtr , x] ) ) )

dados<-rbind(Waterdata,Vegetationdata,LowUrbandata) #datasets
Verd <- factor(c(rep("1",dim(Waterdata)[1]), rep("2",dim(Vegetationdata)[1]),rep("3",dim(LowUrbandata)[1]))) #Converting character variable to categorical

#____________________________________________________#
# Selecting the test datasets

Waterdata <- sapply(1:3, function(x) c( Re( image[y.minote:y.maxote , x.minote:x.maxote , x] ) ) )
Vegetationdata <- sapply(1:3, function(x) c( Re( image[y.minfte:y.maxfte , x.minfte:x.maxfte , x] ) ) )
LowUrbandata <- sapply(1:3, function(x) c( Re( image[y.minldte:y.maxldte , x.minldte:x.maxldte , x] ) ) )

dadostest<-rbind(Waterdata,Vegetationdata,LowUrbandata) #datasets
Verdtest <- factor(c(rep("1",dim(Waterdata)[1]), rep("2",dim(Vegetationdata)[1]), rep("3",dim(LowUrbandata)[1]))) #Converting character variable to categorical
#___________________________________________________#

train<-treino<-data.frame(dados)
cl<-verdadeiro<-Verd
test<-teste<-data.frame(dadostest)
verdadeiro_teste<-Verdtest



#<>---------------------------------------------------------------------------------------------------------------------<>#
# Principal Components
#<>---------------------------------------------------------------------------------------------------------------------<>#
fit <- princomp(dadostest, cor=F)

colors = ifelse(Verd==1,'#453781FF',
		 ifelse(Verd==2, '#287D8EFF','#FDE725FF')) # Color Vector
colors = factor(colors,levels=c('#453781FF', '#287D8EFF', '#FDE725FF'))

regioes = ifelse(Verd==1,"Water",
		 ifelse(Verd==2, "Forest","Urban")) # Regions Vector
regioes = factor(regioes,levels=c("Water","Forest","Urban"))

pch_edit = ifelse(Verd==1,19,
		 ifelse(Verd==2, 17,15)) # Regions Vector

#--------------------------------------------------------------------------- #
# Plot - Principal Components 2D

grDevices::cairo_pdf("../results/San Francisco/principal_component_San_Francisco_2d.pdf", 
                     width=4.1, height=4.1, pointsize = 28)
fviz_pca_ind(fit,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = regioes,#paste("Region",Verdtest), # color by groups
             fill.ind = regioes,#,
			#palette = levels(colors),  
             addEllipses = FALSE, # Concentration ellipses
             legend.title = "",
             title="",
             alpha.ind = 0.6, #transparencia
               pointsize= 2.0,
             labelsize = 9.0) +
viridis::scale_colour_viridis(discrete = TRUE, option = "D")+
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), 
    legend.text = element_text(colour="black"),legend.position = "top", legend.box = "vertical" ) +
 xlab("First component")+
  ylab("Second component") +
xlim(-0.5, 4) +
ylim(-1.5,1.5)
 dev.off()


# Plot - Principal Components 3D
library(lattice)

pdf("../results/San Francisco/principal_component_San_Francisco_3d.pdf",width=7, height=7)
# par.set <-
#   list(axis.line = list(col = "transparent"),
#        clip = list(panel = "off"))
cloud(fit$scores[,3]~fit$scores[,1]*fit$scores[,2], zlab="Third\n  component",
      xlab="First\n component", ylab = "Second\n  component",
      #screen=list(x=0,y=50,z=0),
      # par.settings = par.set,
      screen=list(x=35,y=40,z=25),
      col=as.character(colors),
      pch=pch_edit,
      distance = .35, zoom = .6,
      scales = list( arrows = TRUE, col="darkblue"),
      par.settings = list(axis.line = list(col = 'transparent'))
      )
dev.off()


#-----------------------    Violin by channels ------------------------------ 
library(rgl)
library(viridis)
library(ggplot2)
library(FactoMineR)
library(tidyverse)
library(factoextra)
library(reshape2)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


#-----------------------    Violin by channels ------------------------------ 

#até 3º quartil
HH <- subset(c(image[,,1]), c(Re(image[,,1])) <=  summary(c(Re(image[,,1])))[5]  ) #Pegando amostra ate o valor do terceiro quartil
HV <- subset(c(image[,,2]), c(Re(image[,,2])) <=  summary(c(Re(image[,,2])))[5]  ) #Pegando amostra ate o valor do terceiro quartil
VV <- subset(c(image[,,3]), c(Re(image[,,3])) <=  summary(c(Re(image[,,3])))[5]  ) #Pegando amostra ate o valor do terceiro quartil


#dadostest1 <- test %>% mutate(verdadeiro_teste)
#colnames(dadostest1) <- c("HH","HV","VV","verdadeiro_teste")
#violin1 <- dadostest1[,-4] %>% melt() %>% data.frame()

violin1 <- data.frame(variable = c(rep("HH",length(HH)),rep("HV",length(HV)),rep("VV",length(VV))),value=c(HH,HV,VV))
head(violin1) 

pl <- ggplot(violin1,aes(x = variable, y = value, col="red")) +
  geom_violin(adjust = 1.0, size=2, show.legend = FALSE) +
  #viridis::scale_colour_viridis(discrete = TRUE, option = "D")+
scale_fill_manual(values = c("red", "red", "red"))+   
geom_boxplot(col="black",aes(fill=variable), notch=FALSE, alpha=0.4, width=0.22, 
               outlier.colour = "red", outlier.shape = 1,outlier.size = 5) +
  facet_grid(. ~ variable,  scales='free', space = "free")+
  xlab("")+ggtitle("Channels")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("")+
  theme(
    panel.grid.major = element_line(colour = "gray", 
                                    linetype = "dotted"),
    panel.background = element_rect(fill = "white", 
                                    colour="black"),
    strip.text.x = element_text(size=18, 
                                hjust=0.5, 
                                vjust=0.5,
                                face="bold", lineheight = 0.5),
    strip.text.y = element_text(size=18, 
                                hjust=0.5, 
                                vjust=0.5,
                                face="bold", lineheight = 0.5),
    strip.background = element_rect(colour="black", fill="gray98"),
    axis.text=element_text(size=18, face="bold", colour="gray24"),
    
    axis.text.x  = element_blank(),
    #axis.text.x  = element_text(angle=0, vjust=0.5, size=18,colour="white"),
    axis.title=element_text(size=18,face="bold"),
    plot.title = element_text(size = 18, colour = "black", face="bold", hjust=0.5),
    legend.position ="none", 
    panel.spacing = unit(2, "lines")
  )+guides(color=FALSE)+scale_y_continuous(label=scientific_10)+ylim(c(0,max(violin1$value)))

ggsave("../results/San Francisco/plot_channels_San_Francisco.png", pl, width=10, height =8.5, scale=0.8)



#-----------------------    Violin by Regions ------------------------------ 
library(rgl)
library(viridis)
library(ggplot2)
library(FactoMineR)
library(tidyverse)
library(factoextra)
library(reshape2)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#até 3º quartil
WaterdataV <- subset(c(Waterdata[,1]), c(Waterdata[,1]) <=  summary(c(Waterdata[,1]))[5]  ) #Pegando amostra ate o valor do terceiro quartil
VegetationdataV <- subset(c(Vegetationdata[,1]), c(Vegetationdata[,1]) <=  summary(c(Vegetationdata[,1]))[5]  ) #Pegando amostra ate o valor do terceiro quartil
LowUrbandataV <- subset(c(LowUrbandata[,1]), c(LowUrbandata[,1]) <=  summary(c(LowUrbandata[,1]))[5]  ) #Pegando amostra ate o valor do terceiro quartil



dadostest1 <- data.frame(WaterdataV,VegetationdataV,LowUrbandataV) #datasets
colnames(dadostest1) <- c("Ocean","Forest","Urban")

violin1 <- dadostest1 %>% melt() %>% data.frame()
head(violin1) 

pl <- ggplot(violin1,aes(x = variable, y = value, col=variable)) +
  geom_violin(adjust = 1.0, size=2, show.legend = FALSE) +
  viridis::scale_colour_viridis(discrete = TRUE, option = "D")+
  geom_boxplot(col="darkblue",aes(fill=variable), notch=FALSE, alpha=0.4, width=0.22, 
               outlier.colour = "red", outlier.shape = 1,outlier.size = 5) +
  # geom_jitter(aes(color=variable), alpha=0.6, shape=21, width = 0.03, height =  0.06)+
  stat_summary(fun=median, geom="point", size=1)+ #, color="red")
  viridis::scale_fill_viridis(discrete = TRUE, option = "D")+
  facet_grid(. ~ variable,  scales='free', space = "free")+
  xlab("")+ggtitle("Channels")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("")+
  theme(
    panel.grid.major = element_line(colour = "gray", 
                                    linetype = "dotted"),
    panel.background = element_rect(fill = "white", 
                                    colour="black"),
    strip.text.x = element_text(size=18, 
                                hjust=0.5, 
                                vjust=0.5,
                                face="bold", lineheight = 0.5),
    strip.text.y = element_text(size=18, 
                                hjust=0.5, 
                                vjust=0.5,
                                face="bold", lineheight = 0.5),
    strip.background = element_rect(colour="black", fill="gray98"),
    axis.text=element_text(size=18, face="bold", colour="gray24"),
    
    axis.text.x  = element_blank(),
    #axis.text.x  = element_text(angle=0, vjust=0.5, size=18,colour="white"),
    axis.title=element_text(size=18,face="bold"),
    plot.title = element_text(size = 18, colour = "black", face="bold", hjust=0.5),
    legend.position ="none", 
    panel.spacing = unit(2, "lines")
  )+guides(color=FALSE)+scale_y_continuous(label=scientific_10)+ylim(c(0,max(violin1$value)))

ggsave("../results/San Francisco/plot_San_Francisco_Class.png", pl, width=10, height =8.5, scale=0.8)


#<>----------------------------------------------------------------------------------



#<>---------------------------------------------------------------------------------------------------------------------<>#
# Beginning of Classification
#<>---------------------------------------------------------------------------------------------------------------------<>#

#<>---------------------------------------------------------------------------------------------------------------------<>#
#<>---------------------------------------------------------------------------------------------------------------------<>#
#########
## KNN ##
#########
# k parameter -> number of neighbours
saidafin<-c()
matrizesconfusao <- list() #Classification Matrix

for(k in 1:25){
	print(k)
	saida3<-kNNII(treino, teste, cl=verdadeiro, viz = k)
	c.knn3<-table(verdadeiro_teste,saida3)
	matrizesconfusao[[k]] <- c.knn3
	kapa <- kappaclas(c.knn3)
	diag<-sum(diag(c.knn3))
	(paJ1<-diag/sum(c.knn3)*100)
	saidafin<-rbind(saidafin,c(k,paJ1,kapa$kappa,kapa$VarianciaKappa))
}

# Save Classification Matrix
save(matrizesconfusao,file="../results/San Francisco/Resultado_1-knn_Matrizes_de_confusao.Rdata")

# Save partial results
colnames(saidafin)<-c("k","Accuracy","kappa", "var.kappa")
saidafin1 = saidafin
lresults[[1]] <- saidafin
save(lresults, file = "../results/San Francisco/classification_results_San_Francisco.Rdata")
#<>---------------------------------------------------------------------------------------------------------------------<>#

# Sensitivity analysis plot
pdf("../results/San Francisco/Resultado_1-knn_sensitivity_analysis.pdf", width=7, height=5.5)

minimo = min(saidafin1[,2])-10;
maximo = ifelse(max(saidafin1[,2])+20>100,100,max(saidafin1[,2])+10)
plot(saidafin1[,1],saidafin1[,2], type = "n", ylim=c(minimo,maximo),
     xlab = "k parameter", ylab="Accuracy",pch=20)

abline(h=max(saidafin1[,2]),col="lightgray",lty=2)
abline(v=saidafin1[which.max(saidafin1[,2]),1],col="lightgray",lty=2)

points(saidafin1[,1],saidafin1[,2],pch=20,col="black",type="b")

points(saidafin1[which.max(saidafin1[,2]),1],max(saidafin1[,2]), pch=20, col="green")
points(saidafin1[which.max(saidafin1[,2]),1],max(saidafin1[,2]), pch=21, col="black")

text(saidafin1[which.max(saidafin1[,2]),1],max(saidafin1[,2])+0.5, paste(round(max(saidafin1[,2]),digits=2),"%") )
dev.off()
#
#--------------------------------------------------------------------------- #
#<>---------------------------------------------------------------------------------------------------------------------<>#

#<>---------------------------------------------------------------------------------------------------------------------<>#
#<>---------------------------------------------------------------------------------------------------------------------<>#
#######################
## Kernel radial KNN ##
#######################

# k parameter -> number of neighbours
a<-c(0.0001,0.001,0.01,0.1,1,10,100,1000) # sigma parameter variation

saidafin<-c()
saidafin1<-c()
matrizesconfusao <- list() #Classification Matrix
matrizesconfusao_1 <- list() #Classification Matrix
for(k in 1:25){
  for(l in 1:length(a)){
	print(c(k,l))
		saida4<-kernel_kNN(treino, teste, cl=verdadeiro, viz = k, kernel="radial", sigma=a[l])
		saida4 <- factor(saida4,levels=1:3)
		c.knn4<-table(verdadeiro_teste,saida4)
		matrizesconfusao_1[[l]] <- c.knn3 
		kapa <- kappaclas(c.knn4)
		diag<-sum(diag(c.knn4))
		(paJ2<-diag/sum(c.knn4)*100)
		saidafin1<-rbind(saidafin1,c(k,a[l],paJ2,kapa$kappa,kapa$VarianciaKappa))
		  }
		matrizesconfusao[[k]] <- matrizesconfusao_1
}

# Save Classification Matrix
save(matrizesconfusao,file="../results/San Francisco/Resultado_2-Kernel_knn_Matrizes_de_confusao.Rdata")


#--------------------------------------------------------------------------- #
# Sensitivity analysis plot
pdf("../results/San Francisco/Resultado_2-kernel_knn_sensitivity_analysis.pdf", width=7, height=5.5)

minimo = min(saidafin1[,3])-10;
maximo = ifelse(max(saidafin1[,3])+20>100,100,max(saidafin1[,3])+10)
plot(unique(saidafin1[,1]),saidafin1[1:25,3], type = "n", ylim=c(minimo,maximo),
     xlab = "k parameter", ylab="Accuracy",pch=20)

for(i in 1:length(a)){

subconj = subset(saidafin1,saidafin1[,2]==a[i]) 

points(subconj[,1],subconj[,3],pch=20,col=i,type="b")

}
legend("bottomright",legend=unique(a),col=1:length(a),pch=20)


points(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3]), pch=20, col="green")
points(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3]), pch=21, col="black")

text(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3])+0.5, paste(round(max(saidafin1[,3]),digits=2),"%") )


dev.off()
#--------------------------------------------------------------------------- #



saidafin<-saidafin1
or<-order(saidafin[,3], decreasing = T)
saidafin<-saidafin[or,]

#<>-------------------------------------------<>#
####     Taking the highest ranks of each k #####
saidafin2<-saidafin[1,]
for(i in 2:dim(saidafin)[1]){
if(saidafin[i,1]==saidafin[(i-1),1]){}else{ saidafin2<-rbind(saidafin2,saidafin[i,]) }
}


or2<-order(as.numeric(saidafin2[,1]), decreasing = F)
saidafin2<-saidafin2[or2,]

saidafin3<-saidafin2[1,]
for(i in 2:dim(saidafin2)[1]){
  if(saidafin2[i,1]==saidafin2[(i-1),1]){}else{ saidafin3<-rbind(saidafin3,saidafin2[i,]) }
}
#<>----------------------------------------------<>#

# Save partial results
colnames(saidafin3)<-c("k","sigma","Accuracy","kappa", "var.kappa")
saidafin3
lresults[[2]] <- saidafin3
save(lresults, file = "../results/San Francisco/classification_results_San_Francisco.Rdata")
#<>---------------------------------------------------------------------------------------------------------------------<>#

#<>---------------------------------------------------------------------------------------------------------------------<>#
###############
## Fuzzy KNN ##
###############

# k parameter -> number of neighbours
a<-c(2,3,4,5) # m parameter variation (fuzzy)


saidafin<-c()
saidafin1<-c()
matrizesconfusao <- list() #Classification Matrix
matrizesconfusao_1 <- list() #Classification Matrix

for(k in 1:25){
  for(l in 1:length(a)){
    print(c(k,a[l]))
    ## Fuzzy KNN ##
    saida6<-fuzzy_kNN(treino, teste, cl=verdadeiro, viz = k, m=a[l]) 
    colnames(saida6)<-c("1","2","3")
    saida7<-colnames(saida6)[apply(saida6,1,which.max)]
	saida7 <- factor(as.numeric(saida7),levels=1:3)
    c.knn6<-table(verdadeiro_teste,saida7)
	matrizesconfusao_1[[l]] <- c.knn6
	kapa <- kappaclas(c.knn6) 
    diag<-sum(diag(c.knn6))
    (paJ4<-diag/sum(c.knn6)*100)
    saidafin1<-rbind(saidafin1,c(k,a[l],paJ4,kapa$kappa,kapa$VarianciaKappa))
  }
matrizesconfusao[[k]] <- matrizesconfusao_1
}

# Save Classification Matrix
save(matrizesconfusao,file="../results/San Francisco/Resultado_3-Fuzzy_knn_Matrizes_de_confusao.Rdata")

#--------------------------------------------------------------------------- #
# Sensitivity analysis plot
pdf("../results/San Francisco/Resultado_3-fuzzy_knn_sensitivity_analysis.pdf", width=7, height=5.5)

minimo = min(saidafin1[,3])-10;
maximo = ifelse(max(saidafin1[,3])+20>100,100,max(saidafin1[,3])+10)
plot(unique(saidafin1[,1]),saidafin1[1:25,3], type = "n", ylim=c(minimo,maximo),
     xlab = "k parameter", ylab="Accuracy",pch=20)


for(i in 1:length(a)){

subconj = subset(saidafin1,saidafin1[,2]==a[i]) 

points(subconj[,1],subconj[,3],pch=20,col=i,type="b")

}
legend("bottomright",legend=unique(a),col=1:length(a),pch=20)


points(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3]), pch=20, col="green")
points(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3]), pch=21, col="black")

text(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3])+0.5, paste(round(max(saidafin1[,3]),digits=2),"%") )


dev.off()
#--------------------------------------------------------------------------- #



saidafin<-saidafin1
or<-order(saidafin[,3], decreasing = T)
saidafin<-saidafin[or,]

#<>-------------------------------------------<>#
####     Taking the highest ranks of each k #####

saidafin2<-saidafin[1,]
for(i in 2:dim(saidafin)[1]){
  if(saidafin[i,1]==saidafin[(i-1),1]){}else{ saidafin2<-rbind(saidafin2,saidafin[i,]) }
}


or2<-order(as.numeric(saidafin2[,1]), decreasing = F)
saidafin2<-saidafin2[or2,]

saidafin3<-saidafin2[1,]
for(i in 2:dim(saidafin2)[1]){
  if(saidafin2[i,1]==saidafin2[(i-1),1]){}else{ saidafin3<-rbind(saidafin3,saidafin2[i,]) }
}

colnames(saidafin3)<-c("k","sigma","Accuracy","kappa", "var.kappa")
saidafin3
lresults[[3]] <- saidafin3
save(lresults, file = "../results/San Francisco/classification_results_San_Francisco.Rdata")
#<>---------------------------------------------------------------------------------------------------------------------<>#

#<>---------------------------------------------------------------------------------------------------------------------<>#
######################
## Kernel Fuzzy KNN ##
######################

a<-c(0.001,0.01,0.1,1,2)

saidafin<-c()
saidafin1<-c()
matrizesconfusao <- list() #Classification Matrix
matrizesconfusao_1 <- list() #Classification Matrix


for(k in 1:25){
  for(l in 1:5){
    print(c(k,a[l]))
    
    saida8<-kernel_f_kNN(treino, teste, cl=verdadeiro, viz = k, kernel="radial",sigma=a[l])
    colnames(saida8)<-c("1","2","3")
    saida9<-colnames(saida8)[apply(saida8,1,which.max)]
    c.knn8<-table(verdadeiro_teste,saida9)
	matrizesconfusao_1[[l]] <- c.knn8
	kapa <- kappaclas(c.knn8) 
    diag<-sum(diag(c.knn8))
    (paJ5<-diag/sum(c.knn8)*100)
    saidafin1<-rbind(saidafin1,c(k,a[l],paJ5,kapa$kappa,kapa$VarianciaKappa))
  }
matrizesconfusao[[k]] <- matrizesconfusao_1
}

# Save Classification Matrix
save(matrizesconfusao,file="../results/San Francisco/Resultado_4-kernel_Fuzzy_knn_Matrizes_de_confusao.Rdata")

#--------------------------------------------------------------------------- #
# Sensitivity analysis plot
pdf("../results/San Francisco/Resultado_4-kernel_fuzzy_knn_sensitivity_analysis.pdf", width=7, height=5.5)

minimo = min(saidafin1[,3])-10;
maximo = ifelse(max(saidafin1[,3])+20>100,100,max(saidafin1[,3])+10)
plot(unique(saidafin1[,1]),saidafin1[1:25,3], type = "n", ylim=c(minimo,maximo),
     xlab = "k parameter", ylab="Accuracy",pch=20)

for(i in 1:length(a)){

subconj = subset(saidafin1,saidafin1[,2]==a[i]) 

points(subconj[,1],subconj[,3],pch=20,col=i,type="b")

}
legend("bottomright",legend=unique(a),col=1:length(a),pch=20)


points(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3]), pch=20, col="green")
points(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3]), pch=21, col="black")

text(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3])+0.5, paste(round(max(saidafin1[,3]),digits=2),"%") )


dev.off()
#--------------------------------------------------------------------------- #


saidafin<-saidafin1
or<-order(saidafin[,3], decreasing = T)
saidafin<-saidafin[or,]

#<>-------------------------------------------<>#
####     Taking the highest ranks of each k #####
saidafin2<-saidafin[1,]
for(i in 2:dim(saidafin)[1]){
  if(saidafin[i,1]==saidafin[(i-1),1]){}else{ saidafin2<-rbind(saidafin2,saidafin[i,]) }
}


or2<-order(as.numeric(saidafin2[,1]), decreasing = F)
saidafin2<-saidafin2[or2,]

saidafin3<-saidafin2[1,]
for(i in 2:dim(saidafin2)[1]){
  if(saidafin2[i,1]==saidafin2[(i-1),1]){}else{ saidafin3<-rbind(saidafin3,saidafin2[i,]) }
}

colnames(saidafin3)<-c("k","sigma","Accuracy","kappa", "var.kappa")
saidafin3
lresults[[4]] <- saidafin3
save(lresults, file = "../results/San Francisco/classification_results_San_Francisco.Rdata")
#<>---------------------------------------------------------------------------------------------------------------------<>#



#<>---------------------------------------------------------------------------------------------------------------------<>#
####################
##	 Naive Bayes  ##
####################

train<-treino<-data.frame(dados)
cl<-verdadeiro<-Verd
test<-teste<-data.frame(dadostest)
verdadeiro_teste<-Verdtest

	dados_treino <- cbind(treino,verdadeiro)

	fit_NB <- naiveBayes(verdadeiro~., data=dados_treino)
	predictions <- predict(fit_NB,newdata=test, type = "class")
	predictions_table <- table( factor(verdadeiro,levels=c(1,2,3)), factor(predictions,levels=c(1,2,3)) )
	kapa <- kappaclas(predictions_table)
	diag<-sum(diag(predictions_table))
	(paJ4<-diag/sum(predictions_table)*100)
	saidafin<-c(paJ4,kapa$kappa,kapa$VarianciaKappa)

save(predictions_table,file="../results/San Francisco/Resultado_5-Naive-Bayes_Matriz_de_confusao.Rdata")

names(saidafin) <- c("Accuracy","kappa", "var.kappa")
lresults[[5]] <- saidafin
save(lresults, file = "../results/San Francisco/classification_results_San_Francisco.Rdata")
#<>---------------------------------------------------------------------------------------------------------------------<>#

#<>---------------------------------------------------------------------------------------------------------------------<>#
#########
## SVM ##
#########
#Support Vector Machine
#install.packages("e1071")
library(e1071)

set.seed(04011991)
a<-c(1,10,100,1000)
b<-c(3,30,300,3000)

saidafin<-c()
saidafin1<-c()
matrizesconfusao <- list() #Classification Matrix
matrizesconfusao_1 <- list() #Classification Matrix


for(k in 1:length(a)){
  for(l in 1:length(b)){
    print(c(a[k],b[l]))
    model<-svm(dados, Verd, type="C-classification",cost = a[k], gamma = b[l],kernel="radial")
    saida8 <- predict(model, dadostest)
    c.knn8<-table(verdadeiro_teste,saida8)
	kapa <- kappaclas( c.knn8)
	matrizesconfusao_1[[l]] <- c.knn8
    diag<-sum(diag(c.knn8))
    (paJ5<-diag/sum(c.knn8)*100)
    saidafin1<-rbind(saidafin1,c(a[k],b[l],paJ5,kapa$kappa,kapa$VarianciaKappa))
  }
matrizesconfusao[[k]] <- matrizesconfusao_1
}

# Save Classification Matrix
save(matrizesconfusao,file="../results/San Francisco/Resultado_6-SVM_Matrizes_de_confusao.Rdata")


#--------------------------------------------------------------------------- #
# Sensitivity analysis plot
pdf("../results/San Francisco/Resultado_6-SVM_sensitivity_analysis.pdf", width=7, height=5.5)

minimo = min(saidafin1[,3])-10;
maximo = ifelse(max(saidafin1[,3])+20>100,100,max(saidafin1[,3])+10)
plot(unique(saidafin1[,1]),unique(saidafin1[,1]), type = "n", ylim=c(minimo,maximo),
     xlab = "k parameter", ylab="Accuracy",pch=20)

for(i in 1:length(b)){

subconj = subset(saidafin1,saidafin1[,2]==b[i]) 

points(subconj[,1],subconj[,3],pch=20,col=i,type="b")

}

legend("bottomright",legend=unique(a),col=1:length(a),pch=20)


points(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3]), pch=20, col="green")
points(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3]), pch=21, col="black")

text(saidafin1[which.max(saidafin1[,3]),1],max(saidafin1[,3])+0.5, paste(round(max(saidafin1[,3]),digits=2),"%") )


dev.off()
#--------------------------------------------------------------------------- #


saidafin<-saidafin1
or<-order(saidafin[,3], decreasing = T)
saidafin<-saidafin[or,]

#<>-------------------------------------------<>#
####     Taking the highest ranks of each k #####
saidafin2<-saidafin[1,]
for(i in 2:dim(saidafin)[1]){
  if(saidafin[i,1]==saidafin[(i-1),1]){}else{ saidafin2<-rbind(saidafin2,saidafin[i,]) }
}

or2<-order(as.numeric(saidafin2[,1]), decreasing = F)
saidafin2<-saidafin2[or2,]

saidafin3<-saidafin2[1,]
for(i in 2:dim(saidafin2)[1]){
  if(saidafin2[i,1]==saidafin2[(i-1),1]){}else{ saidafin3<-rbind(saidafin3,saidafin2[i,]) }
}

#<>----------------------------------------------<>#
# Save Classification Matrix
colnames(saidafin3)<-c("constant","learning parameter","% accuracy","kappa", "var.kappa")
saidafin3
lresults[[6]] <- saidafin3
save(lresults, file = "../results/San Francisco/classification_results_San_Francisco.Rdata")
#<>---------------------------------------------------------------------------------------------------------------------<>#

#<>---------------------------------------------------------------------------------------------------------------------<>#
######################
## Xboost Classifier ##
######################
library(xgboost)

train<-treino<-data.frame(dados)
cl<-verdadeiro<-Verd
test<-teste<-data.frame(dadostest)
verdadeiro_teste<-Verdtest


dados_treino <- cbind(train,verdadeiro)
dadostestec <- cbind(test,verdadeiro_teste)

set.seed(978)
data_train <- dados_treino[sample(nrow(dados_treino)),]
data_test <- dadostestec[sample(nrow(dadostestec)),]


x_treino <- data_train[,-4] %>% as.matrix()
y_treino <- data_train[,4]  %>% as.numeric()-1
x_teste <- data_test[,-4] %>% as.matrix()
y_teste <- data_test[,4] %>% as.numeric()-1


etaa <- c(0.3,0.2,0.5,0.7)
lam <- c(0.6,0.9,1.2,2.5,3.5)

results_xgb <-c()

for(i in 1:length(etaa)){
    for(n in 1:length(lam)){
      
      xgb.train <- xgb.DMatrix(data=x_treino,label=y_treino)
      xgb.test  <-  xgb.DMatrix(data=x_teste,label=y_teste)
      
      
      num_class = 3
      params = list( # usando os parâmetros de defoult primiero 
        booster="gbtree",
        eta=etaa[i],
        max_depth=6,
        gamma=0, # funciona bem com valores menores de max_depht
        subsample=0.1,
        colsample_bytree=1,
        #min_child_weight=1,
        lambda=lam[n],
        objective="multi:softmax",  
        eval_metric="mlogloss",   
        num_class=num_class
      )
      
      
      xgbcv <- xgb.cv(params = params, 
                      data = xgb.train, 
                      nrounds = 200, 
                      nthread = 25,
                      nfold = 2, 
                      showsd = T, 
                      stratified = T, 
                      # print.every.n = 100,
                      early.stop.round = 70,
                      watchlist = list(val=xgb.test,train=xgb.train),
                      seed=956,
                      maximize = F)
      
      
      results_xgb <- rbind(results_xgb,c(etaa[i],lam[n],xgbcv$best_iteration,xgbcv$evaluation_log[xgbcv$best_iteration,-1])) 
    }
  }

results_xgb <- results_xgb %>% data.frame()

results_xgb3 <- results_xgb[,c(-5,-7)] %>% 
  mutate(Dif=as.numeric(train_mlogloss_mean)-as.numeric(test_mlogloss_mean))

results_xgb1 <- results_xgb[,1:3] %>% data.frame()
colnames(results_xgb1) <- c("Eta", "Lambda","n_rounds") #"train_mean","train_sd", "test_mean","test_sd")

matrizesconfusao <- list()
results_fim <- c()
for(k in 1:dim(results_xgb1)[1]){

num_class = 3

 params1 = list( 
  booster="gbtree",
  eta=as.numeric(results_xgb1[k,1]),
  max_depth=6,
  subsample=0.5,
  colsample_bytree=1,
  gamma=0,
  lambda=as.numeric(results_xgb1[k,2]),
  objective="multi:softmax",
  eval_metric="mlogloss",   
  num_class=num_class
)



fit <- xgb.train (params = params1, 
                  data = xgb.train, 
                  nrounds = as.numeric(results_xgb1[k,3]), 
                  watchlist = list(val=xgb.test,train=xgb.train),
                  maximize = F) 


fit$evaluation_log
xgbpred<- predict(fit,xgb.test)
xgbpred1 <- factor(xgbpred, levels = 0:2)
classes <- factor(y_teste, levels = 0:2)

matriz <- table(classes, xgbpred1) 
accuracy <- (confusionMatrix(matriz)$overall[1]) *100
kapa <- kappaclas(matriz)$kappa
varkapa <- kappaclas(matriz)$VarianciaKappa
print(matriz)
matrizesconfusao[[k]] <- matriz

results_fim <- rbind(results_fim,c(results_xgb1[k,1],results_xgb1[k,2],results_xgb1[k,3],accuracy,kapa,varkapa))

}

colnames(results_fim) <- c("eta", "lambda", "n_tree","% accuracy","kappa", "var.kappa")
results_fim
lresults[[7]] <- results_fim

save(lresults, file = "../results/San Francisco/classification_results_San_Francisco.Rdata")
save(matrizesconfusao,file="../results/San Francisco/Resultado_7-Xgboost_Matriz_de_confusao.Rdata")
#<>---------------------------------------------------------------------------------------------------------------------<>#




#<>---------------------------------------------------------------------------------------------------------------------<>#
######################################
## Kullback-Leibler Distance PolSAR ##
######################################

# Considering that test samples have the same window size in all regions

#____________________________________________________#
#reading the training datasets

#<>------------------------------------------------<>
image_PolSAR <- array(NA,dim = c(dim(image)[1:2],6))
image_PolSAR[,,1] <- image[,,1]+0i
image_PolSAR[,,2] <- image[,,2]+0i
image_PolSAR[,,3] <- image[,,3]+0i
image_PolSAR[,,4] <- complex(real = image[,,4], imaginary = image[,,7])
image_PolSAR[,,5] <- complex(real = image[,,5], imaginary = image[,,8])
image_PolSAR[,,6] <- complex(real = image[,,6], imaginary = image[,,9])

WaterdataPol <- sapply(1:6, function(x) c(  image_PolSAR[y.minotr:y.maxotr , x.minotr:x.maxotr , x]  ) )
VegetationdataPol <- sapply(1:6, function(x) c( image_PolSAR[y.minftr:y.maxftr , x.minftr:x.maxftr , x] ) )
LowUrbandataPol <- sapply(1:6, function(x) c( image_PolSAR[y.minldtr:y.maxldtr , x.minldtr:x.maxldtr , x] ) )

dadosPol <-rbind(WaterdataPol,VegetationdataPol,LowUrbandataPol) #datasets
VerdPol <- factor(c(rep("1",dim(WaterdataPol)[1]), rep("2",dim(VegetationdataPol)[1]),rep("3",dim(LowUrbandataPol)[1]))) #Converting character variable to categorical
#____________________________________________________#
#reading the test datasets
WaterdataPol <- sapply(1:6, function(x) c(  image_PolSAR[y.minote:y.maxote , x.minote:x.maxote , x]  ) )
VegetationdataPol <- sapply(1:6, function(x) c( image_PolSAR[y.minfte:y.maxfte , x.minf:xte.maxfte , x] ) )
LowUrbandataPol <- sapply(1:6, function(x) c( image_PolSAR[y.minldte:y.maxldte , x.minldte:x.maxldte , x] ) )

dadostestPol <-rbind(WaterdataPol,VegetationdataPol,LowUrbandataPol) #datasets
VerdtestPol <- factor(c(rep("1",dim(WaterdataPol)[1]), rep("2",dim(VegetationdataPol)[1]), rep("3",dim(LowUrbandataPol)[1]))) #Converting character variable to categorical
#___________________________________________________#

trainPol <- treinoPol <- data.frame(dadosPol)
clPol <- verdadeiroPol <- VerdPol
testPol <- testePol <- data.frame(dadostestPol)
verdadeiro_testePol <- VerdtestPol


# ------------------------------------------ #
#  Finding the parameters of the regions     #
# ------------------------------------------ #
saida1<-gerarpar(trainPol,clPol)

L=4

# Format data.frame -> type = 2; Format array => type = 1
resultado <- kl.distance(test=testPol,cltest = verdadeiro_testePol,matrizes_par = saida1$matriz,qg = saida1$ngroups,nlooks=L,nlinhas=length(y.minote:y.maxote),ncolunas=length(x.minote:x.maxote),k=1,type=2)

	resultado_table <- table( factor(verdadeiro_testePol,levels=c(1,2,3)), factor(resultado,levels=c(1,2,3)) )
	kapa <- kappaclas(resultado_table)
	diag<-sum(diag(resultado_table))
	(paJ4<-diag/sum(resultado_table)*100)
	saidafin<-c(paJ4,kapa$kappa,kapa$VarianciaKappa)

save(resultado_table,file="../results/San Francisco/Resultado_8-Kullback-Leibler_Matriz_de_confusao.Rdata")

names(saidafin) <- c("Accuracy","kappa", "var.kappa")
lresults[[8]] <- saidafin
save(lresults, file = "../results/San Francisco/classification_results_San_Francisco.Rdata")


#<>---------------------------------------------------------------------------------------------------------------------<>#
#<>---------------------------------------------------------------------------------------------------------------------<>#
#<>---------------------------------------------------------------------------------------------------------------------<>#
# image full classification
# load(file = "../results/San Francisco/classification_results_San_Francisco.Rdata")
#<>---------------------------------------------------------------------------------------------------------------------<>#
# Lendo os resultados
load(file = "../results/San Francisco/classification_results_San_Francisco.Rdata")

library(ggplot2)
library(reshape2)

# the "test" data is now the full image
teste <- data.frame(c(Re(image[,,1])),c(Re(image[,,2])),c(Re(image[,,3])) )

#<>-------<>
# KNN
#<>-------<>
# Determining the best k

k = lresults[[1]][which.max(lresults[[1]][,3]),][1] #melhor k

# KNN
saida3<-kNNII(treino, teste, cl=verdadeiro, viz = k)

graf.class2<-saida3

saida_g_m <- matrix(as.numeric(graf.class2),ncol=ncol(image))
longData<-melt(apply((saida_g_m), 2, rev))

longData[,3] <- as.factor(longData[,3])

png("../results/San Francisco/1 - kNN - San Francisco.png",width=500, height=400)
  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_manual(labels = c("Ocean", "Forest","Urban"), values = c('#440144', '#21908c', '#FDE725'))+ 
    theme_linedraw() +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
dev.off()

#<>-----------<>
# Kernel KNN
#<>-----------<>
# Determining the best k and sigma
# 23 - 5
k1=lresults[[2]][which.max(lresults[[2]][,3]),][1] # k Best
a1=lresults[[2]][which.max(lresults[[2]][,3]),][2] # value referring to the best sigma
#kernel KNN
saida4<-kernel_kNN(treino, teste, cl=verdadeiro, viz = k1, kernel="radial", sigma=a1) 

graf.class2<-saida4

saida_g_m <- matrix(as.numeric(graf.class2),ncol=ncol(image))
longData<-melt(apply((saida_g_m), 2, rev))
longData[,3] <- as.factor(longData[,3])

png("../results/San Francisco/2 - Kernel kNN - San Francisco.png",width=500, height=400)
  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_manual(labels = c("Ocean", "Forest","Urban"), values = c('#440144', '#21908c', '#FDE725'))+
    theme_linedraw() +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
dev.off()


#<>-----------<>
# Fuzzy KNN
#<>-----------<>
# Determining the best k and m
# 14 - 4
k1=lresults[[3]][which.max(lresults[[3]][,3]),][1] # k Best
m1=lresults[[3]][which.max(lresults[[3]][,3]),][2] # m referring to the best k
#fuzzy KNN
saida6<-fuzzy_kNN(treino, teste, cl=verdadeiro, viz = k1, m=m1) 
#colnames(saida6)<-c("1","2","3")
#saida7<-colnames(saida6)[apply(saida6,1,which.max)]

saida62<-array(NA,dim=c(3,dim(image)[1]*dim(image)[2])) #dim => numero de observacoes x numero de grupos

for(i in 1:dim(saida6)[2]){

saida62[i] <- saida6[[i]]
print(i);

}

saida62 <- t(saida62)
colnames(saida62)<-c("1","2","3")
saida7<-colnames(saida62)[apply(saida62,1,which.max)]


graf.class2<-saida7

saida_g_m <- matrix(as.numeric(graf.class2),ncol=ncol(image))
longData<-melt(apply((saida_g_m), 2, rev))

longData[,3] <- as.factor(longData[,3])

png("../results/San Francisco/3 - Fuzzy kNN - San Francisco.png",width=500, height=400)

  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_manual(labels = c("Ocean", "Forest","Urban"), values = c('#440144', '#21908c', '#FDE725'))+
    theme_linedraw() +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())

dev.off()


#<>--------------------<>
# Fuzzy Kernel KNN
#<>--------------------<>
# Determining the best k and sigma
k1=as.numeric(lresults[[4]][which.max(lresults[[4]][,3]),][1]) #melhor k
a1=as.numeric(lresults[[4]][which.max(lresults[[4]][,3]),][2]) # sigma referente ao melhor k
#fuzzy kernel KNN
saida8<-kernel_f_kNN(treino, teste, cl=verdadeiro, viz = k1, kernel="radial",sigma=a1)

colnames(saida8) <- c("1","2","3")
saida9 <- colnames(saida8)[apply(saida8,1,which.max)]

graf.class2 <- saida9

saida_g_m <- matrix(as.numeric(graf.class2),ncol=ncol(image))
longData<-melt(apply((saida_g_m), 2, rev))

longData[,3] <- as.factor(longData[,3])

png("../results/San Francisco/4 - Fuzzy Kernel kNN - San Francisco.png",width=500, height=400)
  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_manual(labels = c("Ocean", "Forest","Urban"), values = c('#440144', '#21908c', '#FDE725'))+
    theme_linedraw() +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
dev.off()

#<>----------------------<>
# Naive Bayes
#<>----------------------<>
colnames(treino) <- c("V1","V2","V3")
dados_treino <- cbind(treino,verdadeiro)
colnames(teste) <- c("V1","V2","V3")

fit_NB <- naiveBayes(verdadeiro~., data=dados_treino)
predictions <- predict(fit_NB,newdata=teste)

graf.class2 <- predictions

saida_g_m <- matrix(as.numeric(graf.class2),ncol=ncol(image))
longData<-melt(apply((saida_g_m), 2, rev))

longData[,3] <- as.factor(longData[,3])

png("../results/San Francisco/5 - Naive Bayes - San Francisco.png",width=500, height=400)

  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_manual(labels = c("Ocean", "Forest","Urban"), values = c('#440144', '#21908c', '#FDE725'))+ 
    theme_linedraw() +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())

dev.off()


#<>----------------------<>
# Support Vector Machine
#<>----------------------<>
library(e1071)
library(parallel)
# Reading results from the best parameters
a1=lresults[[6]][which.max(lresults[[6]][,3]),][1] # cost best
b1=lresults[[6]][which.max(lresults[[6]][,3]),][2] # gamma best
# SVM
model<-svm(dados, Verd, type="C-classification",cost = a1, gamma = b1,kernel="radial")
saida8 <- predict(model, teste)

graf.class2 <- saida8

saida_g_m <- matrix(as.numeric(graf.class2),ncol=ncol(image))
longData<-melt(apply((saida_g_m), 2, rev))

longData[,3] <- as.factor(longData[,3])

png("../results/San Francisco/6 - SVM - San Francisco.png",width=500, height=400)
  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_manual(labels = c("Ocean", "Forest","Urban"), values = c('#440144', '#21908c', '#FDE725'))+ 
    theme_linedraw() +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
dev.off()

#<>----------------------<>
# Xboost Classifier
#<>----------------------<>
library(xgboost)

teste <- data.frame(c(as.numeric(Re(image[,,1]))),c(as.numeric(Re(image[,,2]))),
                    c(as.numeric(Re(image[,,3]))) )
colnames(teste) <- c("X1","X2","X3")
teste <- teste %>% as.matrix()
xgb.train1 <- xgb.DMatrix(data=x_treino,label=y_treino)
xgb.test1 <-  xgb.DMatrix(data=teste)


eta=lresults[[7]][which.max(lresults[[7]][,4]),][1] # eta
lam=lresults[[7]][which.max(lresults[[7]][,4]),][2] # lambda
ntree=lresults[[7]][which.max(lresults[[7]][,4]),][3] # n_tree

params2 = list( 
  booster="gbtree",
  eta=eta,
  max_depth=6,
  gamma=0,
  subsample=0.5,
  colsample_bytree=1,
  lambda=lam,
  objective="multi:softmax",  
  eval_metric="mlogloss",   
  num_class=num_class
)



fit2 <- xgb.train (params = params2, 
                  data = xgb.train1, 
                  nrounds = as.numeric(ntree), 
                  maximize = F) 


xgbpred2 <- predict (fit,newdata=xgb.test1)

graf.class2 <- xgbpred2+1

saida_g_m <- matrix(as.numeric(graf.class2),ncol=ncol(image))
longData<-melt(apply((saida_g_m), 2, rev))
longData[,3] <- as.factor(longData[,3])


png("../results/San Francisco/7 - Xgboost - San Francisco.png",width=500, height=400)
  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_manual(labels = c("Ocean", "Forest","Urban"), values = c('#440144', '#21908c', '#FDE725'))+
    theme_linedraw() +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
dev.off()


#<>----------------------------------<>
# Kullback-Leibler Distance PolSAR
#<>----------------------------------<>


#<>------------------------------------------------<>
image_PolSAR <- array(NA,dim = c(dim(image)[1:2],6))
image_PolSAR[,,1] <- image[,,1]+0i
image_PolSAR[,,2] <- image[,,2]+0i
image_PolSAR[,,3] <- image[,,3]+0i
image_PolSAR[,,4] <- complex(real = image[,,4], imaginary = image[,,7])
image_PolSAR[,,5] <- complex(real = image[,,5], imaginary = image[,,8])
image_PolSAR[,,6] <- complex(real = image[,,6], imaginary = image[,,9])

WaterdataPol <- sapply(1:6, function(x) c(  image_PolSAR[y.minotr:y.maxotr , x.minotr:x.maxotr , x]  ) )
VegetationdataPol <- sapply(1:6, function(x) c( image_PolSAR[y.minftr:y.maxftr , x.minftr:x.maxftr , x] ) )
LowUrbandataPol <- sapply(1:6, function(x) c( image_PolSAR[y.minldtr:y.maxldtr , x.minldtr:x.maxldtr , x] ) )

dadosPol <-rbind(WaterdataPol,VegetationdataPol,LowUrbandataPol) #datasets
VerdPol <- factor(c(rep("1",dim(WaterdataPol)[1]), rep("2",dim(VegetationdataPol)[1]),rep("3",dim(LowUrbandataPol)[1]))) #Converting character variable to categorical

trainPol <- treinoPol <- data.frame(dadosPol)
clPol <- verdadeiroPol <- VerdPol
saida1<-gerarpar(trainPol,clPol)

L=4

# Format data.frame -> type = 2; Format array => type = 1
resultado <- kl.distance(test=image_PolSAR,matrizes_par = saida1$matriz,qg = saida1$ngroups,nlooks=L,nlinhas=dim(image_PolSAR)[1],ncolunas=dim(image_PolSAR)[2],k=1,type=1)

saida_g_m <- matrix(as.numeric(resultado),ncol=ncol(image))
longData<-melt(apply((saida_g_m), 2, rev))

longData[,3] <- as.factor(longData[,3])

png("../results/San Francisco/8 - KL distance - San Francisco.png",width=500, height=400)
  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_manual(labels = c("Ocean", "Forest","Urban"), values = c('#440144', '#21908c', '#FDE725'))+
    theme_linedraw() +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
dev.off()

