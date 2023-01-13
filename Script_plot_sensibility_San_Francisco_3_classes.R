
#-------------------------------------------------------------
#     Plot sensibility ---------------------------------------
#-------------------------------------------------------------

library(ggplot2)
library(plyr)
library(tidyverse)


load("classification_results_San_Francisco.Rdata")

lresults[[1]]

minimo <- min(lresults[[1]][,2])-10
maximo = ifelse(max(lresults[[1]][,2])+20>100,100,max(lresults[[1]][,2])+10)

pdf("Resultado_1-knn_sensitivity_analysis.pdf", width=7, height=5.5)

 ggplot() + geom_line(aes(y = lresults[[1]][,2], x = (lresults[[1]][,1])), size=1.5,
                           data = data.frame(lresults[[1]]) , stat="identity")+#,linetype="dashed") +
                            ylim(minimo,maximo)+ theme_bw() +
                            xlab("k parameter")+
                            ylab("Accuracy")+
theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),
      axis.text.x  = element_text(colour="black",size=16),
      axis.text.y  = element_text(colour="black",size=16),
      strip.background = element_rect(colour="black", fill="gray98"),
      strip.text.x = element_text(size=18, 
                                  hjust=0.5, 
                                  vjust=0.5,
                                  face="bold", lineheight = 0.5),
      strip.text.y = element_text(size=18, 
                                  hjust=0.5, 
                                  vjust=0.5,
                                  face="bold", lineheight = 0.5),
      axis.title=element_text(size=16,face="bold")
     
      ) +
  annotate(
    geom = "point",
     x = lresults[[1]][which.max(lresults[[1]][,2]),][1],
     y = lresults[[1]][which.max(lresults[[1]][,2]),][2],
     size = 6,
     colour = "red"
   ) 

dev.off()


minimo <- min(lresults[[2]][,3])-5
maximo = ifelse(max(lresults[[2]][,3])+20>100,100,max(lresults[[2]][,3])+10)

pdf("Resultado_2-Kernel-knn_sensitivity_analysis.pdf", width=7, height=5.5)
ggplot() + geom_line(aes(y = lresults[[2]][,3], x = (lresults[[2]][,1])), size=1.5,
                           data = data.frame(lresults[[2]]) , stat="identity")+#,linetype="dashed") +
  ylim(minimo,maximo)+ theme_bw() +
  xlab("k parameter")+
  ylab("Accuracy")+
theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),
         axis.text.x  = element_text(colour="black",size=16),
         axis.text.y  = element_text(colour="black",size=16),
         strip.background = element_rect(colour="black", fill="gray98"),
         strip.text.x = element_text(size=18, 
                                     hjust=0.5, 
                                     vjust=0.5,
                                     face="bold", lineheight = 0.5),
         strip.text.y = element_text(size=18, 
                                     hjust=0.5, 
                                     vjust=0.5,
                                     face="bold", lineheight = 0.5),
         axis.title=element_text(size=16,face="bold")
         
) +
  
  annotate(
    geom = "point",
    x = lresults[[2]][which.max(lresults[[2]][,3]),][1],
    y = lresults[[2]][which.max(lresults[[2]][,3]),][3],
    size = 7,
    colour = "red"
  ) 
dev.off()



minimo <- min(lresults[[3]][,3])-10
maximo = ifelse(max(lresults[[3]][,3])+20>100,100,max(lresults[[3]][,3])+10)

pdf("Resultado_3-Fuzzy-knn_sensitivity_analysis.pdf", width=7, height=5.5)
ggplot() + geom_line(aes(y = lresults[[3]][,3], x = (lresults[[3]][,1])), size=1.5,
                     data = data.frame(lresults[[3]]) , stat="identity")+#,linetype="dashed") +
  ylim(minimo,maximo)+ theme_bw() +
  xlab("k parameter")+
  ylab("Accuracy")+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),
        axis.text.x  = element_text(colour="black",size=16),
        axis.text.y  = element_text(colour="black",size=16),
        strip.background = element_rect(colour="black", fill="gray98"),
        strip.text.x = element_text(size=18, 
                                    hjust=0.5, 
                                    vjust=0.5,
                                    face="bold", lineheight = 0.5),
        strip.text.y = element_text(size=18, 
                                    hjust=0.5, 
                                    vjust=0.5,
                                    face="bold", lineheight = 0.5),
        axis.title=element_text(size=16,face="bold")
        
  ) +
  
  annotate(
    geom = "point",
    x = lresults[[3]][which.max(lresults[[3]][,3]),][1],
    y = lresults[[3]][which.max(lresults[[3]][,3]),][3],
    size = 7,
    colour = "red"
  ) 
dev.off()




minimo <- min(lresults[[4]][,3])-10
maximo = ifelse(max(lresults[[4]][,3])+20>100,100,max(lresults[[4]][,3])+10)

pdf("Resultado_4-Kernel-Fuzzy-knn_sensitivity_analysis.pdf", width=7, height=5.5)
ggplot() + geom_line(aes(y = lresults[[4]][,3], x = (lresults[[4]][,1])), size=1.5,
                     data = data.frame(lresults[[4]]) , stat="identity")+#,linetype="dashed") +
  ylim(minimo,maximo)+ theme_bw() +
  xlab("k parameter")+
  ylab("Accuracy")+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),
        axis.text.x  = element_text(colour="black",size=16),
        axis.text.y  = element_text(colour="black",size=16),
        strip.background = element_rect(colour="black", fill="gray98"),
        strip.text.x = element_text(size=18, 
                                    hjust=0.5, 
                                    vjust=0.5,
                                    face="bold", lineheight = 0.5),
        strip.text.y = element_text(size=18, 
                                    hjust=0.5, 
                                    vjust=0.5,
                                    face="bold", lineheight = 0.5),
        axis.title=element_text(size=16,face="bold")
        
  ) +
  
  annotate(
    geom = "point",
    x = lresults[[4]][which.max(lresults[[4]][,3]),][1],
    y = lresults[[4]][which.max(lresults[[4]][,3]),][3],
    size = 7,
    colour = "red"
  ) 
dev.off()

minimo <- min(lresults[[6]][,3])-5
maximo = ifelse(max(lresults[[6]][,3])+20>100,100,max(lresults[[6]][,3])+10)

pdf("Resultado_6-SVM_sensitivity_analysis.pdf", width=7, height=5.5)
ggplot() + geom_line(aes(y = lresults[[6]][,3], x = (lresults[[6]][,1])), size=1.5,
                     data = data.frame(lresults[[6]]) , stat="identity")+#,linetype="dashed") +
  ylim(minimo,maximo)+ theme_bw() +
  xlab("k parameter")+
  ylab("Accuracy")+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),
        axis.text.x  = element_text(colour="black",size=16),
        axis.text.y  = element_text(colour="black",size=16),
        strip.background = element_rect(colour="black", fill="gray98"),
        strip.text.x = element_text(size=18, 
                                    hjust=0.5, 
                                    vjust=0.5,
                                    face="bold", lineheight = 0.5),
        strip.text.y = element_text(size=18, 
                                    hjust=0.5, 
                                    vjust=0.5,
                                    face="bold", lineheight = 0.5),
        axis.title=element_text(size=16,face="bold")
        
  ) +
  
  annotate(
    geom = "point",
    x = lresults[[6]][which.max(lresults[[6]][,3]),][1],
    y = lresults[[6]][which.max(lresults[[6]][,3]),][3],
    size = 7,
    colour = "red"
  ) 
dev.off()
  

dados1 <- cbind(unlist(lresults[[7]][,3]), unlist( lresults[[7]][,4])) %>% data.frame()



minimo <- min(dados1[,2])-10
maximo = ifelse(max(dados1[,2])+20>100,100,max(dados1[,2])+10)

pdf("Resultado_8-XGboost_sensitivity_analysis.pdf", width=7, height=5.5)
ggplot() + geom_line(aes(y = dados1[,2], x = dados1[,1]), size=1.5,
                     data = dados1, stat="identity")+#,linetype="dashed") +
  scale_x_continuous(breaks=c(4,14,28,34,48))+
  ylim(minimo,maximo)+ theme_bw() +
  xlab("n Tree")+
  ylab("Accuracy")+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),
        axis.text.x  = element_text(colour="black",size=16),
        axis.text.y  = element_text(colour="black",size=16),
        strip.background = element_rect(colour="black", fill="gray98"),
        strip.text.x = element_text(size=18, 
                                    hjust=0.5, 
                                    vjust=0.5,
                                    face="bold", lineheight = 0.5),
        strip.text.y = element_text(size=18, 
                                    hjust=0.5, 
                                    vjust=0.5,
                                    face="bold", lineheight = 0.5),
        axis.title=element_text(size=16,face="bold")
        
  ) +
  
  annotate(
    geom = "point",
    x = unlist(dados1[which.max(dados1[,2]),][1]),
    y = unlist(dados1[which.max(dados1[,2]),][2]),
    size = 7,
    colour = "red"
  ) 
dev.off()

  

