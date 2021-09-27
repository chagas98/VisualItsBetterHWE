#Bibliotecas utilizadas
library(ggplot2)
library(tidyr)
library(dplyr)
library(colorspace)
library(gmodels)
library(scales)
library(reshape)
library(HardyWeinberg)
library(ggrepel)
library(rlist)


rm(list = ls())

#-----------Variaveis utilizadas no Pipeline-------------------#
File <- readline(prompt="Nome do Arquivo: ")
Dir <- readline(prompt="Local/Diretorio: ")
VarMod <- readline(prompt = "Condicao alterada (nome da coluna): ")
IntMod <- as.numeric(unlist(strsplit((readline(prompt = "Condicao alterada (numero da coluna):")),",")))
Exp <- readline(prompt = "Nome do Experimento: ")


#Zerar variaveis, achar diretorio
work_dir <- "/home/samuel/Documentos/Evolucao/Artigo"   #
setwd(work_dir)
getwd()

#Abrir dados
data0 <- read.csv(
  "Dados.csv",
  stringsAsFactors = F,
  header = T,
  sep = ","
)
data0$Name
#Filtrar pelo experimento
data1 <- dplyr::filter(data0,Name == "Fit")
data1
pre <- data1[,IntMod]

listat<-lapply(pre, `[`)

listat 
#Agrupar e realizar a media
data_mean <- aggregate(data1[, c("Population_Size", 
                                 "Prop_R_Al", "Prop_r_Al",
                                 "Prop_rr_Geno", "Prop_Rr_Geno",
                                 "Prop_RR_Geno")], 
                       listat, mean)
data_mean
#Gerar tabelas
data2 <- data1 %>%
  gather(., key = "Var", value = "Values", Prop_R_Al:Prop_RR_Geno) %>%
  filter(., Var != "Population_Size")
data2
data_mean_gather <- gather(data_mean, key = "Var", value = "Values", Prop_rr_Geno:Prop_RR_Geno)
data_mean_gather
hw_cn <- filter(data2, Var == "Prop_R_Al" | Var == "Prop_r_Al")
hw_cn

#--------------Valores Esperados-------------------#
data_test <- as.data.frame(seq(0, 1, by = 0.01))
#Adicionar teorico p
colnames(data_test) <- "p"
data_test
#Teorico q
data_test$q <- 1 - data_test$p
#Calcilar freq teoricas para genotipos
data_test$AA <- data_test$p ^ 2
data_test$Aa <- 2 * data_test$p * data_test$q
data_test$aa <- data_test$q ^ 2
#Calacular freq total teorica
data_test$total <-
  (data_test$p ^ 2) + (2 * data_test$p * data_test$q) + (data_test$q ^ 2)
#Gerar tabela com valores
data_teo <- gather(data_test, key = "Var", value = "Values", 3:6)
data_teo
#----------------------------------------------------#



###################Calculos#######################
#Calculo do N a partir da freq genotipica
data_stat <- data_mean
data_stat$RRN <-
  round(data_stat$Population_Size * data_stat$Prop_RR_Geno)
data_stat$RrN <-
  round(data_stat$Population_Size * data_stat$Prop_Rr_Geno)
data_stat$rrN <-
  round(data_stat$Population_Size * data_stat$Prop_rr_Geno)

#Define matrix para rodar o HWE
data_stat1 <- data.matrix(data_stat[, c("RRN", "RrN", "rrN")])
colnames(data_stat1) <- c("RR", "Rr", "rr")
data_stat
Results <- HWChisqMat(data_stat1)
Output <- cbind(data_stat, Results$chisqvec, Results$pvalvec)
Output$Resultado <- ifelse(Output$`Results$chisqvec` > 3.84,"Rejeita" ,
                    ifelse(Output$`Results$chisqvec` < 3.84, 'H0', 'Igual'))
write.csv(Output, file = "Resultados_fit.csv")
print(Output)

#----------------------------Graficos------------------------------#
Fig1 <-
  ggplot() +
  geom_point(data1, mapping = aes(x = Prop_R_Al, y = Prop_r_Al, 
                                  color = Generation, size = Tratamento))+
  labs(
    title = "",
    x = "\n Frequência Alélica R (p)",
    y = "Frequência Alélica r (q) \n",
    color = "Geracões"
  )+ 
  #geom_smooth(method = lm, color = "black") +
  theme_bw()
Fig1

Fig2 <- ggplot(hw_cn, aes(x = Generation, y = Values, color = Var), size = 10) +
  geom_point()+
  labs(
    title = "",
    x = "Geracões",
    y = "Frequência Alélica \n",
    color = ""
  ) +
  scale_color_manual(
    labels = c(
      "Freq r (q)",
      "Freq R (p)"
    ),
    values = c(
      "#3399FF",
      "#003366")
  ) +
  theme_bw() +
  geom_smooth(method = lm)+
  facet_wrap(~Tratamento)
Fig2

Fig3 <- ggplot() +
  geom_line(data_teo, mapping = aes(x = p, y = Values, color = Var)) +
  geom_point(data_mean_gather,
             mapping = aes(
               x = Prop_R_Al,
               y = Values,
               shape = Var
             )) +
  scale_color_manual(
    labels = c(
      "q² teórico",
      "2pq teórico",
      "p² teórico",
      "p²+2pq+q²"
    ),
    values = c(
      "#3399FF",
      "#003366",
      "#3399CC",
      "#33CCFF"
    )
  ) +
  xlab("Frequência Alélica p") +
  ylab("Frequência Genotípica") +
  labs(color = "", shape = "")+
  theme_bw() +
  ggrepel::geom_text_repel(data = data_mean_gather, aes(x = Prop_R_Al, y = Values,
                                                         label = Generation))+
  facet_wrap(~Tratamento )
Fig3


