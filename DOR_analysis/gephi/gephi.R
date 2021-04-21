## see which cytokines can predict DOR

rm(list=ls())

library(pROC)
library(glmnet)
library(ggpubr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(aod)

org.data = read.csv('../../cleaned.infertile.data.csv')
org.data$Infertility.diagnosis = toupper(org.data$Infertility.diagnosis)

org.data = org.data[-grep('DOR, ENDO', org.data$Infertility.diagnosis), ]
org.data = org.data[org.data$CD.of.blood.draw == 'R', ]

LOD = org.data[org.data$X == 'LOD', , drop = F]

all.res = data.frame()
# for (diag in c('RPL', 'DOR', 'PCOS', 'ENDOMETRIOSIS', 'UNEXPLAINED')){
# for (diag in c('DOR', 'PCOS', 'ENDOMETRIOSIS', 'UNEXPLAINED')){
for (diag in c('DOR')){
  # for (diag in c('ENDOMETRIOSIS')){
  
  ### filter out missing values
  infertility.data = org.data[org.data$Missing != 1,]
  
  ### keep only certain diseases
  
  infertility.data = infertility.data[infertility.data$Fertile == 1 |grepl(diag, infertility.data$Infertility.diagnosis),]
  #infertility.data %>% View()
  infertility.data$Infertility.diagnosis
  
  ### retain only values from R collection date
  #infertility.data = infertility.data[grep('R$', infertility.data$CD.of.blood.draw),]
  # infertility.data = infertility.data[grep('R', infertility.data$CD.of.blood.draw),]
  
  ### remove NA row
  infertility.data = infertility.data[!is.na(infertility.data$X),]
  #infertility.data %>% View()
  #rownames(infertility.data))
  
  ### extract cytokines data
  cytokines.data = infertility.data[grep('IL8', colnames(infertility.data)):grep('CSF1', colnames(infertility.data))]
  
  ### get rid of cytokines with more than 50% missing values
  # n.missing = sapply(colnames(cytokines.data), function(x){
  #   LOD = cytokines.data[nrow(cytokines.data),x]
  #   sum(cytokines.data[,x] < LOD)
  # })
  # cytokines.data = cytokines.data[,n.missing/nrow(cytokines.data) < 0.5]
  
  ### turn all values lower than LOD to LOD
  for (col.i in colnames(cytokines.data)){
    LOD_x = as.numeric(LOD[col.i])
    cytokines.data[which(cytokines.data[,col.i] < LOD_x), col.i] = LOD_x
  }
  
  cytokines.data = data.frame(cytokines.data)
  
  # fertile.data = bind_cols(infertility.data[, c('Fertile','AMH', 'AFC', 'CD.of.blood.draw', 'Age', 'BMI')], cytokines.data)
  fertile.data = bind_cols(infertility.data[, c('Fertile','AMH', 'AFC',  'CD.of.blood.draw', 'Infertility.diagnosis',
                                                # 'CORTISOL', 'ESTRADIOL',	'Progestrone',	'TESTOSTERONE',	'T3',	'T4',
                                                'Age', 'BMI')], cytokines.data)
  ## switch fertile
  fertile.data$DOR = factor(1 - fertile.data$Fertile)
  
  ## log variables
  fertile.data$Age = log(fertile.data$Age, 2)
  fertile.data$BMI = log(fertile.data$BMI, 2)
  fertile.data$AMH = log(fertile.data$AMH, 2)
  fertile.data$AFC = log(fertile.data$AFC, 2)
  # fertile.data$Age = rnorm(nrow(fertile.data), mean = log(24,2), sd = 3)
  # fertile.data$BMI = rnorm(nrow(fertile.data), mean = log(24,2), sd = 3)
  
}

fertile.data$Fertile
fertile.data$Infertility.diagnosis
fertile.data$DOR

index.1st.cyto = grep('IL8', colnames(fertile.data))
index.last.cyto = grep('CSF1', colnames(fertile.data))


generate.for.gephi = function(dor, outname){
  dor.cor = cor(fertile.data[fertile.data$DOR == dor,index.1st.cyto:index.last.cyto])
  dor.cor[lower.tri(dor.cor)] = NA
  dor.cor
  dor.cor[,!rowSums(is.na(dor.cor)) == ncol(dor.cor)]
  
  dor.cor.dummy = dor.cor
  dor.cor.id = data.frame(Id = 1:length(dor.cor.dummy %>% colnames()), label = colnames(dor.cor.dummy))
  dor.cor.id %>% head()
  write.csv(dor.cor.id, paste0(outname, '.id.csv'), row.names = F)
  
  colnames(dor.cor.dummy) = 1:length(colnames(dor.cor.dummy))
  rownames(dor.cor.dummy) = 1:length(colnames(dor.cor.dummy))
  xy <- t(combn(colnames(dor.cor.dummy), 2))
  xy
  dor.cor.list = data.frame(xy, cor=dor.cor.dummy[xy])
  dor.cor.list$cor = abs(dor.cor.list$cor)
  dor.cor.list%>%head()
  # dor.cor.list = dor.cor.list[abs(dor.cor.list$cor) > 0.8,]
  colnames(dor.cor.list)[1:3] = c('Source', 'Target', 'Weight')
  dor.cor.list = dor.cor.list[!is.na(dor.cor.list$Source),]
  write.csv(dor.cor.list, paste0(outname, '.edges.csv'), row.names = F)
}

generate.for.gephi(0, 'CTRL2.weighted')
generate.for.gephi(1, 'DOR2.weighted')
