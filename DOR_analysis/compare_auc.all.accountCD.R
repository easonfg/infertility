rm(list=ls())

library(pROC)
library(glmnet)
library(ggpubr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(aod)

org.data = read.csv('../cleaned.infertile.data.csv')
org.data$Infertility.diagnosis = toupper(org.data$Infertility.diagnosis)
# org.data = org.data[org.data$CD.of.blood.draw == 'R',]

org.data = org.data[-grep('DOR, ENDO', org.data$Infertility.diagnosis), ]

LOD = org.data[org.data$X == 'LOD', , drop = F]

all.res = data.frame()

### filter out missing values
infertility.data = org.data[org.data$Missing != 1,]

diag = 'DOR'
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


### turn all values lower than LOD to LOD
for (col.i in colnames(cytokines.data)){
  LOD_x = as.numeric(LOD[col.i])
  cytokines.data[which(cytokines.data[,col.i] < LOD_x), col.i] = LOD_x
}

cytokines.data = data.frame(cytokines.data)

# fertile.data = bind_cols(infertility.data[, c('Fertile','AMH', 'AFC', 'CD.of.blood.draw', 'Age', 'BMI')], cytokines.data)
fertile.data = bind_cols(infertility.data[, c('Fertile','AMH', 'AFC',  'CD.of.blood.draw',
                                              'CORTISOL', 'ESTRADIOL',	'Progestrone',	'TESTOSTERONE',	'T3',	'T4',
                                              'Age', 'BMI')], cytokines.data)
fertile.data$Fertile = factor(1 - fertile.data$Fertile)
fertile.data$Age = log(fertile.data$Age, 2)
fertile.data$BMI = log(fertile.data$BMI, 2)
fertile.data$AMH = log(fertile.data$AMH, 2)
fertile.data$AFC = log(fertile.data$AFC, 2)

for (col in colnames(fertile.data)){
  fertile.data[is.na(fertile.data[,col]),col] = mean(fertile.data[,col], na.rm = T)
}

fertile.data = Filter(function(x)(length(unique(x))>1), fertile.data)

### load coefficients
coef.means = read.csv('cytokines_coef_accountCD.csv')
amh.afc.coef.means = read.csv('AMH_AFC.csv')

coef.means = coef.means[-grep('CD_of_', coef.means$names),]

coef.means
get_auc = function(score.coef){
  cyto.scores = as.matrix(fertile.data[,score.coef$names])%*% as.matrix(score.coef$coef_means)%>% as.vector()
  # scores = fertile.data[,rownames(score.coef)]%>%rowSums()
  print(score.coef)
  
  plot(roc(fertile.data$Fertile, cyto.scores))
  
  print(auc(roc(fertile.data$Fertile, cyto.scores)))
}


coef.means
# sub.coef = coef.means[coef.means$freq > 3,]
sub.coef = coef.means[coef.means$freq > 0,]

score.coef = sub.coef[-c(1,2),,drop = F] # sans age bmi
get_auc(score.coef)
score.coef = sub.coef[-c(2),,drop = F] # keep age
# score.coef
get_auc(score.coef)
score.coef = sub.coef[-c(1),,drop = F] # keep bmi
score.coef
get_auc(score.coef)
score.coef = sub.coef[,,drop = F] # all
score.coef
get_auc(score.coef)
# score.coef = sub.coef
score.coef

#cyto.scores = as.matrix(fertile.data[,score.coef$names])%*% as.matrix(score.coef$coef_means)%>% as.vector()
## scores = fertile.data[,rownames(score.coef)]%>%rowSums()
#cyto.scores
#
#auc(roc(fertile.data$Fertile, cyto.scores))

colnames(amh.afc.coef.means)[1] = 'names'
hormone.score.coef = amh.afc.coef.means[c(2,3,4,5),,drop = F] ## all
get_auc(hormone.score.coef)
hormone.score.coef = amh.afc.coef.means[c(2,4),,drop = F] #afc, amh
get_auc(hormone.score.coef)
hormone.score.coef = amh.afc.coef.means[c(2,3,4),,drop = F] #afc, amh, age
get_auc(hormone.score.coef)
hormone.score.coef = amh.afc.coef.means[c(2,3,5),,drop = F] #afc, amh, bmi
get_auc(hormone.score.coef)
hormone.score.coef = amh.afc.coef.means[c(4),,drop = F] # amh
get_auc(hormone.score.coef)
hormone.score.coef = amh.afc.coef.means[c(2),,drop = F] # afc
get_auc(hormone.score.coef)
hormone.score.coef = amh.afc.coef.means[c(3),,drop = F] # age
get_auc(hormone.score.coef)
hormone.score.coef = amh.afc.coef.means[c(3,5),,drop = F] # age, bmi
get_auc(hormone.score.coef)
hormone.score.coef = amh.afc.coef.means[c(5),,drop = F] # bmi
get_auc(hormone.score.coef)
# hormone.score.coef = coef.means[-c(1,5),,drop = F]
hormone.score.coef
scores = as.matrix(fertile.data[,hormone.score.coef$names])%*% as.matrix(hormone.score.coef$coef_means)%>% as.vector()
# scores = fertile.data[,rownames(score.coef)]%>%rowSums()
scores

auc(roc(fertile.data$Fertile, scores))

# roc.test(auc(roc(fertile.data$Fertile, scores)), auc(roc(fertile.data$Fertile, cyto.scores)))
# roc.test(auc(roc(fertile.data$Fertile, scores)), auc(roc(fertile.data$Fertile, cyto.scores)), method = 'bootstrap')


#roc_obj.test <- roc(fertile.data$Fertile, as.matrix(fertile.data[,sub.coef[,,drop = F]$names])%*% as.matrix(sub.coef[,,drop = F]$coef_means)%>% as.vector())
plot.test = plot(roc(fertile.data$Fertile, as.matrix(fertile.data[,sub.coef[,,drop = F]$names])%*% as.matrix(sub.coef[,,drop = F]$coef_means)%>% as.vector()),
                 print.auc = F, col = "red")
plot.test = plot(roc(fertile.data$Fertile, as.matrix(fertile.data[,amh.afc.coef.means[c(2,3,4,5),,drop = F]$names])%*% as.matrix(amh.afc.coef.means[c(2,3,4,5),,drop = F]$coef_means)%>% as.vector()),
                 print.auc = F, col = "green", lty = 2, print.auc.y = .4,  add = TRUE)
plot.test = plot(roc(fertile.data$Fertile, as.matrix(fertile.data[,amh.afc.coef.means[c(3),,drop = F]$names])%*% as.matrix(amh.afc.coef.means[c(3),,drop = F]$coef_means)%>% as.vector()),
                 print.auc = F, col = "blue",print.auc.y = .4,  add = TRUE)
legend("bottomright", (c('Cytokines+Age+BMI','AMH+AFC+Age+BMI',  'Age')), lty=1, lwd = 3,
       bty="n", col = c('red', 'green', 'blue'))
#dev.off()

graph.data = coef.means
print(
  ggplot(graph.data, aes(coef_means, freq, label = names)) +
    # ggplot(graph.data, aes(freq, coef_means)) +
    geom_point(aes(color = abs(coef_means))) +
    # geom_text(aes(label=names),size = 3, hjust=-0.5, vjust=0.5) +
    # geom_text(aes(label=ifelse(freq>75,as.character(names),'')),size = 3, hjust=-0.5, vjust=0.5) +
    geom_text_repel() +
    # nudge_y       = 32 - subset(graph.data, freq > 25)$freq,
    # size          = 4,
    # box.padding   = 1.5,
    # point.padding = 0.5,
    # force         = 100,
    # segment.size  = 0.2,
    # segment.color = "grey50",
    # direction     = "x") +
    xlab("Average Coefficients") + 
    ylab("Frequencies") +
    labs(color='Absolute Average\ 
         Coefficients') + ## legend title
    ggtitle('DOR Variables')+
    theme_bw()
)

plot.data = coef.means
plot.data
plot.data$names = factor(plot.data$names, levels = plot.data[order(plot.data$coef_means), 'names'])
plot.data$names = factor(plot.data$names, levels = plot.data[order(abs(plot.data$coef_means)), 'names'])
ggplot(plot.data, aes(x = names, y = coef_means))+
  geom_bar(stat = 'Identity') +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

