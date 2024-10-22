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

org.data = read.csv('../cleaned.infertile.data.csv')
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
  # infertility.data = infertility.data[infertility.data$Fertile == 1 |grepl('RPL', infertility.data$Infertility.diagnosis),]
  # infertility.data = infertility.data[infertility.data$Fertile == 1 |grepl('DOR', infertility.data$Infertility.diagnosis),]
  # infertility.data = infertility.data[infertility.data$Fertile == 1 |grepl('PCOS', infertility.data$Infertility.diagnosis),]
  # infertility.data = infertility.data[infertility.data$Fertile == 1 |grepl('endometriosis', infertility.data$Infertility.diagnosis),]
  # infertility.data = infertility.data[infertility.data$Fertile == 1 |grepl('explained', infertility.data$Infertility.diagnosis),]
  
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
  fertile.data = bind_cols(infertility.data[, c('Fertile','AMH', 'AFC',  'CD.of.blood.draw',
                                                # 'CORTISOL', 'ESTRADIOL',	'Progestrone',	'TESTOSTERONE',	'T3',	'T4',
                                                'Age', 'BMI')], cytokines.data)
  fertile.data$Fertile = factor(1 - fertile.data$Fertile)
  fertile.data$Age = log(fertile.data$Age, 2)
  fertile.data$BMI = log(fertile.data$BMI, 2)
  fertile.data$AMH = log(fertile.data$AMH, 2)
  fertile.data$AFC = log(fertile.data$AFC, 2)
  # fertile.data$Age = rnorm(nrow(fertile.data), mean = log(24,2), sd = 3)
  # fertile.data$BMI = rnorm(nrow(fertile.data), mean = log(24,2), sd = 3)
  
  
  ## remove columsn with single unique value ##
  # fertile.data = Filter(function(x)(length(unique(x))>1), fertile.data)
  # str(fertile.data)
  
  lm.model = lm(Fertile~AMH + BMI + Age, data = fertile.data)
  lm.model = lm(AMH~Fertile + BMI + Age, data = fertile.data)
  summary(lm.model)
  glm.model = glm(Fertile~AMH + BMI + Age, data = fertile.data, family = binomial(link="logit"))
  #glm.model = glm(as.formula(paste0('Fertile~ ', paste(colnames(fertile.data)[4:ncol(fertile.data)], collapse = '+'))),
  #                data = fertile.data, family = binomial(link="logit"))
  summary(glm.model)
  
  for (col in colnames(fertile.data)){
    fertile.data[is.na(fertile.data[,col]),col] = mean(fertile.data[,col], na.rm = T)
  }
  
  
  # cv.amh <- cv.glmnet(as.matrix(fertile.data[,c('AMH', 'AFC', 'Age', 'BMI')]),
  cv.amh <- cv.glmnet(as.matrix(fertile.data[,c('AMH', 'AFC', 'Age', 'BMI')]),
                      fertile.data$Fertile,
                      family = "binomial", nfold = nrow(fertile.data), type.measure = "auc", paralle = TRUE)
  
  # plot(cv.amh)
  
  return_features = function(coeff){
    top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
    top_features = top_features[order(top_features$coefficient),]
    # print(top_features)
    return(top_features)
  }
  
  tvec.amh = return_features(coef(cv.amh, s = 'lambda.min'))$coefficient
  names(tvec.amh) = return_features(coef(cv.amh, s = 'lambda.min'))$name
  (tvec.amh)
  # 
  # cv <- cv.glmnet(as.matrix(fertile.data[,c(5:ncol(fertile.data))]),
  # # cv <- cv.glmnet(as.matrix(fertile.data[,c(3:ncol(fertile.data))]),
  # # cv <- cv.glmnet(as.matrix(fertile.data[,c(4:10)]),
  #                 fertile.data$Fertile,
  #                 family = "binomial", nfold = nrow(fertile.data), type.measure = "auc", paralle = TRUE)
  
  # plot(cv)
  
  
  # tvec = return_features(coef(cv, s = 'lambda.min'))$coefficient
  # names(tvec) = return_features(coef(cv, s = 'lambda.min'))$name
  # (tvec)
  
  res.coef = c()
  max.auc = c()
  res.coef.values = c()
  for (all.rep in 1:100){
    for (j in 1:10){
      #print(paste('I: ', all.rep, 'J: ', j))
      print(paste('I: ', all.rep))
      md3cv <- cv.glmnet(as.matrix(scale(fertile.data[,c(4:ncol(fertile.data))])),
                         fertile.data$Fertile,
                         family = "binomial", nfold = 3, type.measure = "auc", paralle = TRUE, alpha = 0.8)
      # plot(md3cv)
      
      return_features = function(coeff){
        top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
        top_features = top_features[order(top_features$coefficient),]
        # print(top_features)
        return(top_features)
      }
      
      tvec = return_features(coef(md3cv, s = 'lambda.min'))$coefficient
      names(tvec) = return_features(coef(md3cv, s = 'lambda.min'))$name
      
      res.coef.values = c(res.coef.values, tvec)
      
      res.coef = c(res.coef, return_features(coef(md3cv, s = 'lambda.min'))$name)
      max.auc = c(max.auc, max(md3cv$cvm))
      print(return_features(coef(md3cv, s = 'lambda.min'))$name)
      # print(res.coef)
      # print(max.auc)
      # res.coef$name
    }
    print(mean(max.auc))
    print(sd(max.auc))
    #jpeg(paste('pred_within_subtypes_res/Moderate/lasso/Moderate_binary_glmnet_example.jpeg', sep = ''),
    #     units="in", width=10, height=10, res=500)
    #plot(md3cv)
    #dev.off()
  }
  
  df.res.coef.values = data.frame(res.coef.values)
  df.res.coef.values$names = names(res.coef.values)
  
}

library(ggplot2)
library(dplyr)
library(tibble)

## count from the 100 reps
binary.coef.values = df.res.coef.values
genes.count = table(binary.coef.values$names)
df.genes.count = data.frame(genes.count)
df.genes.count = df.genes.count %>% column_to_rownames(var = 'Var1')
hist(genes.count[order(genes.count, decreasing = T)])
(genes.count[order(genes.count, decreasing = T)])

coef.means = binary.coef.values %>% group_by(names) %>%
  summarize(coef_means = mean(res.coef.values, na.rm = TRUE))

## weighted coef base on number of appearances in the 100 reps
coef.means = coef.means %>% column_to_rownames(var = "names") 
coef.means
# coef.means * 
weighted.coef = coef.means * df.genes.count[rownames(coef.means),]
weighted.coef = weighted.coef[-1,, drop = FALSE]
head(weighted.coef)
weighted.coef = data.frame(scale(weighted.coef))
weighted.coef = weighted.coef[order(weighted.coef, decreasing = T),,drop = FALSE ]
weighted.coef

#write.csv((genes.count[order(genes.count, decreasing = T)])[-1], 'res/all_proteins/binary_ordered_genes_count.csv')
#write.csv(weighted.coef, 'res/all_proteins/binary_weighted_coef.csv')
#write.csv(coef.means, 'res/all_proteins/binary_mean_coef.csv')

graph.data = coef.means
graph.data$freq = genes.count[rownames(graph.data)]
graph.data = graph.data[-1,]

graph.data$names = sapply(strsplit(rownames(graph.data),'\\.'), function(x) x[1])

library(ggrepel)
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

### cytokines scores
# remove intercept
score.coef = coef.means[-1,,drop = F]
score.coef
# remove intercept, age and bmi
score.coef = coef.means[-c(1,2,3),,drop = F]
scores = as.matrix(fertile.data[,rownames(score.coef)])%*% as.matrix(score.coef$coef_means)%>% as.vector()
# scores = fertile.data[,rownames(score.coef)]%>%rowSums()
scores

auc(roc(fertile.data$Fertile, scores))

plot.data = rownames_to_column(coef.means)[-1,]
plot.data$rowname = factor(plot.data$rowname, levels = plot.data[order(plot.data$coef_means), 'rowname'])
ggplot(plot.data, aes(x = rowname, y = coef_means))+
  geom_bar(stat = 'Identity') +
  theme_bw()

write.csv(graph.data, 'cytokines_coef_onlyR.csv')


## sub cytokines scores
graph.data
sub.coef = graph.data[graph.data$freq > 30,]
score.coef = sub.coef[-c(1,2),,drop = F]
scores = as.matrix(fertile.data[,rownames(score.coef)])%*% as.matrix(score.coef$coef_means)%>% as.vector()
# scores = fertile.data[,rownames(score.coef)]%>%rowSums()
scores

auc(roc(fertile.data$Fertile, scores))

sub.coef$names = factor(sub.coef$names, levels = sub.coef[order(sub.coef$coef_means), 'names'])
sub.coef$names = factor(sub.coef$names, levels = sub.coef[order(abs(sub.coef$coef_means)), 'names'])
ggplot(sub.coef, aes(x = names, y = coef_means))+
  geom_bar(stat = 'Identity') +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

###AMH scores
# amh.score.coef = data.frame(tvec.amh)[-1, , drop = F]
#amh.score.coef = data.frame(tvec.amh)[-c(1,3,4), , drop = F]
tvec.amh
# remove intercept, age and bmi
amh.score.coef = data.frame(tvec.amh)[-c(1,4,5), , drop = F]
# amh.score.coef = data.frame(tvec.amh)[-c(1,3,4,5), , drop = F]
amh.score.coef
amh.scores = as.matrix(fertile.data[,rownames(amh.score.coef)])%*% as.matrix(amh.score.coef$tvec.amh)%>% as.vector()
amh.scores
auc(roc(fertile.data$Fertile, amh.scores))

###AMH AFC org scores
org.amh.afc.scores = as.matrix(fertile.data$AMH + fertile.data$AFC)%>% as.vector()
auc(roc(fertile.data$Fertile, org.amh.afc.scores))

##amh
org.amh.scores = as.matrix(fertile.data$AMH)%>% as.vector()
auc(roc(fertile.data$Fertile, org.amh.scores))

## afc
org.afc.scores = as.matrix(fertile.data$AFC)%>% as.vector()
auc(roc(fertile.data$Fertile, org.afc.scores))

## age
org.age.scores = as.matrix(fertile.data$Age)%>% as.vector()
auc(roc(fertile.data$Fertile, org.age.scores))

#weighted scores
weighted.coef
weighted.score.coef = weighted.coef[-c(1,2),,drop = F]
weighted.score.coef
weighted.scores = as.matrix(fertile.data[,rownames(weighted.score.coef)])%*% as.matrix(weighted.score.coef$coef_means)%>% as.vector()
auc(roc(fertile.data$Fertile, weighted.scores))

