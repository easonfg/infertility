rm(list=ls())

library(ggpubr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(aod)

org.data = read.csv('../cleaned.infertile.data.csv')
org.data$Infertility.diagnosis = toupper(org.data$Infertility.diagnosis)
org.data = org.data[-grep('DOR, ENDO', org.data$Infertility.diagnosis), ]

LOD = org.data[org.data$X == 'LOD', , drop = F]


all.res = data.frame()
# for (diag in c('RPL', 'DOR', 'PCOS', 'ENDOMETRIOSIS', 'UNEXPLAINED')){
for (diag in c( 'DOR' )){
  
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
  # infertility.data = infertility.data[grep('R$', infertility.data$CD.of.blood.draw),]
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
  cytokines.data[nrow(cytokines.data),]
  for (col.i in colnames(cytokines.data)){
    LOD_x = as.numeric(LOD[col.i])
    cytokines.data[which(cytokines.data[,col.i] < LOD_x), col.i] = LOD_x
  }
  
  cytokines.data = data.frame(cytokines.data)
  
  fertile.data = bind_cols(infertility.data[, c('Fertile', 'Age',
                                                'CORTISOL', 'ESTRADIOL',	'Progestrone',	'TESTOSTERONE',	'T3',	'T4',
                                                'Cycle.Date', 'BMI')], cytokines.data)
  fertile.data$Fertile = factor(1-fertile.data$Fertile)
  fertile.data$Age = log(fertile.data$Age, 2)
  fertile.data$BMI = log(fertile.data$BMI, 2)
  
  
  ## remove columsn with single unique value ##
  fertile.data = Filter(function(x)(length(unique(x))>1), fertile.data)
  
  res = sapply(c('CORTISOL', 'ESTRADIOL',	'Progestrone',	'TESTOSTERONE',	'T3',	'T4'), function(x){
    lm.model = lm(as.formula(paste(x, '~Fertile + Age + BMI + Cycle.Date + BMI')), data = fertile.data)
    # lm.model = lm(as.formula(paste(x, '~Fertile + Age + Cycle.Date')), data = fertile.data)
    # lm.model = lm(as.formula(paste(x, '~Fertile + Age')), data = fertile.data)
    temp.wt.res = wald.test(b=coef(lm.model), Sigma=vcov(lm.model), Terms=c(2))
    temp.wt.res$result$chi2[3]
    summary(lm.model)
    # c(summary(lm.model)$coefficients[2,c(1,4)])
    # c(summary(lm.model)$coefficients[2,c(1,4)], temp.wt.res$result$chi2[3])
    c(summary(lm.model)$coefficients[2,c(1)], temp.wt.res$result$chi2[3])
  })
  
  res = t(res) %>% data.frame()
  res$diag = diag
  #res[,'padj'] = p.adjust(res[,2], method = 'fdr')
  
  #print(res[res[,2] < 0.05,])
  
  all.res = bind_rows(all.res, res)
  
}

all.res
all.res[,'padj'] = p.adjust(all.res[,2], method = 'fdr')
print(all.res[all.res$padj < 0.1,])
# all.res[grep('CXCL9', rownames(all.res)),]

### dot plot of p values
dot.plot.data = all.res
dot.plot.data = dot.plot.data %>% rownames_to_column('Cytokines')
dot.plot.data$Cytokines = sapply(strsplit(dot.plot.data$Cytokines, '\\.\\.\\.'), function(x) x[1])
dot.plot.data

gg.res = ggplot(subset(dot.plot.data, P < 0.05), aes(diag, Cytokines)) + 
  geom_point(aes(size = -log10(padj), fill = V1), 
             colour = 'black', shape = 21) + 
  geom_text(aes(label=round(V1,2)),size = 3, hjust=-0.5, vjust=0.5) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  labs(size = '-log10(FDR)', fill = 'coefficients') +
  theme_bw()
print(gg.res)
### dot plot of p values


plot.data = org.data[grepl('RPL|DOR|PCOS|ENDOMETRIOSIS|UNEXPLAINED',org.data$Infertility.diagnosis) | org.data$Fertile == 1,]
# plot.data = plot.data[grep('R$', plot.data$CD.of.blood.draw),]
plot.data = plot.data[grep('R', plot.data$CD.of.blood.draw),]
# plot.data %>% View()

for (x in c('RPL','ENDOMETRIOSIS', 'DOR', 'PCOS',  'UNEXPLAINED')){
  plot.data[grep(x, plot.data$Infertility.diagnosis), 'Infertility.diagnosis'] = x
}
plot.data[plot.data$Fertile == 1, 'Infertility.diagnosis'] = 'FERTILE'
plot.data$Infertility.diagnosis = factor(plot.data$Infertility.diagnosis,
                                         levels = c('FERTILE', 'RPL', 'DOR', 'PCOS', 'ENDOMETRIOSIS', 'UNEXPLAINED'))

# ### plotting the significant R only cytokines (CXCL9)
# my_comparisons <- list( #c('FERTILE', 'DOR'),  c('FERTILE', 'PCOS'),
#   #c('FERTILE', 'RPL'),  c('FERTILE', 'UNEXPLAINED'),
#   c('FERTILE', 'ENDOMETRIOSIS'))
# ggboxplot(plot.data, x = "Infertility.diagnosis", y = "CXCL9",
#           color = "Infertility.diagnosis", palette = "jco", add = "jitter") + 
#   stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
#   theme(legend.position="none")
# # stat_compare_means(method = "t.test") 
# 
# plot.data %>% 
#   ggplot(aes(x=Infertility.diagnosis,y=CXCL9, label = CD.of.blood.draw, color=Infertility.diagnosis))+
#   geom_boxplot(width=.5)+
#   # jittered text with geom_text
#   geom_text(
#     position=position_jitter(width=0.15))+
#   theme_classic()+
#   theme(legend.position="none")
# ### plotting the significant R only cytokines (CXCL9)

### plotting the significant R+n cytokines 
my_comparisons <- list( c('FERTILE', 'DOR'),  c('FERTILE', 'PCOS'),
                        c('FERTILE', 'RPL'),  c('FERTILE', 'UNEXPLAINED'),
                        c('FERTILE', 'ENDOMETRIOSIS'))
my_comparisons <- list( c('FERTILE', 'DOR'))
ggboxplot(plot.data[grep('FERTILE|DOR', plot.data$Infertility.diagnosis),], x = "Infertility.diagnosis", y = "TNFSF14",
          color = "Infertility.diagnosis", palette = "jco", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(legend.position="none")

plot.data[grep('FERTILE|DOR', plot.data$Infertility.diagnosis),] %>% 
  ggplot(aes(x=Infertility.diagnosis,y=TNFSF14, label = CD.of.blood.draw, color=Infertility.diagnosis))+
  geom_boxplot(width=.5)+
  theme_classic()+
  # jittered text with geom_text
  geom_text(
    position=position_jitter(width=0.15))+
  theme(legend.position="none")

plot.data[grep('FERTILE|DOR', plot.data$Infertility.diagnosis),] %>% 
  ggplot(aes(x=Infertility.diagnosis,y=CD40, label = CD.of.blood.draw, color=Infertility.diagnosis))+
  geom_boxplot(width=.5)+
  theme_classic()+
  # jittered text with geom_text
  geom_text(
    position=position_jitter(width=0.15))+
  theme(legend.position="none")

plot.data[grep('FERTILE|DOR', plot.data$Infertility.diagnosis),] %>% 
  ggplot(aes(x=Infertility.diagnosis,y=ST1A1, label = CD.of.blood.draw, color=Infertility.diagnosis))+
  geom_boxplot(width=.5)+
  theme_classic()+
  # jittered text with geom_text
  geom_text(
    position=position_jitter(width=0.15))+
  theme(legend.position="none")

plot.data[grep('FERTILE|PCOS', plot.data$Infertility.diagnosis),] %>% 
  ggplot(aes(x=Infertility.diagnosis,y=FGF21, label = CD.of.blood.draw, color=Infertility.diagnosis))+
  geom_boxplot(width=.5)+
  theme_classic()+
  # jittered text with geom_text
  geom_text(
    position=position_jitter(width=0.15))+
  theme(legend.position="none")
### plotting the significant R+n cytokines 

sub.plot.data = plot.data[,c('CD.of.blood.draw', 'CXCL9', 'Infertility.diagnosis')]
sub.plot.data[sub.plot.data$Infertility.diagnosis == 'ENDOMETRIOSIS',]
#res = sapply(colnames(fertile.data)[3:ncol(fertile.data)], function(x){
#  lm.model = glm(as.formula(paste('Fertile ~', x,  '+Age')), family=binomial(link='logit'), data = fertile.data)
#  summary(lm.model)
#  c(summary(lm.model)$coefficients[2,c(1,4)])
#})
#
#res = t(res) %>% data.frame()
#res
#res[,'padj'] = p.adjust(res[,2], method = 'fdr')
#
#res[res$Pr...t.. < 0.05,]


