### paired POFI

rm(list=ls())

library(ggpubr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(aod)

org.data = read.csv('cleaned.infertile.data.csv')
org.data$Infertility.diagnosis = toupper(org.data$Infertility.diagnosis)
org.data

LOD = org.data[org.data$X == 'LOD', , drop = F]
LOD


org.data = org.data[grep('POFI', org.data$X),]
org.data$Sample.ID = sapply(strsplit(org.data$Sample.ID, ' '), function(x)x[1])
org.data$Sample.ID



sub.org.data = org.data[-nrow(org.data),]

## get rid of NA AMH
sub.org.data = sub.org.data[!is.na(sub.org.data$AMH),]

### filter out missing values
infertility.data = sub.org.data[sub.org.data$Missing != 1,]

### remove NA row
infertility.data = infertility.data[!is.na(infertility.data$X),]
infertility.data
#infertility.data %>% View()
#rownames(infertility.data))

### extract cytokines data
cytokines.data = infertility.data[grep('IL8', colnames(infertility.data)):grep('CSF1', colnames(infertility.data))]

### turn all values lower than LOD to LOD
cytokines.data[nrow(cytokines.data),]
for (col.i in colnames(cytokines.data)){
  LOD_x = as.numeric(LOD[col.i])
  cytokines.data[which(cytokines.data[,col.i] < LOD_x), col.i] = LOD_x
}

cytokines.data = data.frame(cytokines.data)
# cytokines.data %>% View()

infertility.data[infertility.data$Fertile == 1, 'Infertility.diagnosis'] = 'FERTILE'
infertility.data$Infertility.diagnosis = gsub('MALE FACTOR, ', '', infertility.data$Infertility.diagnosis)
infertility.data$Infertility.diagnosis = gsub('MALE, ', '', infertility.data$Infertility.diagnosis)
infertility.data$Infertility.diagnosis = gsub(', MALE FACTOR', '', infertility.data$Infertility.diagnosis)
infertility.data[grep('UNEXPLAINED', infertility.data$Infertility.diagnosis), 'Infertility.diagnosis'] = 'UNEXPLAINED'
infertility.data[grep('PCOS', infertility.data$Infertility.diagnosis), 'Infertility.diagnosis'] = 'PCOS'
infertility.data[grep('DOR', infertility.data$Infertility.diagnosis), 'Infertility.diagnosis'] = 'DOR'
infertility.data[grep('ENDOM', infertility.data$Infertility.diagnosis), 'Infertility.diagnosis'] = 'ENDOMETRIOSIS'
infertility.data[grep('UTERINE', infertility.data$Infertility.diagnosis), 'Infertility.diagnosis'] = 'UTERINE'
infertility.data$Infertility.diagnosis
colnames(infertility.data)

fertile.data = bind_cols(infertility.data[, c('Cycle.Date', 'Cycle.Type', 'Sample.ID',
                                              'Age', 'BMI', 'Fertile', 'Infertility.diagnosis')], cytokines.data)
fertile.data
fertile.data$Age = log(fertile.data$Age, 2)
fertile.data$BMI = log(fertile.data$BMI, 2)

library(lme4)
library(lmerTest)


# library(lmerTest)
### natural
sub.org.data = fertile.data[grep('IVF|natural', fertile.data$Cycle.Type),]
sub.org.data$Cycle.Type

res.natural = sapply(colnames(sub.org.data)[grep('IL8', colnames(sub.org.data)):ncol(sub.org.data)], function(x){
  
  # model = lmer(as.formula(paste(x, '~ Cycle.Type + Age+ (1|Sample.ID)')), data=sub.org.data,REML=TRUE)
  # model = lmer(as.formula(paste(x, '~ Cycle.Type + Age +BMI + Fertile+  (1|Sample.ID)')), data=sub.org.data,REML=TRUE)
  model = lmer(as.formula(paste(x, '~ Cycle.Type + Age +BMI + Infertility.diagnosis+ (1|Sample.ID)')), data=sub.org.data,REML=TRUE)
  # model = lmer(as.formula(paste(x, '~ Cycle.Type + Age +BMI + Fertile+ Infertility.diagnosis+ (1|Sample.ID)')), data=sub.org.data,REML=TRUE)
  c(summary(model)$coefficients[2,c(1, 5)])
})

res.natural = t(res.natural) %>% data.frame()
res.natural[,'padj'] = p.adjust(res.natural[,2], method = 'fdr')
res.natural[res.natural$padj<0.1,]
res.natural[res.natural$Pr...t..<0.05,]

### medicated
sub.org.data = fertile.data[grep('IVF|medicated', fertile.data$Cycle.Type),]
sub.org.data$Cycle.Type

res.medicated = sapply(colnames(sub.org.data)[grep('IL8', colnames(sub.org.data)):ncol(sub.org.data)], function(x){
  print(x)
  if (length(unique(sub.org.data[,x])) > 1){
    
    model = lmer(as.formula(paste(x, '~ Cycle.Type + Age+ BMI + Infertility.diagnosis + (1|Sample.ID)')), data=sub.org.data,REML=TRUE)
    c(summary(model)$coefficients[2,c(1, 5)])
  }else{return()}
})

res.medicated = Filter(Negate(is.null), res.medicated)
res.medicated
res.medicated = do.call(rbind, res.medicated) %>% data.frame()


head(res.medicated)
res.medicated[,'padj'] = p.adjust(res.medicated[,2], method = 'fdr')
res.medicated[res.medicated$padj<0.1,]
res.medicated[res.medicated$Pr...t..<0.05,]


### medicated vs natural
sub.org.data = fertile.data[grep('natural|medicated', fertile.data$Cycle.Type),]
sub.org.data$Cycle.Type

res.medicated.natural = sapply(colnames(sub.org.data)[grep('IL8', colnames(sub.org.data)):ncol(sub.org.data)], function(x){
  print(x)
  if (length(unique(sub.org.data[,x])) > 1){
    
    model = lmer(as.formula(paste(x, '~ Cycle.Type + Age+BMI + Infertility.diagnosis +  (1|Sample.ID)')), data=sub.org.data,REML=TRUE)
    c(summary(model)$coefficients[2,c(1, 5)])
  }else{return()}
})

res.medicated.natural = Filter(Negate(is.null), res.medicated.natural)
res.medicated.natural
res.medicated.natural = do.call(rbind, res.medicated.natural) %>% data.frame()


head(res.medicated.natural)
res.medicated.natural[,'padj'] = p.adjust(res.medicated.natural[,2], method = 'fdr')
res.medicated.natural[res.medicated.natural$padj<0.1,]
res.medicated.natural[res.medicated.natural$Pr...t..<0.05,]


res.natural$diag ='natural_IVF'
res.medicated$diag ='medicated_IVF'
res.medicated.natural$diag ='natural_medicated'
res = bind_rows(res.natural, res.medicated, res.medicated.natural)
res%>%head()
colnames(res)[2] = 'P'

## change rownames to column
res = res %>% rownames_to_column('Cytokines')
res$Cytokines = sapply(strsplit(res$Cytokines, '\\.\\.\\.'), function(x) x[1])

### dot plot of p values
dot.plot.data = res
dot.plot.data%>%head()
# dot.plot.data$padj = dot.plot.data$padj + 0.0000001

gg.res = ggplot(subset(dot.plot.data, P < 0.05), aes(diag, Cytokines)) +
# gg.res = ggplot(subset(dot.plot.data, padj < 0.1), aes(diag, Cytokines)) +
  geom_point(aes(size = -log10(padj), fill = Estimate), 
             colour = 'black', shape = 21) + 
  geom_text(aes(label=round(Estimate,2)),size = 3, hjust=-0.5, vjust=0.5) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  labs(size = '-log10(Raw P values)', fill = 'coefficients') +
  theme_bw()
print(gg.res)

res%>%head()
sig.all.res = res[res$P < 0.05,]
sig.all.res%>%head()
long.fertile.data = gather(fertile.data, cytokines, measurements, IL8:CSF1, factor_key=TRUE)
cytokines.int =  unique(sig.all.res$Cytokines)
cytokines.int
# long.fertile.data = long.fertile.data[long.fertile.data$cytokines %in% cytokines.int,]
long.fertile.data

# long.fertile.data$cytokines = factor(long.fertile.data$cytokines, levels = dot.plot.data$cytokines[order(dot.plot.data$Estimate)])




for (hormone in c(sig.all.res$diag%>%unique())){
  cyto.int = sig.all.res[sig.all.res$diag %in% hormone, 'Cytokines']
  # cyto.int = long.fertile.data$cytokines %>% unique()
  cyto.chunks.cd = split(cyto.int, ceiling(seq_along(cyto.int)/10))
  print(hormone)
  # tmp.long.fertile.data = long.fertile.data[long.fertile.data$Cycle.Type %in% strsplit(hormone, '_')[[1]],]
  pdf(paste0('ivf_med_natural_comp/', hormone, '_regression', '.pdf'), height = 10, width = 20, onefile = TRUE)
  for (i in cyto.chunks.cd){
    temp.fertile.data.diag = long.fertile.data[long.fertile.data$cytokines %in% i, ]
    print(i)
    print(
      ggboxplot(temp.fertile.data.diag, x = "cytokines", y = 'measurements',
                color = "Cycle.Type", palette = "jco", add = "jitter") + 
        # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
        #theme(legend.position="none")+
        facet_wrap(vars(cytokines), scales = 'free', nrow = 2)
    )
  }
  dev.off()
}

#anova(model)
#
#rand(model)
#
################
#
#library(car)
#
#scatterplot(IL8 ~ BMI | Cycle.Type, data = sub.org.data, smooth=F, reg.line=F)
## scatterplot(BDNF ~ BMI | Cycle.Type, data = sub.org.data, smooth=F, reg.line=F)
#
################
#
#library(multcompView)
#
#library(lsmeans)
#
#leastsquare = lsmeans(model,
#                      pairwise ~ Cycle.Type,
#                      adjust="tukey")
#
#leastsquare
#
#CLD = cld(leastsquare,
#          alpha=0.05,
#          Letters=letters,      ### Use lower-case letters for .group
#          adjust="tukey")       ### Tukey-adjusted comparisons
#
#CLD
#
################################
#
#library(ggplot2)
#
#qplot(x    = Time ,
#      y    = lsmean,
#      data = CLD) +
#  
#  geom_errorbar(aes(
#    ymin  = lower.CL,
#    ymax  = upper.CL,
#    width = 0.15))
#
################################
#
#
#T1  = sub.org.data$IL8[sub.org.data$Cycle.Type=="IVF"]
#T2  = sub.org.data$IL8[sub.org.data$Cycle.Type=="natural"]
#
#
#Difference = T2 - T1
#
#X = Data$Individual[Data$Time=="T1"]
#
#
#barplot(Difference,
#        names.arg = X, 
#        col="dark gray",
#        xlab="Individual",
#        ylab="Difference (T2 ñ T1)")
#
##################################
#
#print(data.frame(Time = c("T1", "T2"), Mean = c(mean(T1), mean(T2))))
#
###################################################
#
#
#hist(residuals(model), col="darkgray")
#
#plot(predict(model), residuals(model))
#