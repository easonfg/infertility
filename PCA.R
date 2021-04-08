rm(list=ls())

library(ggpubr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(aod)

org.data = read.csv('cleaned.infertile.data.csv')
org.data$Infertility.diagnosis = toupper(org.data$Infertility.diagnosis)
org.data = org.data[!is.na(org.data$X),]

LOD = org.data[org.data$X == 'LOD', , drop = F]

org.data = org.data[-nrow(org.data),]

org.data[org.data$Fertile == 1, 'Infertility.diagnosis'] = 'FERTILE'

###### INFERTILITY DIAGNOSIS #####
infertility.data = org.data[org.data$Missing != 1,]
infertility.data = infertility.data[!is.na(infertility.data$X),]
infertility.data = infertility.data[grep('R$', infertility.data$CD.of.blood.draw),]
infertility.data$Infertility.diagnosis
infertility.data = infertility.data[grep('FERTILE|PCOS|DOR|ENDOMETRIOSIS', infertility.data$Infertility.diagnosis), ]
infertility.data

for (x in c('ENDOMETRIOSIS', 'DOR', 'PCOS')){
  infertility.data[grep(x, infertility.data$Infertility.diagnosis), 'Infertility.diagnosis'] = x
}

cytokines.data = infertility.data[grep('IL8', colnames(infertility.data)):grep('CSF1', colnames(infertility.data))]


cytokines.data[nrow(cytokines.data),]
### turn all values lower than LOD to LOD
for (col.i in colnames(cytokines.data)){
  LOD_x = as.numeric(LOD[col.i])
  cytokines.data[cytokines.data[,col.i] < LOD_x, col.i] = LOD_x
}

pca.res = prcomp(cytokines.data)
library("factoextra")
fviz_eig(pca.res)
pca.plot = fviz_pca_ind(pca.res,
                        col.ind = as.factor(infertility.data$Infertility.diagnosis) # Color by the quality of representation
                        # geom="point", ## don't label
                        #repel = TRUE     # Avoid text overlapping
)
print(pca.plot)

print(
ggplot(data.frame(pca.res$x)) + geom_point(aes(x=PC1, y=PC2,
                                                  #color = as.factor(infertility.data$Cycle.Date),
                                                  color = as.factor(infertility.data$Infertility.diagnosis),
                                                  size = as.factor(infertility.data$Infertility.diagnosis == 'FERTILE'))) +
  scale_shape_manual(values=1:nlevels(as.factor(infertility.data$Infertility.diagnosis)))+
  stat_ellipse(aes(x=PC1, y=PC2,color=as.factor(infertility.data$Infertility.diagnosis), group=as.factor(infertility.data$Infertility.diagnosis)),type = "norm")+
  theme_classic()
)
###### INFERTILITY DIAGNOSIS #####


###### Cycle Date #####
cd.data = org.data[org.data$Missing != 1,]
cd.data = cd.data[!is.na(cd.data$X),]
date.counts = table(cd.data$Cycle.Date)
cd.data = cd.data[cd.data$Cycle.Date %in% names(date.counts[date.counts > 2]),]
#cd.data$Cycle.Date %>% table()
cd.data$Infertility.diagnosis = gsub(', MALE FACTOR', '', cd.data$Infertility.diagnosis)
cd.data$Infertility.diagnosis = gsub('MALE, ', '', cd.data$Infertility.diagnosis)
cd.data$Infertility.diagnosis = gsub('MALE FACTOR, ', '', cd.data$Infertility.diagnosis)
cd.data$Infertility.diagnosis[grep('SIS, TUB', cd.data$Infertility.diagnosis)] = 'ENDOMETRIOSIS'
cd.data$Infertility.diagnosis[grep('PCOS, AZOO', cd.data$Infertility.diagnosis)] = 'PCOS'
cd.data$Infertility.diagnosis[grep('RECURRENT P', cd.data$Infertility.diagnosis)] = 'RPL'
cd.data$Infertility.diagnosis[grep('^ENDO', cd.data$Infertility.diagnosis)] = 'ENDOMETRIOSIS'
cd.data$Infertility.diagnosis[grep('^UNEXPLAINED', cd.data$Infertility.diagnosis)] = 'UNEXPLAINED'

cd.cytokines.data = cd.data[grep('IL8', colnames(cd.data)):grep('CSF1', colnames(cd.data))]

cd.pca.res = prcomp(cd.cytokines.data)
library("factoextra")
fviz_eig(cd.pca.res)
cd.pca.plot = fviz_pca_ind(cd.pca.res,
                        col.ind = as.factor(cd.data$Cycle.Date), # Color by the quality of representation
                        geom="point" ## don't label
                        #repel = TRUE     # Avoid text overlapping
)+ xlim(-10,10) 
print(cd.pca.plot)

#cd.pca.plot = fviz_pca_ind(cd.pca.res,
#                           geom="point",  pointsize = 1, 
#                           habillage = cd.data$Cycle.Date, addEllipses=TRUE, ellipse.level=0.95
#) + xlim(-10,15) 

ggplot(data.frame(cd.pca.res$x)) + geom_point(aes(x=PC1, y=PC2,
                                                         color = as.factor(cd.data$Cycle.Date),
                                                         shape = as.factor(cd.data$Infertility.diagnosis),
                                                         size = as.factor(cd.data$Infertility.diagnosis == 'FERTILE'))) +
  scale_shape_manual(values=1:nlevels(as.factor(cd.data$Infertility.diagnosis)))+
  stat_ellipse(aes(x=PC1, y=PC2,color=as.factor(cd.data$Cycle.Date), group=as.factor(cd.data$Cycle.Date)),type = "norm")+
  theme_classic()
###### Cycle Date #####

###### Cycle Date sans R#####
cd.data.sans14 = org.data[org.data$Missing != 1,]
cd.data.sans14 = cd.data.sans14[!is.na(cd.data.sans14$X),]
date.counts = table(cd.data.sans14$Cycle.Date)
cd.data.sans14 = cd.data.sans14[cd.data.sans14$Cycle.Date %in% names(date.counts[date.counts > 2]),]
cd.data.sans14 = cd.data.sans14[cd.data.sans14$Cycle.Date != 14,]
cd.data.sans14$Cycle.Date %>% table()
cd.cytokines.data.sans14 = cd.data.sans14[grep('IL8', colnames(cd.data.sans14)):grep('CSF1', colnames(cd.data.sans14))]

cd.pca.res.sans14 = prcomp(cd.cytokines.data.sans14)
library("factoextra")
fviz_eig(cd.pca.res.sans14)
cd.pca.plot.sans14 = fviz_pca_ind(cd.pca.res.sans14,
                                  #habillage = cd.data.sans14$Infertility.diagnosis,
                        col.ind = as.factor(cd.data.sans14$Cycle.Date),
                        pointsize = 3, # Color by the quality of representation
                        geom="point" ## don't label
                        #repel = TRUE     # Avoid text overlapping
)+ xlim(-10,10) 
print(cd.pca.plot.sans14)

ggplot(data.frame(cd.pca.res.sans14$x)) + geom_point(aes(x=PC1, y=PC2,
                                                         color = as.factor(cd.data.sans14$Cycle.Date),
                                                         shape = as.factor(cd.data.sans14$Infertility.diagnosis),
                                                     size = as.factor(cd.data.sans14$Infertility.diagnosis == 'FERTILE'))) +
  scale_shape_manual(values=1:nlevels(as.factor(cd.data.sans14$Infertility.diagnosis)))+
  stat_ellipse(aes(x=PC1, y=PC2,color=as.factor(cd.data.sans14$Cycle.Date), group=as.factor(cd.data.sans14$Cycle.Date)),type = "norm")+
  theme_classic()
  


#cd.pca.plot = fviz_pca_ind(cd.pca.res,
#                           geom="point",  pointsize = 1, 
#                           habillage = cd.data$Cycle.Date, addEllipses=TRUE, ellipse.level=0.95
#) + xlim(-10,15) 
###### Cycle Date #####