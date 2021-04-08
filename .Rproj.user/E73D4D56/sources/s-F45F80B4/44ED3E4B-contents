## only using R

rm(list=ls())

library(ggpubr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(aod)

org.data = read.csv('cleaned.infertile.data.csv')
org.data$Infertility.diagnosis = toupper(org.data$Infertility.diagnosis)

org.data = org.data[-grep('DOR, ENDO', org.data$Infertility.diagnosis), ]


LOD = org.data[org.data$X == 'LOD', , drop = F]

all.res = data.frame()
# for (diag in c('RPL', 'DOR', 'PCOS', 'ENDOMETRIOSIS', 'UNEXPLAINED')){
# for (diag in c('DOR', 'PCOS', 'ENDOMETRIOSIS', 'UNEXPLAINED')){
for (diag in c('DOR', 'PCOS', 'ENDOMETRIOSIS')){
  
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
  infertility.data = infertility.data[grep('R$', infertility.data$CD.of.blood.draw),]
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
    cytokines.data[cytokines.data[,col.i] < LOD_x, col.i] = LOD_x
  }
  
  cytokines.data = data.frame(cytokines.data)
  
  fertile.data = bind_cols(infertility.data[, c('Fertile', 'Age', 'BMI')], cytokines.data)
  fertile.data$Fertile = factor(1 - fertile.data$Fertile)
  fertile.data$Age = log(fertile.data$Age, 2)
  fertile.data$BMI = log(fertile.data$BMI, 2)
  # fertile.data$Age = rnorm(nrow(fertile.data), mean = log(24,2), sd = 3)
  # fertile.data$BMI = rnorm(nrow(fertile.data), mean = log(24,2), sd = 3)
  
  
  ## remove columsn with single unique value ##
  fertile.data = Filter(function(x)(length(unique(x))>1), fertile.data)
  
  res = sapply(colnames(fertile.data)[4:ncol(fertile.data)], function(x){
    # print(x)
    # lm.model = lm(as.formula(paste(x, '~Fertile')), data = fertile.data)
    # lm.model = lm(as.formula(paste(x, '~Fertile + BMI')), data = fertile.data)
    # lm.model = lm(as.formula(paste(x, '~Fertile + Age')), data = fertile.data)
    lm.model = lm(as.formula(paste(x, '~Fertile + BMI + Age')), data = fertile.data)
    # print(summary(lm.model))
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
print(all.res[all.res$P < 0.05,])

sub_padj = c()
for (diag in all.res$diag %>% unique()){
  sub_padj = c(sub_padj, p.adjust(all.res[all.res$diag == diag, 'P'], method = 'fdr'))
}
all.res$sub.padj = sub_padj
all.res[all.res$sub.padj < 0.05,]


## change rownames to column
all.res = all.res %>% rownames_to_column('Cytokines')
all.res$Cytokines = sapply(strsplit(all.res$Cytokines, '\\.\\.\\.'), function(x) x[1])

### dot plot of p values
dot.plot.data = all.res
dot.plot.data$padj = dot.plot.data$padj + 0.0000001

# gg.res = ggplot(subset(dot.plot.data, P < 0.05), aes(diag, Cytokines)) + 
gg.res = ggplot(subset(dot.plot.data, padj < 0.1), aes(diag, Cytokines)) + 
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

### plotting the significant R only cytokines (CXCL9)
# plot.data = org.data[grepl('RPL|DOR|PCOS|ENDOMETRIOSIS|UNEXPLAINED',org.data$Infertility.diagnosis) | org.data$Fertile == 1,]
plot.data = org.data[grepl('DOR|PCOS|ENDOMETRIOSIS',org.data$Infertility.diagnosis) | org.data$Fertile == 1,]
plot.data = plot.data[grep('R$', plot.data$CD.of.blood.draw),]
# plot.data = plot.data[grep('R', plot.data$CD.of.blood.draw),]
# plot.data %>% View()

for (x in c('RPL','ENDOMETRIOSIS', 'DOR', 'PCOS',  'UNEXPLAINED')){
  plot.data[grep(x, plot.data$Infertility.diagnosis), 'Infertility.diagnosis'] = x
}
plot.data[plot.data$Fertile == 1, 'Infertility.diagnosis'] = 'FERTILE'
plot.data$Infertility.diagnosis = factor(plot.data$Infertility.diagnosis,
                                         levels = c('FERTILE', 'RPL', 'DOR', 'PCOS', 'ENDOMETRIOSIS', 'UNEXPLAINED'))

### plotting the significant R only cytokines (CXCL9)
my_comparisons <- list( #c('FERTILE', 'DOR'),  c('FERTILE', 'PCOS'),
                        #c('FERTILE', 'RPL'),  c('FERTILE', 'UNEXPLAINED'),
                        # c('FERTILE', 'ENDOMETRIOSIS'))
                        c('FERTILE', 'DOR'))
# ggboxplot(plot.data, x = "Infertility.diagnosis", y = "CXCL9",
ggboxplot(plot.data, x = "Infertility.diagnosis", y = "CASP8",
# ggboxplot(plot.data, x = "Infertility.diagnosis", y = "IL22RA1",
# ggboxplot(plot.data, x = "Infertility.diagnosis", y = "IL1B",
          color = "Infertility.diagnosis", palette = "jco", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
   theme(legend.position="none")
  # stat_compare_means(method = "t.test") 

plot.data %>% 
  ggplot(aes(x=Infertility.diagnosis,y=CXCL9, label = CD.of.blood.draw, color=Infertility.diagnosis))+
  geom_boxplot(width=.5)+
  # jittered text with geom_text
  geom_text(
            position=position_jitter(width=0.15))+
  theme_classic()+
  theme(legend.position="none")
### plotting the significant R only cytokines (CXCL9)


###plotting by diagnosis###
infertility.data = org.data[org.data$Missing != 1,]

### retain only values from R collection date
infertility.data = infertility.data[grep('R$', infertility.data$CD.of.blood.draw),]
# infertility.data = infertility.data[grep('R', infertility.data$CD.of.blood.draw),]

### remove NA row
infertility.data = infertility.data[!is.na(infertility.data$X),]

### extract cytokines data
cytokines.data = infertility.data[grep('IL8', colnames(infertility.data)):grep('CSF1', colnames(infertility.data))]
cytokines.data

### turn all values lower than LOD to LOD
for (col.i in colnames(cytokines.data)){
  LOD_x = as.numeric(LOD[col.i])
  cytokines.data[cytokines.data[,col.i] < LOD_x, col.i] = LOD_x
}

fertile.data = bind_cols(infertility.data[, c('Fertile', 'Infertility.diagnosis')], cytokines.data)

long.fertile.data = gather(fertile.data, cytokines, measurements, IL8:CSF1, factor_key=TRUE)
long.fertile.data
long.fertile.data[long.fertile.data$Fertile == 1, 'Infertility.diagnosis'] = 'FERTILE'
long.fertile.data%>%head()

for (diag_i in unique(all.res$diag)){
  sig.res = all.res[all.res$P < 0.05, ]
  unique.cyto = sig.res[sig.res$diag == diag_i,'Cytokines'] %>% unique()
  sub.plot.data.infertile = long.fertile.data[grep(diag_i, long.fertile.data$Infertility.diagnosis),]
  sub.plot.data.infertile = sub.plot.data.infertile[sub.plot.data.infertile$cytokines %in% unique.cyto ,]
  sub.plot.data.infertile$Infertility.diagnosis = diag_i

  sub.plot.data.FERTILE = long.fertile.data[grep('FERTILE', long.fertile.data$Infertility.diagnosis),]
  sub.plot.data.FERTILE = sub.plot.data.FERTILE[sub.plot.data.FERTILE$cytokines %in% unique.cyto ,]
  sub.plot.data.FERTILE$Infertility.diagnosis = 'FERTILE'
  
  sub.plot.data = bind_rows(sub.plot.data.infertile, sub.plot.data.FERTILE)
  sub.plot.data$Infertility.diagnosis = factor(sub.plot.data$Infertility.diagnosis, levels = c('FERTILE', diag_i))
  
  cyto.chunks.diag = split(unique(sub.plot.data$cytokines), ceiling(seq_along(unique(sub.plot.data$cytokines))/20))
  
  pdf(paste0('dis_diag_res/', diag_i, '.pdf'), height = 10, width = 20, onefile = TRUE)
  for (i in cyto.chunks.diag){
    temp.fertile.data.diag = sub.plot.data[sub.plot.data$cytokines %in% i, ]
    print(
      ggboxplot(temp.fertile.data.diag, x = "Infertility.diagnosis", y = 'measurements',
                color = "Infertility.diagnosis", palette = "jco", add = "jitter") + 
        # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
        theme(legend.position="none")+
        facet_wrap(vars(cytokines), scales = 'free', nrow = 2)
    )
  }
  dev.off()
}
###plotting by  diagnosis###
