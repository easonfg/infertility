### hormones correlation
rm(list=ls())

library(ggpubr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(aod)

org.data = read.csv('../cleaned.infertile.data.csv')
org.data$Infertility.diagnosis = toupper(org.data$Infertility.diagnosis)

# org.data = org.data[grep('IVF', org.data$Cycle.Type),]


LOD = org.data[org.data$X == 'LOD', , drop = F]
LOD


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

fertile.data = bind_cols(infertility.data[, c('Cycle.Date', 'Cycle.Type', 'AMH',
                                              # "ESTRADIOL", "Progestrone", "TESTOSTERONE", 'T4', 
                                              'CORTISOL', 'ESTRADIOL',	'Progestrone',	'TESTOSTERONE',	'T3',	'T4',
                                              'Age', 'BMI', 'Fertile', 'Infertility.diagnosis')], cytokines.data)
fertile.data
fertile.data$Age = log(fertile.data$Age, 2)
fertile.data$BMI = log(fertile.data$BMI, 2)
fertile.data$AMH = log(fertile.data$AMH, 2)


## remove columsn with single unique value ##
fertile.data = Filter(function(x)(length(unique(x))>1), fertile.data)

## factorize the cycle dates and infertility diagnosis
# fertile.data$Cycle.Date = factor(fertile.data$Cycle.Date)
# fertile.data$Cycle.Type = factor(fertile.data$Cycle.Type)
fertile.data$Cycle.Type
fertile.data$Infertility.diagnosis = factor(fertile.data$Infertility.diagnosis)
fertile.data$Cycle.Date

quartiles = quantile(fertile.data$TESTOSTERONE)
lower.quart.testosterone = fertile.data[fertile.data$TESTOSTERONE < quartiles[2],]
upper.quart.testosterone = fertile.data[fertile.data$TESTOSTERONE > quartiles[4],]
lower.quart.testosterone$quartile.testosterone = 'lower'
upper.quart.testosterone$quartile.testosterone = 'upper'

sub.fertile.data = bind_rows(lower.quart.testosterone, upper.quart.testosterone)
sub.fertile.data

all.res = data.frame()
# for (diag in c('AMH')){
# for (diag in c("ESTRADIOL", "Progestrone", "TESTOSTERONE", 'T4')){
#for (diag in c('CORTISOL', 'ESTRADIOL',	'Progestrone',	'TESTOSTERONE',	'T3',	'T4')){
for (diag in c('quartile.testosterone')){
  
  res = sapply(colnames(sub.fertile.data)[grep('IL8', colnames(sub.fertile.data)):grep('CSF1', colnames(sub.fertile.data))], function(x){
    
    # lm.model = lm(as.formula(paste(x, paste0('~', diag, ' + Age +BMI + Fertile+ Cycle.Date + Infertility.diagnosis'))), data = sub.fertile.data)
    lm.model = lm(as.formula(paste(x, paste0('~', diag, ' + Age +BMI + Cycle.Date + Cycle.Type+ Infertility.diagnosis'))), data = sub.fertile.data)
    # lm.model = lm(as.formula(paste(x, paste0('~', diag, ' + Age +BMI + Cycle.Date' ))), data = sub.fertile.data)
    c(summary(lm.model)$coefficients[2,c(1, 4)])
  })
  
  res = t(res) %>% data.frame()
  res$hormone = diag
  
  all.res = bind_rows(all.res, res)
}

all.res
all.res[,'padj'] = p.adjust(all.res[,2], method = 'fdr')
print(all.res[all.res$padj < 0.1,])
print(all.res[all.res$padj < 0.05,])
print(all.res[all.res$P < 0.05,])

# sub_padj = c()
# for (diag in all.res$hormone %>% unique()){
#   sub_padj = c(sub_padj, p.adjust(all.res[all.res$hormone == diag, 'P'], method = 'fdr'))
# }
# all.res$sub.padj = sub_padj
# all.res[all.res$sub.padj < 0.05,]


## change rownames to column
all.res = all.res %>% rownames_to_column('Cytokines')
all.res$Cytokines = sapply(strsplit(all.res$Cytokines, '\\.\\.\\.'), function(x) x[1])


#res = res[!is.na(res$Pr...t..),]
#res$padj = p.adjust(res$Pr...t.., method = 'fdr')




# all.res$cytokines = sapply(str_split(all.res$cytokines, '\\.\\.\\.'), function(x) x[1])

# sig.adj.all.res = (all.res[all.res$padj < 0.1,])
# sig.raw.all.res = (all.res[all.res$Pr...t.. < 0.05,])


### dot plot of p values
dot.plot.data = all.res
colnames(dot.plot.data)[3] = 'raw.P'
# dot.plot.data$Cytokines = factor(dot.plot.data$Cytokines, levels = dot.plot.data$Cytokines[order(dot.plot.data$Estimate)])

# gg.res = ggplot(subset(dot.plot.data, padj < 0.1), aes(hormone, Cytokines)) + 
gg.res = ggplot(subset(dot.plot.data, raw.P < 0.05), aes(hormone, Cytokines)) + 
  geom_point(aes(size = -log10(padj), fill = Estimate), 
             colour = 'black', shape = 21) + 
  geom_text(aes(label=round(Estimate,2)),size = 3, hjust=-0.5, vjust=0.5) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  labs(size = '-log10(raw.P.val)', fill = 'coefficients') +
  theme_bw()
print(gg.res)

gg.res = ggplot(subset(dot.plot.data, padj < 0.1), aes(hormone, Cytokines)) + 
  geom_point(aes(size = -log10(padj), fill = Estimate), 
             colour = 'black', shape = 21) + 
  geom_text(aes(label=round(Estimate,2)),size = 3, hjust=-0.5, vjust=0.5) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  labs(size = '-log10(FDR)', fill = 'coefficients') +
  theme_bw()
print(gg.res)
### dot plot of p values



sig.all.res = (all.res[all.res$Pr...t.. < 0.05,])
sig.all.res%>%head()
long.fertile.data = gather(sub.fertile.data, cytokines, measurements, IL8:CSF1, factor_key=TRUE)
cytokines.int =  unique(sig.all.res$Cytokines)
cytokines.int
# long.fertile.data = long.fertile.data[long.fertile.data$cytokines %in% cytokines.int,]
long.fertile.data

# long.fertile.data$cytokines = factor(long.fertile.data$cytokines, levels = dot.plot.data$cytokines[order(dot.plot.data$Estimate)])



sub.plot.data = long.fertile.data
sub.plot.data

# cyto.chunks.cd = split(levels(sub.plot.data$cytokines), ceiling(seq_along(levels(sub.plot.data$cytokines))/11))

# for (hormone in c("ESTRADIOL", "Progestrone", "TESTOSTERONE", 'T4')){
# for (hormone in c('CORTISOL', 'ESTRADIOL',	'Progestrone',	'TESTOSTERONE',	'T3',	'T4')){
for (hormone in c('quartile.testosterone')){
  # cyto.int = sig.all.res[sig.all.res$hormone %in% hormone, 'Cytokines']
  cyto.int = long.fertile.data$cytokines %>% unique()
  cyto.chunks.cd = split(cyto.int, ceiling(seq_along(cyto.int)/10))
  print(hormone)
  pdf(paste0(hormone, '_regression', '.pdf'), height = 10, width = 20, onefile = TRUE)
  for (i in cyto.chunks.cd){
    temp.fertile.data.diag = long.fertile.data[long.fertile.data$cytokines %in% i, ]
    print(i)
    print(
      ggboxplot(temp.fertile.data.diag, x = "cytokines", y = 'measurements',
                color = "quartile.testosterone", palette = "jco", add = "jitter") + 
        # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
        #theme(legend.position="none")+
        facet_wrap(vars(cytokines), scales = 'free', nrow = 2)
    )
  }
  dev.off()
}
# dev.off()

