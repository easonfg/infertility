## looking to compare different time points
## using lm without outliers
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

date.counts = table(org.data$Cycle.Date)
sub.org.data = org.data[org.data$Cycle.Date %in% names(date.counts[date.counts > 2]),]
sub.org.data$Cycle.Date %>% table()

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
  cytokines.data[cytokines.data[,col.i] < LOD_x, col.i] = LOD_x
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

# fertile.data = bind_cols(infertility.data[, c('Cycle.Date', 'Age', 'Fertile', 'Infertility.diagnosis')], cytokines.data)
fertile.data = bind_cols(infertility.data[, c('Cycle.Date', 'Age', 'BMI', 'Fertile', 'Infertility.diagnosis')], cytokines.data)
fertile.data
fertile.data$Age = log(fertile.data$Age, 2)
fertile.data$BMI = log(fertile.data$BMI, 2)


## remove columsn with single unique value ##
fertile.data = Filter(function(x)(length(unique(x))>1), fertile.data)

## factorize the cycle dates and infertility diagnosis
fertile.data$Cycle.Date = factor(fertile.data$Cycle.Date)
fertile.data$Infertility.diagnosis = factor(fertile.data$Infertility.diagnosis)
fertile.data$Cycle.Date

timepoint.combos = combn(fertile.data$Cycle.Date %>% unique() %>% sort(), 2)


all.res = data.frame()
for (col.i in 1:ncol(timepoint.combos)) {
# for (col.i in 6:ncol(timepoint.combos)) {
  sub.fertile.data = fertile.data[fertile.data$Cycle.Date == timepoint.combos[1, col.i]|fertile.data$Cycle.Date == timepoint.combos[2,col.i],]
  
  res = sapply(colnames(sub.fertile.data)[grep('IL8', colnames(sub.fertile.data)):ncol(sub.fertile.data)], function(x){
    # lm.model = lm(as.formula(paste(x, '~Cycle.Date + Age + Fertile')), data = sub.fertile.data)
    # temp.data1 = sub.fertile.data[sub.fertile.data$Cycle.Date == timepoint.combos[1, col.i],c('Cycle.Date', 'Age', 'Fertile', 'Infertility.diagnosis', x)]
    temp.data1 = sub.fertile.data[sub.fertile.data$Cycle.Date == timepoint.combos[1, col.i],c('Cycle.Date', 'Age', 'BMI', 'Fertile', 'Infertility.diagnosis', x)]
    outliers <- boxplot(temp.data1[,x], plot=FALSE)$out
    temp.data1<- temp.data1[!(temp.data1[,x] %in% outliers),]
    
    # temp.data2 = sub.fertile.data[sub.fertile.data$Cycle.Date == timepoint.combos[2, col.i],c('Cycle.Date', 'Age', 'Fertile', 'Infertility.diagnosis', x)]
    temp.data2 = sub.fertile.data[sub.fertile.data$Cycle.Date == timepoint.combos[2, col.i],c('Cycle.Date', 'Age', 'BMI', 'Fertile', 'Infertility.diagnosis', x)]
    outliers <- boxplot(temp.data1[,x], plot=FALSE)$out
    temp.data2<- temp.data2[!(temp.data2[,x] %in% outliers),]
    
    temp.data = bind_rows(temp.data1, temp.data2)
    lm.model = lm(as.formula(paste(x, '~Cycle.Date + Age +BMI + Fertile+ Infertility.diagnosis')), data = temp.data)
    # lm.model = lm(as.formula(paste(x, '~Cycle.Date + Age + Fertile')), data = temp.data)
    # lm.model = lm(as.formula(paste(x, '~Cycle.Date + Age + Infertility.diagnosis')), data = temp.data)
    c(summary(lm.model)$coefficients[2,c(1, 4)])
  })
  res = t(res) %>% data.frame()
  res$timepoint = paste(timepoint.combos[,col.i], collapse = '-')
  all.res = bind_rows(all.res, res)
}

all.res = all.res[!is.na(all.res$Pr...t..),]
all.res$padj = p.adjust(all.res$Pr...t.., method = 'fdr')
all.res = all.res[order(all.res$timepoint),]

# print(all.res[all.res$padj < 0.05,])

all.res = all.res %>% rownames_to_column('cytokines')
all.res

all.res$cytokines = sapply(str_split(all.res$cytokines, '\\.\\.\\.'), function(x) x[1])

sig.all.res = (all.res[all.res$padj < 0.05,])
table(sig.all.res$cytokines) %>% sort()
print(sig.all.res)


### dot plot of p values
dot.plot.data = sig.all.res
dot.plot.data

gg.res = ggplot(subset(dot.plot.data, padj < 0.05), aes(timepoint, cytokines)) + 
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

# fertile.data %>% 
#   ggplot(aes(x=Cycle.Date,y=CCL28, color=Cycle.Date))+
#   geom_boxplot(width=.5)+
#   theme_classic()+
#   # jittered text with geom_text
#   theme(legend.position="none")

my_comparisons <- list(
  c('14', '12'),
  c('14', '19'),
  c('14', '20'),
  c('14', '21'),
  c('14', '28'))


table(sig.all.res$cytokines) %>% sort()
table(sig.all.res$timepoint)
cyto.int = 'CXCL10'
cyto.int = 'TRAIL'
cyto.int = 'CD244'
cyto.int = 'LIFR'
cyto.int = 'TNFRSF9'

sig.all.res[sig.all.res$cytokines == 'MCP3',]

long.fertile.data = gather(fertile.data, cytokines, measurements, IL8:CSF1, factor_key=TRUE)
cytokines.int =  unique(sig.all.res$cytokines)
long.fertile.data = long.fertile.data[long.fertile.data$cytokines %in% cytokines.int,]
long.fertile.data

## keep only the timepoints with signficance to graph
for (cyto in unique(long.fertile.data$cytokines)){
  unique.time = unique(unlist(str_split(sig.all.res[sig.all.res$cytokines == cyto, 'timepoint'], '-')))
  long.fertile.data = long.fertile.data[!(long.fertile.data$cytokines == cyto & !long.fertile.data$Cycle.Date %in% unique.time), ]
  
  # get rid of outliers
  for (time.i in unique.time){
   one.time.one.cyto = long.fertile.data[(long.fertile.data$cytokines == cyto & long.fertile.data$Cycle.Date %in% time.i), ]
   one.time.one.cyto
   outliers <- boxplot(one.time.one.cyto$measurements, plot=FALSE)$out
   long.fertile.data<- long.fertile.data[!(long.fertile.data$measurements %in% outliers),]

  }
  
}

###plotting by cytokines###
# cyto.chunks = split(cytokines.int, ceiling(seq_along(cytokines.int)/20))
# 
# pdf(paste0('traj_res/without_outliers.pdf'), height = 20, width = 20, onefile = TRUE)
# # pdf(paste0('traj_res/with_outliers.pdf'), height = 20, width = 20, onefile = TRUE)
# for (i in cyto.chunks){
#   temp.fertile.data = long.fertile.data[long.fertile.data$cytokines %in% i, ]
#   print(
#   ggboxplot(temp.fertile.data, x = "Cycle.Date", y = 'measurements',
#             color = "Cycle.Date", palette = "jco", add = "jitter") + 
#     # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
#     theme(legend.position="none")+
#     facet_wrap(vars(cytokines), scales = 'free', nrow = 2)
#   )
# }
# dev.off()
###plotting by cytokines###

###plotting by cycle date###
for (cd.i in unique(sig.all.res$timepoint)){
  unique.time.cd = unlist(str_split(cd.i, '-'))
  unique.cyto = sig.all.res[sig.all.res$timepoint == cd.i,'cytokines'] %>% unique()
  sub.plot.data = long.fertile.data[long.fertile.data$cytokines %in% unique.cyto & long.fertile.data$Cycle.Date %in% unique.time.cd, ]
  
  cyto.chunks.cd = split(unique(sub.plot.data$cytokines), ceiling(seq_along(unique(sub.plot.data$cytokines))/20))
  
  pdf(paste0('traj_res/', cd.i, '.pdf'), height = 10, width = 20, onefile = TRUE)
  for (i in cyto.chunks.cd){
    temp.fertile.data.cd = sub.plot.data[sub.plot.data$cytokines %in% i, ]
    print(
      ggboxplot(temp.fertile.data.cd, x = "Cycle.Date", y = 'measurements',
                color = "Cycle.Date", palette = "jco", add = "jitter") + 
        # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
        theme(legend.position="none")+
        facet_wrap(vars(cytokines), scales = 'free', nrow = 2)
    )
  }
  dev.off()
}
###plotting by cycle date###


### plotting frequencies of data
org.freq = table(org.data$Cycle.Date) %>% data.frame()
org.freq$Var1 = as.character(org.freq$Var1)

plot(org.freq$Var1, org.freq$Freq, type = 'h', ylim=c(0,50), xlab = 'Cycle Date', ylab = 'Count')
axis(side=1, at=org.freq$Var1, labels=org.freq$Var1)
text(org.freq$Var1, org.freq$Freq,  org.freq$Freq,
     cex=1.05, pos=3,col="red")

fertile.data.freq = table(fertile.data$Cycle.Date) %>% data.frame()
fertile.data.freq$Var1 = as.character(fertile.data.freq$Var1)

plot(fertile.data.freq$Var1, fertile.data.freq$Freq, type = 'h', ylim=c(0,50), xlab = 'Cycle Date', ylab = 'Count')
axis(side=1, at=fertile.data.freq$Var1, labels=fertile.data.freq$Var1)
text(fertile.data.freq$Var1, fertile.data.freq$Freq,  fertile.data.freq$Freq,
     cex=1.05, pos=3,col="red")

sig.all.res.freq = table(sig.all.res$timepoint) %>% data.frame()
# sig.all.res.freq$Var1 = as.character(sig.all.res.freq$Var1)
sig.all.res.freq

plot(sig.all.res.freq$Var1, sig.all.res.freq$Freq, type = 'h', ylim = c(0,60), xlab = 'Cycle Date Comparison', ylab = 'Count')
axis(side=1, at=sig.all.res.freq$Var1, labels=sig.all.res.freq$Var1)
text(sig.all.res.freq$Var1, sig.all.res.freq$Freq,  sig.all.res.freq$Freq,
     cex=1.05, pos=3,col="red")
### plotting frequencies of data

