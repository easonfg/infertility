## using lm with outliers
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
cytokines.data = sapply(colnames(cytokines.data), function(x){
  LOD = cytokines.data[nrow(cytokines.data),x]
  cytokines.data[cytokines.data[,x] < LOD,x] = LOD
  cytokines.data[,x]
})

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

fertile.data = bind_cols(infertility.data[, c('Cycle.Date', 'Age', 'Fertile', 'Infertility.diagnosis')], cytokines.data)
fertile.data
fertile.data$Age = log(fertile.data$Age, 2)


## remove columsn with single unique value ##
fertile.data = Filter(function(x)(length(unique(x))>1), fertile.data)

## factorize the cycle dates and infertility diagnosis
fertile.data$Cycle.Date = factor(fertile.data$Cycle.Date)
fertile.data$Infertility.diagnosis = factor(fertile.data$Infertility.diagnosis)
fertile.data$Cycle.Date

timepoint.combos = combn(fertile.data$Cycle.Date %>% unique() %>% sort(), 2)


all.res = data.frame()
for (col.i in 1:ncol(timepoint.combos)) {
  sub.fertile.data = fertile.data[fertile.data$Cycle.Date == timepoint.combos[1, col.i]|fertile.data$Cycle.Date == timepoint.combos[2,col.i],]
  
  res = sapply(colnames(sub.fertile.data)[4:ncol(sub.fertile.data)], function(x){
    # lm.model = lm(as.formula(paste(x, '~Cycle.Date + Age + Fertile')), data = sub.fertile.data)
    lm.model = lm(as.formula(paste(x, '~Cycle.Date + Age + Fertile+ Infertility.diagnosis')), data = sub.fertile.data)
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
ggboxplot(fertile.data, x = "Cycle.Date", y = cyto.int,
          color = "Cycle.Date", palette = "jco", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(legend.position="none")+
  facet_wrap(vars(subtype), scales = 'free', nrow = 2)


sig.all.res[sig.all.res$cytokines == 'MCP3',]

long.fertile.data = gather(fertile.data, cytokines, measurements, IL8:CSF1, factor_key=TRUE)
cytokines.int =  unique(sig.all.res$cytokines)
long.fertile.data = long.fertile.data[long.fertile.data$cytokines %in% cytokines.int,]
long.fertile.data

for (cyto in unique(long.fertile.data$cytokines)){
  unique.time = unique(unlist(str_split(sig.all.res[sig.all.res$cytokines == cyto, 'timepoint'], '-')))
  long.fertile.data = long.fertile.data[!(long.fertile.data$cytokines == cyto & !long.fertile.data$Cycle.Date %in% unique.time), ]
}
long.fertile.data1 = long.fertile.data[long.fertile.data$cytokines %in% cytokines.int[1:18], ]
long.fertile.data2 = long.fertile.data[long.fertile.data$cytokines %in% cytokines.int[19:length(cytokines.int)], ]



ggboxplot(long.fertile.data1, x = "Cycle.Date", y = 'measurements',
          color = "Cycle.Date", palette = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(legend.position="none")+
  facet_wrap(vars(cytokines), scales = 'free', nrow = 2)

ggboxplot(long.fertile.data2, x = "Cycle.Date", y = 'measurements',
          color = "Cycle.Date", palette = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(legend.position="none")+
  facet_wrap(vars(cytokines), scales = 'free', nrow = 2)
