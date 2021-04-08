## testing out annova and tukey
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

fertile.data = bind_cols(infertility.data[, c('Cycle.Date', 'Age')], cytokines.data)
fertile.data
fertile.data$Age = log(fertile.data$Age, 2)


## remove columsn with single unique value ##
fertile.data = Filter(function(x)(length(unique(x))>1), fertile.data)

## factorize the cycle dates
fertile.data$Cycle.Date = factor(fertile.data$Cycle.Date)
fertile.data$Cycle.Date

res = sapply(colnames(fertile.data)[3:ncol(fertile.data)], function(x){
  lm.model = lm(as.formula(paste(x, '~Cycle.Date + Age')), data = fertile.data)
  
  # temp.wt.res = wald.test(b=coef(lm.model), Sigma=vcov(lm.model), Terms=c(2))
  # temp.wt.res$result$chi2[3]
  # summary(lm.model)
  # c(summary(lm.model)$coefficients[2,c(1,4)])
  # c(summary(lm.model)$coefficients[2,c(1,4)], temp.wt.res$result$chi2[3])
  c(summary(lm.model)$coefficients[2,c(1, 4)])
})

res = t(res) %>% data.frame()
res[,'padj'] = p.adjust(res[,2], method = 'fdr')

print(res[res[,3] < 0.05,])



annova.res = sapply(colnames(fertile.data)[3:ncol(fertile.data)], function(x){
  res.aov = aov(as.formula(paste(x, '~Cycle.Date + Age')), data = fertile.data)
  c(summary(res.aov)[[1]][1,5])
})

annova.res = annova.res %>% data.frame()
annova.res$padj =  p.adjust(annova.res[,1], method = 'fdr')
print(annova.res[annova.res$padj < 0.05,])

tukey.res = sapply(colnames(fertile.data)[3:ncol(fertile.data)], function(x){
  res.aov = aov(as.formula(paste(x, '~Cycle.Date + Age')), data = fertile.data)
  list(TukeyHSD(res.aov, which = 'Cycle.Date')$Cycle.Date[,4, drop = F])
})

tukey.res= do.call(cbind, tukey.res)
colnames(tukey.res) = colnames(fertile.data)[3:ncol(fertile.data)]
tukey.res = tukey.res %>% data.frame()
tukey.res = tukey.res %>% rownames_to_column(var = 'Date')
tukey.res = gather(tukey.res, cytokines, pvalues, IL8:CSF1, factor_key=TRUE)
tukey.res
tukey.res$padj =  p.adjust(tukey.res$pvalues, method = 'fdr')
sig.tukey = (tukey.res[tukey.res$padj < 0.05,])
table(sig.tukey$Date)

count.date.tukey = table(sig.tukey$Date)
sub.sig.tukey = sig.tukey[sig.tukey$Date %in% names(count.date.tukey[count.date.tukey > 2]),]
sub.sig.tukey = sub.sig.tukey[order(sub.sig.tukey$Date),]
sub.sig.tukey

fertile.data %>% 
  ggplot(aes(x=Cycle.Date,y=CCL28, color=Cycle.Date))+
  geom_boxplot(width=.5)+
  theme_classic()+
  # jittered text with geom_text
  theme(legend.position="none")

my_comparisons <- list(
                        c('14', '19'),
                        c('14', '20'),
                        c('14', '21'),
                        c('14', '28'))
ggboxplot(fertile.data[fertile.data$Cycle.Date != '12', ], x = "Cycle.Date", y = "CCL28",
          color = "Cycle.Date", palette = "jco", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(legend.position="none")

ggboxplot(fertile.data, x = "Cycle.Date", y = "TWEAK",
          color = "Cycle.Date", palette = "jco", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(legend.position="none")
