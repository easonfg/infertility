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

org.data = org.data[grep('IVF', org.data$Cycle.Type),]


LOD = org.data[org.data$X == 'LOD', , drop = F]
LOD

## add cycle name according to cycle date
org.data = add_column(org.data, Cycle.Name = 'O', .after = 'Cycle.Date')
org.data[(org.data$Cycle.Date <= 5) & (!is.na(org.data$Cycle.Date)), 'Cycle.Name'] = 'M'
org.data[(org.data$Cycle.Date > 14) & (!is.na(org.data$Cycle.Date)), 'Cycle.Name'] = 'L'
org.data[(org.data$Cycle.Date <= 13) & (org.data$Cycle.Date > 5) & (!is.na(org.data$Cycle.Date)), 'Cycle.Name'] = 'F'
org.data$Cycle.Name

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

fertile.data = bind_cols(infertility.data[, c('Cycle.Date', 'Cycle.Type', 'AMH', 'Cycle.Name',
                                              # "ESTRADIOL", "Progestrone", "TESTOSTERONE", 'T4', 
                                              'CORTISOL', 'ESTRADIOL',	'Progestrone',	'TESTOSTERONE',	'T3',	'T4',
                                              'Age', 'BMI', 'Fertile', 'Infertility.diagnosis')])
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


long.fertile.data = gather(fertile.data, hormones, measurements, CORTISOL:T4, factor_key=TRUE)
# long.fertile.data = long.fertile.data[long.fertile.data$cytokines %in% cytokines.int,]
long.fertile.data$Cycle.Name = factor(long.fertile.data$Cycle.Name, levels = c('M', 'F', 'O', 'L'))

# long.fertile.data$cytokines = factor(long.fertile.data$cytokines, levels = dot.plot.data$cytokines[order(dot.plot.data$Estimate)])



sub.plot.data = long.fertile.data
sub.plot.data

# cyto.chunks.cd = split(levels(sub.plot.data$cytokines), ceiling(seq_along(levels(sub.plot.data$cytokines))/11))

# for (hormone in c("ESTRADIOL", "Progestrone", "TESTOSTERONE", 'T4')){

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- long.fertile.data$Cycle.Date %>% unique() %>% length()
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

cycle.var = 'Cycle.Date'
print(
ggboxplot(long.fertile.data, x = cycle.var, y = 'measurements',
          color = cycle.var, palette = "jco", add = "jitter") + 
          scale_color_manual(values = mycolors) +
  # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  #theme(legend.position="none")+
  facet_wrap(vars(hormones), scales = 'free', nrow = 2)

)