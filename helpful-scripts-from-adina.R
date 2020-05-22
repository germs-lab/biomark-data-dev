#devtools::install_github("sjackman/uniqtag")
library(plyr)
library(uniqtag)
library(ggplot2)

setwd("~/Box Sync/biomark")
raw_data_file = "A1.1 Raw Data.csv"
raw_data <- read.csv(raw_data_file, sep=",", row.names=1, header = TRUE)
dim(raw_data)
head(raw_data)
raw_data$full_assay <- paste(raw_data$Name,"-", raw_data$Name.1, sep="")
raw_data$unique_id <- make_unique(raw_data$full_assay)

# Metadata fie has to have a header names "sample_type" where you specify standards or env sample
# You should have metadata for both samples and standards
# It would be better if the "standard sample" had the appropriate dilution
meta_data_file = "metadata.csv"
meta <- read.csv(meta_data_file, header=TRUE)

# Check the NTCs
ct_threshold_ntc = 20 #Should be set by you
ntc <- subset(raw_data, Name == "NTC")
summary(ntc$Value)
check_ntc <- subset(ntc, Value > 20 & Value < 999)
check_ntc$Value
data_nontc <- subset(raw_data, Name != "NTC")

merged <- merge(data_nontc, meta, by = "Name")


#Standards Workflow

standards <- subset(merged, sample_type == "standard")
gene_list <- unique(standards$Name)

for (x in 1:length(gene_list)) {
  each_gene <- subset(standards, Name == gene_list[x])
  standard_only = paste(gene_list[x], '-', gene_list[x], sep='')
  each_gene_filtered = subset(each_gene, full_assay == standard_only)
  #print(each_gene_filtered)
  each_gene_badfiltered = subset(each_gene, full_assay != standard_only)
  each_gene_badfiltered = subset(each_gene_badfiltered, Value < 999) # This can be used to filter
  #print(each_gene_badfiltered)
  f <- ddply(each_gene_filtered, .(rConc, full_assay), summarise, MEAN = mean(Value), SE=sd(Value)/sqrt(length(Value)))
  limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
  f$log <- log10(f$rConc)
  p = ggplot(f, aes_string(x="log", y="MEAN", color="full_assay"))
  p = p + geom_point() + geom_errorbar(position=position_dodge(), limits)+theme(axis.text.x = element_text(angle = 90)) 
  p
  ggsave(paste0(x,'standard.pdf'), device="pdf")
  f <- ddply(each_gene_badfiltered, .(rConc, full_assay), summarise, MEAN = mean(Value), SE=sd(Value)/sqrt(length(Value)))
  limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
  f$log <- log10(f$rConc)
  p = ggplot(f, aes_string(x="log", y="MEAN", color="full_assay"))
  p = p + geom_point() + geom_errorbar(position=position_dodge(), limits)+theme(axis.text.x = element_text(angle = 90)) 
  p
  ggsave(paste0(x,'nonstandard.pdf'), device="pdf")
         
}

#Samples Workflow

samples <- subset(merged, sample_type == "env") #all samples
#average across technical replicates
f <- ddply(samples, .(full_assay, location, plot, time, matrix, Name, Name.1), summarise, MEAN = mean(Value), SE=sd(Value)/sqrt(length(Value)))
#genes that are detected #stats are done prior to removing bad measurements
f_detects <- subset(f, MEAN < 999)
print("Number of Genes with detects:  ")
length(unique(f_detects$Name.1))
unique(f_detects$Name.1)
#Analysis likely useful by gene or by sample
gene_list <- unique(f$Name.1)
for (x in 1:length(gene_list)) {
  each_gene <- subset(f, Name.1 == gene_list[x])
  limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
  f2 <- subset(each_gene, MEAN < 25) #Can set as appropriate - filter by MEAN technical rep Ct
  p = ggplot(f2, aes_string(x="full_assay", y="MEAN", color="matrix"))
  p = p + geom_point() + geom_errorbar(position=position_dodge(), limits)+theme(axis.text.x = element_text(angle = 90)) 
  p
  write.csv(file=paste0(x,".sample.csv"), each_gene, quote=FALSE, row.names=FALSE)
  ggsave(paste0(x,'sample.pdf'), device="pdf")
}
for (x in 1:length(gene_list)) {
  each_gene <- subset(f, Name.1 == gene_list[x])

  #Note this is going to take the average of the averages, not sure if this is to be done
  f_plot <-ddply(each_gene, .(plot, Name.1), summarise, MEAN2 = mean(MEAN), SE=sd(MEAN)/sqrt(length(MEAN)))
  f2 <- subset(f_plot, MEAN2 < 25) #Can set as appropriate - filter by MEAN technical rep Ct
  limits<-aes(ymin=MEAN2-SE, ymax=MEAN2+SE)
  p = ggplot(f2, aes_string(x="plot", y="MEAN2", color="Name.1"))
  p = p + geom_point() + geom_errorbar(position=position_dodge(), limits)+theme(axis.text.x = element_text(angle = 90)) 
  p
  write.csv(file=paste0(x,".byplot.csv"), each_gene, quote=FALSE, row.names=FALSE)
  ggsave(paste0(x,'byplot.pdf'), device="pdf")
}

each_gene <- subset(f, Name.1 == gene_list[x])

#EVERYTHING!
#Note this is going to take the average of the averages, not sure if this is to be done
f_plot <-ddply(f, .(plot, Name.1), summarise, MEAN2 = mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN)))
f2 <- subset(f_plot, MEAN2 < 25) #Can set as appropriate - filter by MEAN technical rep Ct
limits2<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
p = ggplot(f2, aes_string(x="plot", y="MEAN2", color="Name.1"))
p + geom_point() + geom_errorbar(width=0.25, limits2)+theme(axis.text.x = element_text(angle = 90)) 


  