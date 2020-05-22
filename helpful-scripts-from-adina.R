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
ct_threshold_ntc = 20
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
  each_gene_badfiltered = subset(each_gene_badfiltered, Value < 999)
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
