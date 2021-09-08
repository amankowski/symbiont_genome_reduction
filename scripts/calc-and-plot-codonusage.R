library(seqinr)
library(readr)
library(janitor)
library(tidyverse)
library(reshape2)
library(Biostrings)
library(patchwork)
library(coRdon)

setwd("../with_genes")

codon.usage<-""
file.names<-dir(".", pattern=".fasta")
for (i in 1:length(file.names)){
  LIB<-(gsub("(.*)_.*.fasta", "\\1", file.names[i]))
  fasta<-readSet(file=file.names[i])
  fasta.names<-names(fasta)
  ct<-codonTable(fasta)
  cc<-codonCounts(ct)
  codon.tmp<-as.data.frame(rbind(cc, colSums(cc)))
  codon.tmp$LIB<-rep(LIB, nrow(codon.tmp))
  codon.tmp$GENE<-c(fasta.names, "SUM")
  codon.usage<-rbind(codon.usage, codon.tmp)
}

codon.usage<-codon.usage[-1,]
codon.usage.melt<-melt(codon.usage, id.vars=c("LIB", "GENE"))
codon.usage.melt$ORGANISM<-ifelse(grepl("Gamma4", codon.usage.melt$LIB), "Gamma4", "Kentron")

codon.usage.melt$LIB<-as.factor(codon.usage.melt$LIB)
codon.usage.melt$GENE<-as.factor(codon.usage.melt$GENE)
codon.usage.melt$ORGANISM<-as.factor(codon.usage.melt$ORGANISM)
codon.usage.melt$variable<-as.factor(codon.usage.melt$variable)
codon.usage.melt$value<-as.numeric(codon.usage.melt$value)


codon.usage.plot<-ggplot(subset(codon.usage.melt, GENE %in% "SUM"), aes(x=variable, y=value, colour=ORGANISM)) + 
  geom_point() + 
  geom_line() + 
  #geom_boxplot() +
  #scale_fill_manual(values=c("Kentron"="#069995", "Gamma4"="#9ababa")) +
  scale_colour_manual(values=c("Kentron"="#069995", "Gamma4"="#9ababa")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))
  
ggsave(codon.usage.plot, file="codon.usage.eps",width=180, height=210, units="mm")