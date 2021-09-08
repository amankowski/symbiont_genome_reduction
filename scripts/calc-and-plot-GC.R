#!/usr/bin/env Rscript

library(seqinr)
library(readr)
library(janitor)
library(tidyverse)
library(reshape2)
library(Biostrings)
library(patchwork)
library(plyr)
library(dplyr)

clades<-read.csv("../clades.csv", sep="\t", header=T)
clade.list<-c("Gamma4", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11")

setwd("./data/combined")

file.names<-dir(".", pattern=".fasta")
GC.combined<-""

for (i in 1:length(file.names)){
  LIB<-(gsub("(.*)_.*.fasta", "\\1", file.names[i]))
  CODE<-(gsub(".*_(.*).fasta", "\\1", file.names[i]))
  mysequence<-s2c(read_file(file.names[i]))
  GC<-GC(mysequence)
  GC.tmp<-data.frame(LIB, CODE, GC)
  GC.combined<-rbind(GC.combined, GC.tmp)
} 

GC.combined<-GC.combined[-1,]

GC.combined$ORGANISM<-ifelse(grepl("Gamma4", GC.combined$LIB), "Gamma4", "Kentron")
GC.combined$LIB<-as.factor(GC.combined$LIB)
GC.combined$CODE<-as.factor(GC.combined$CODE)
GC.combined$ORGANISM<-as.factor(GC.combined$ORGANISM)
GC.combined$GC<-as.numeric(GC.combined$GC)

write.table(GC.combined, file="GC.combined.csv", quote=F, sep=",", row.names=T, col.names=T)

GC.combined.host<-ggplot(GC.combined, aes(x=ORGANISM, y=GC, fill=CODE)) + 
  geom_boxplot(outlier.size=1) +
  scale_fill_manual(values=c("coding"="#069995", "non-coding"="#9ababa")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))

ggsave(GC.combined.host, file="../../plots/GC.c-vs-nc.host.eps", width=180, height=210, units="mm")

GC.combined.clades<-merge(GC.combined, clades, by="LIB")
GC.combined.clades$CLADE<-as.factor(GC.combined.clades$CLADE)
GC.combined.clades$CLADE<-factor(GC.combined.clades$CLADE, levels=clade.list)

GC.combined.clade<-ggplot(GC.combined.clades, aes(x=CLADE, y=GC, fill=CODE)) + 
  geom_boxplot(outlier.size=1) +
  scale_fill_manual(values=c("coding"="#069995", "non-coding"="#9ababa")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))

ggsave(GC.combined.clade, file="../../plots/GC.c-vs-nc.clade.eps", width=180, height=210, units="mm")

GC.combined.plot<-(GC.combined.host | GC.combined.clade)


setwd("../with_genes")

file.names<-dir(".", pattern=".fasta")
GC.single_genes<-data.frame("LIB", "GENE", "CODE", 0.5, 1)
colnames(GC.single_genes)<-c("LIB", "GENE", "CODE", "GC", "LENGTH")

for (i in 1:length(file.names)){
  LIB<-(gsub("(.*)_.*.fasta", "\\1", file.names[i]))
  fasta<-readDNAStringSet(file.names[i])
  seqs<-paste(fasta)
  fasta.names<-names(fasta)
  for (j in 1:length(seqs)){
    GENE<-fasta.names[j]
    mysequence<-s2c(seqs[j])
    GC<-GC(mysequence)
    GC1<-GC1(mysequence)
    GC2<-GC2(mysequence)
    GC3<-GC3(mysequence)
    LENGTH<-length(mysequence)
    GC.single_genes<-GC.single_genes %>% add_row(LIB=LIB, GENE=GENE, CODE="GC", GC=GC, LENGTH=LENGTH)
    GC.single_genes<-GC.single_genes %>% add_row(LIB=LIB, GENE=GENE, CODE="GC1", GC=GC1, LENGTH=LENGTH)
    GC.single_genes<-GC.single_genes %>% add_row(LIB=LIB, GENE=GENE, CODE="GC2", GC=GC2, LENGTH=LENGTH)
    GC.single_genes<-GC.single_genes %>% add_row(LIB=LIB, GENE=GENE, CODE="GC3", GC=GC3, LENGTH=LENGTH)
}}

GC.single_genes<-GC.single_genes[-1,]
GC.single_genes$ORGANISM<-ifelse(grepl("Gamma4", GC.single_genes$LIB), "Gamma4", "Kentron")
GC.single_genes$LIB<-as.factor(GC.single_genes$LIB)
GC.single_genes$GENE<-as.factor(GC.single_genes$GENE)
GC.single_genes$CODE<-as.factor(GC.single_genes$CODE)
GC.single_genes$ORGANISM<-as.factor(GC.single_genes$ORGANISM)
GC.single_genes$GC<-as.numeric(GC.single_genes$GC)
GC.single_genes$LENGTH<-as.numeric(GC.single_genes$LENGTH)

write.table(GC.single_genes, file="GC.single_genes.csv", quote=F, sep=",", row.names=T, col.names=T)

GC.single_genes.host<-ggplot(GC.single_genes, aes(x=ORGANISM, y=GC, fill=CODE)) + 
  geom_boxplot(outlier.size=1) +
  scale_fill_manual(values=c("GC"="#CCCCCC", "GC1"="#9ababa", "GC2"="#069995", "GC3"="#004542")) +
  theme_minimal()  +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))

ggsave(GC.single_genes.host, file="../../plots/GC.positions.host.eps", width=180, height=210, units="mm")

GC.single_genes.clades<-merge(GC.single_genes, clades, by="LIB")
GC.single_genes.clades$CLADE<-as.factor(GC.single_genes.clades$CLADE)
GC.single_genes.clades$CLADE<-factor(GC.single_genes.clades$CLADE, levels=clade.list)

GC.single_genes.clade<-ggplot(GC.single_genes.clades, aes(x=CLADE, y=GC, fill=CODE)) + 
  geom_boxplot(outlier.size=0.5) +
  scale_fill_manual(values=c("GC"="#CCCCCC", "GC1"="#9ababa", "GC2"="#069995", "GC3"="#004542"))  +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))

ggsave(GC.single_genes.clade, file="../../plots/GC.positions.clade.eps", width=180, height=210, units="mm")

GC.over.length<-ggplot(subset(GC.single_genes, CODE %in% c("GC", "GC3")), aes(x=LENGTH, y=GC, colour=ORGANISM, shape=CODE), size=5) +
  geom_point() +
  scale_colour_manual(values=c("Kentron"="#9ababa", "Gamma4"="#069995"))  +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))

ggsave(GC.over.length, file="../../plots/GC.over.length.eps", height=210, width=180, units="mm")

GC.overall<-rbind(GC.combined, subset(GC.single_genes, select=-c(GENE, LENGTH)))
GC.overall.clades<-merge(GC.overall, clades, by="LIB")

codes<-c("non-coding", "coding", "GC", "GC1", "GC2", "GC3")
GC.overall.clades$CDOE<-as.factor(GC.overall.clades$CODE)
GC.overall.clades$CODE<-factor(GC.overall.clades$CODE, levels=codes)

GC.overall.clades$CLADE<-as.factor(GC.overall.clades$CLADE)
GC.overall.clades$CLADE<-factor(GC.overall.clades$CLADE, levels=clade.list)

GC.overall.host<-ggplot(subset(GC.overall.clades, CODE %in% c("non-coding", "GC", "GC3")), aes(x=ORGANISM, y=GC, fill=CODE)) + 
  geom_boxplot(outlier.size=1) +
  scale_fill_manual(values=c("non-coding"="#CCCCCC", "GC"="#9ababa", "GC3"="#069995")) +
  theme_classic() #+
  #theme(axis.text.x=element_text(angle=90, hjust=1)) +
  #theme(text=element_text(size=5))

ggsave(GC.overall.host, file="GC.overall.host.eps", width=140, height=140, units="mm")

GC.overall.clade<-ggplot(subset(GC.overall.clades, CODE %in% c("non-coding", "GC", "GC3")), aes(x=CLADE, y=GC, fill=CODE)) + 
  geom_boxplot(outlier.size=1) +
  scale_fill_manual(values=c("non-coding"="#CCCCCC", "GC"="#9ababa", "GC3"="#069995")) +
  theme_classic() #+
#theme(axis.text.x=element_text(angle=90, hjust=1)) +
#theme(text=element_text(size=5))

ggsave(GC.overall.clade, file="GC.overall.clade.eps", width=140, height=140, units="mm")

MEAN<-""
for (i in 1:length(clade.list)){
  clade<-clade.list[i]
  NC.mean<-mean(subset(GC.overall.clades, (CODE %in% "non-coding") & (CLADE %in% clade))$GC)
  GC.mean<-mean(subset(GC.overall.clades, (CODE %in% "GC") & (CLADE %in% clade))$GC)
  GC3.mean<-mean(subset(GC.overall.clades, (CODE %in% "GC3") & (CLADE %in% clade))$GC)
  mean.tmp<-data.frame(clade, NC.mean, GC.mean, GC3.mean)
  MEAN<-rbind(MEAN, mean.tmp)
}

MEAN<-MEAN[-1,]
MEAN$clade<-as.factor(MEAN$clade)
MEAN$NC.mean<-as.numeric(MEAN$NC.mean)
MEAN$GC.mean<-as.numeric(MEAN$GC.mean)
MEAN$GC3.mean<-as.numeric(MEAN$GC3.mean)

MEAN$clade<-factor(MEAN$clade, levels=clade.list)
MEAN$diff.NCGC<-100*(MEAN$NC.mean-MEAN$GC.mean)
MEAN$diff.GCGC3<-100*(MEAN$GC3.mean-MEAN$GC.mean)

mean.plot<-ggplot(MEAN) + 
  geom_point(aes(x=clade, y=diff.GCGC3), colour="#9ababa", size=5) + 
  geom_point(aes(x=clade, y=diff.NCGC), colour="#069995", size=5) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black") + 
  theme(axis.line.y=element_line(colour="black", size=.6, linetype="solid")) +
  theme(panel.grid=element_blank())

ggsave(mean.plot, file="mean.plot.eps", width=140, height=140, units="mm")

gained<-read.csv("../../../../gain-and-loss/Gamma4.gained-genes")
GC.single_genes$GAIN<-ifelse(GC.single_genes$GENE %in% gained$V1, "yes", "no")

ggplot(subset(single, (CODE %in% c("GC", "GC3")) & (ORGANISM %in% "Gamma4")),  aes(x=CODE, y=GC)) + 
  geom_boxplot(outlier.size=0.5, fill="#CCCCCC") +
  geom_point(data=subset(single,(CODE %in% c("GC", "GC3")) & (ORGANISM %in% "Gamma4") & (GAIN %in% "yes")), color="#069995") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))

single.clades$CLADE<-as.factor(single.clades$CLADE)
single.clades$CLADE<-factor(single.clades$CLADE, levels=clade.list)

mu<-ddply(single.clades, "CLADE", summarise, grp.mean=mean(LENGTH))
mu2<-ddply(single.clades, "ORGANISM", summarise, grp.mean=mean(LENGTH))
mu3<-ddply(single.clades, "LIB", summarise, grp.mean=mean(LENGTH))


ggplot(subset(single.clades, (CODE %in% "GC")),  aes(x=LENGTH)) + 
  geom_histogram(aes(fill=ORGANISM), binwidth=100) + 
  scale_fill_manual(values=c("Kentron"="#9ababa", "Gamma4"="#069995"))  +
  facet_grid(LIB ~ .) +
  geom_vline(data=mu3, aes(xintercept=grp.mean), colour="black", linetype="dashed") +
  xlim(0, 3000) + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))

ggplot(single.clades.genome.mu.counts,  aes(x=GENE_COUNTS, y=grp.mean, colour=ORGANISM), size=5) + 
  geom_point() + 
  scale_colour_manual(values=c("Kentron"="#9ababa", "Gamma4"="#069995")) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))

ggplot(single.clades.genome.mu.counts,  aes(y=GENE_COUNTS, x=GENOME, colour=CLADE, size=grp.mean), alpha=.2) + 
  geom_point() + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(text=element_text(size=5))
