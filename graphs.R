# Claire Malley
# Graphs for ASHG2017 poster
# Compiled together on September 14

library(data.table)
library(viridis)
library(ggplot2)
library(ggbio)
library(venneuler)
library(eulerr)
options(scipen = 999)

## EDC coverage graph -----
#edc data import
edc <- fread("ADRN_EDC_799_platypus_maxvar8_minreads5.vcf_PASS.TR7.GQ20.NR7.recode.vcf", sep="\t", header=T, skip=51L, stringsAsFactors = F)

splitInfoCol <- function(x){
  dt <- x
  dt <- dt[FILTER=="PASS"] # PASS only sites please
  dt <- dt[,c("POS", "INFO")]
  dt <- dt[,c("BRF", "FR", "HP", "HapScore", "MGOF", "MMLQ", "MQ", "NF", "NR", "PP", "QD", "SC", "SbPval", "Source", "TC", "TCF", "TCR", "TR", "WE", "WS") := tstrsplit(INFO, ";", fixed=T)]
  for(i in names(dt[,-c(1)])){
    dt[[i]] <- sub('.*=', '', dt[[i]])
  }
  return(dt)
}
edc <- splitInfoCol(edc)

#cpg data import
cpg <- fread("CpG_islands.bed", sep="\t", header=F, stringsAsFactors = F) # obtained from Genome Table Browser for hg19
names(cpg) <- c("chr", "start", "end", "cpginfo")
for (i in names(cpg[,c(4)])){
  cpg[[i]] <- sub('.*_', '', cpg[[i]])
  cpg <- cpg[, cpginfo:=as.numeric(cpginfo)]
}
cpg.edc <- cpg[chr=="chr1"]

#segdup data import
segdup <- fread("hg19_segdup.bed", sep="\t", header=F, stringsAsFactors = F) # Obtained from 1000 Genomes phase 3 release
names(segdup) <- c("chr", "start", "end")
segdup.edc <- segdup[chr=="chr1" & start >= 151973148 & start <= 153642014,]
segdup.edc[,c("start") := as.numeric(start)]
segdup.edc[,c("end") := as.numeric(end)]

#edc track:
edc.plot <- ggplot(data=edc) + geom_point(aes(x=as.numeric(edc$POS), y=as.numeric(edc$TC)), size=1, alpha=0.3) + theme_bw() + labs(y="Total reads per variant")+
  # window mean
  annotate("segment", x=151973148, xend=153642014, y=30960.75, yend=30960.75, alpha=0.5, color="red")+
  annotate("text", hjust=1, vjust=1, label="Post-QC\nVariants: 15636\nMean reads: 30961", x=153642014, y=50000, color="black") +
  ylim(0, 50000)+ xlim(151973148, 153642014)+
  # low quality and FLG annotations 
  annotate("rect", xmin=152555000, xmax=152590000, ymin=0, ymax=50000, color="red", fill="red", alpha=0.2)+
  annotate("text", hjust=0, label="LCE3C,\nLCE3B\n12571", x=152600000, y=50000, vjust=1, color="red") +
  
  annotate("rect", xmin=152235000, xmax=152245000, ymin=0, ymax=50000, color="red", fill="red", alpha=0.2)+
  annotate("text", hjust=1, label="18641", x=152235000, y=50000,vjust=1, color="red")+
  
  annotate("rect", xmin=152759000, xmax=152770000, ymin=0, ymax=50000, color="red", fill="red", alpha=0.2)+
  annotate("text", hjust=0, label="LCE1D\n22445", x=152772000, y=50000, vjust=1, color="red")+
  
  annotate("rect", xmin = 152274685, xmax = 152297517, ymin=0, ymax=50000, color="red", fill="red", alpha=0.2)+
  annotate("text", x=152300000, y=50000, label="FLG, FLG-AS1\n29058", hjust=0, vjust=1, color="red")

#segdup track:
segdup.edc <- segdup[chr=="chr1" & start >= 151973148 & start <= 153642014,]
segdup.edc[,c("start") := as.numeric(start)]
segdup.edc[,c("end") := as.numeric(end)]
segdup.plot <- ggplot(segdup.edc) + annotate("segment", x=segdup.edc[,start], xend=segdup.edc[,end], y=0, yend=0, size=40)+
  theme_bw() +  annotate("text", label="Segmental duplications", x=151973148, hjust=0, vjust=1, y=1) +theme(axis.ticks=element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_blank())+xlim(151973148, 153642014) + ylim(0,1)

#cpg islands track:
cpg.plot <- ggplot(cpg.edc) + annotate("segment", x=cpg.edc[,start], xend=cpg.edc[,end], y=0, yend=cpg.edc[,cpginfo], size=1) + theme_bw()+labs(y="Length")+theme(axis.line.y = element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y = element_blank())+ annotate("text", label="CpG Islands", x=151973148, hjust=0, vjust=1, y=150)

#graph it:
tracks(edc.plot, segdup.plot, cpg.plot, heights=c(5,1,1))

## EDC variant overlap graph-------
#setwd("/Users/claire/Documents/adrneh/EDC/")
illumina <- fread("ADRN_761_plus_EH_43_chr01_EDC_SNPs_multiVCF_vcftools_ingotest.filled.QC.recode.vcf", header=T, skip=1740L)
illumina <- illumina[FILTER=="PASS",] # FILTER == PASS only, DP > 7, GQ > 20
platypus <- fread("ADRN_EDC_799_platypus_maxvar8_minreads5.vcf_PASS.TR7.GQ20.NR7.recode.vcf", header=T, skip=51L) # FILTER == PASS only, TR > 7, GQ > 20, NR > 7

# subset illumina to individuals in common
illumina.IDs <- names(illumina)[10:813]
platypus.IDs <- names(platypus)[10:808]

intersect.IDs <- intersect(illumina.IDs, platypus.IDs) #799 total in common
illumina <- subset(illumina, select=c(names(illumina)[1:9], overlap.IDs))
platypus <- subset(platypus, select=c(names(platypus)[1:9], overlap.IDs)) #unchanged

variant.intersect <- intersect(illumina$POS, platypus$POS)
euler(list(Platypus=c(platypus$POS), Illumina=c(illumina$POS)))
#Platypus=3084, Illumina=3016, "Platypus&Illumina"=12552
draw.pairwise.venn(15568, 15636, 12522, fill="cornflowerblue", col="transparent", alpha = rep(0.5, 2), cex = rep(1.5, 3), cat.cex=rep(2, 2), category = c("Illumina", "Platypus"), cat.dist=rep(0.07, 2), margin=0.1)

## FLG concordance graph------
#setwd("/Users/claire/Documents/adrneh/EDC/ashg2017/")

concordance <- fread("concordance-longform.csv")

# Genotype counts for FLG variants by data source
ggplot(concordance, aes(source,count,fill=genotype))+
  geom_bar(position="stack",stat="identity")+
  facet_wrap(~variant, nrow=1)+
  scale_fill_viridis(discrete=T, direction=-1, name="Genotype")+
  theme_bw()
  ylim(0,766)+
  labs(x="Source", y="Count, individuals")
 
