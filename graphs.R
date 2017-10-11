# Claire Malley
# Graphs for ASHG2017 poster
# Compiled together on October 11

library(data.table)
library(viridis)
library(ggplot2)
library(venneuler)
library(VennDiagram)
library(eulerr)
options(scipen = 999)
library(gridExtra)
library(cowplot)
library(stringi)

# Discussion about the R code to create the main figure can be found at https://cemalley.com/2017/10/05/multi-plotting-biological-data-just-one-solution-with-r-ggplot-cowplot/
# As of Oct 11

# EDC coverage graph -----
#edc data import
setwd("/Users/claire/Documents/adrneh/EDC")
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

# structural variants data import
setwd("/Users/claire/Documents/adrneh/EDC/ashg2017")

#edc track:
edc.plot <- ggplot(data=edc) + geom_point(aes(x=as.numeric(edc$POS), y=as.numeric(edc$TC)), size=1, alpha=0.3) + theme_bw() + labs(y="Total reads per variant")+
  # window mean
  annotate(size=4, "text", hjust=1, vjust=1, label="Variants: 15636\nMean reads: 30961", x=153642014, y=50000, color="black") +
  ylim(0, 52000)+
  xlim(151973148, 153642014)+
  annotate("text", size=4, vjust=1, hjust=0, label="chr1:151973148-153642014", x=151973148, y=50000)+
  theme(text = element_text(size=15), legend.position="none", axis.title.x=element_blank(), axis.text.x = element_text(margin=margin(0,0,0,0,"pt"), hjust=0), axis.title.y=element_blank())

# IL4R coverage graph----
il4r <- fread("ADRN_il4ra_799_platypus_maxvar8_minreads5.QC.vcf", header=T, skip=51L, stringsAsFactors = F)
il4r <- splitInfoCol(il4r)
il4r.plot <- ggplot(data=il4r) + geom_point(aes(x=as.numeric(il4r$POS), y=as.numeric(il4r$TC)), size=1, alpha=0.3 ) + theme_bw()+
  #window mean
  annotate(size=4, "text", hjust=1, vjust=1, label="Variants: 17428\nMean reads: 29990", x=28185184, y=50000, color="black")+
  annotate("text", size=4, vjust=1, hjust=0, label="chr16:26516211-28185184", x=26516211, y=50000)+
  ylim(0, 50000)+
  xlim(26516211, 28185184)+
  theme(text = element_text(size=15), legend.position="none", axis.title.x=element_blank(), axis.text.x = element_text(margin=margin(0,0,0,0,"pt"), hjust=0), axis.title.y=element_blank())

# STAT6 ----
setwd("/Users/claire/Documents/adrneh/EDC/")
stat6 <- fread("ADRN_stat6_799_platypus_maxvar8_minreads5.QC.vcf", header=T, skip=51L, stringsAsFactors = F)
stat6 <- splitInfoCol(stat6)
stat6.plot <- ggplot(data=stat6) + geom_point(aes(x=as.numeric(stat6$POS), y=as.numeric(stat6$TC)), size=1, alpha=0.3) + theme_bw() +
  # window mean
  annotate(size=4, "text", hjust=1, vjust=1, label="Variants: 12736\nMean reads: 28494", x=58331753, y=50000, color="black") +
  annotate("text", size=4, vjust=1, hjust=0, label="chr12:56662638-58331753", x=56662638, y=50000)+
  ylim(0, 50000)+xlim(56662638, 58331753)+
  theme(text = element_text(size=15), legend.position="none", axis.title.x=element_blank(), axis.text.x = element_text(margin=margin(0,0,0,0,"pt"), hjust=0), axis.title.y=element_blank())

# IFNG -----
setwd("/Users/claire/Documents/adrneh/EDC/")
ifng <- fread("ADRN_ifng_799_platypus_maxvar8_minreads5.QC.vcf", header=T, skip=51L, stringsAsFactors = F)
ifng <- splitInfoCol(ifng)
ifng.plot <- ggplot(data=ifng) + geom_point(aes(x=as.numeric(ifng$POS), y=as.numeric(ifng$TC)), size=1, alpha=0.3) + theme_bw() +
  # window mean
  annotate(size=4, "text", hjust=1, vjust=1, label="Variants: 16761\nMean reads: 32307", x=69385582, y=50000, color="black") +
  annotate("text", size=4, vjust=1, hjust=0, label="chr16:67716475-69385582", x=67716475, y=50000)+
  ylim(0, 50000)+ xlim(67716475, 69385582)+
  theme(text = element_text(size=15), legend.position="none", axis.title.x=element_blank(), axis.text.x = element_text(margin=margin(0,0,0,0,"pt"), hjust=0), axis.title.y=element_blank())


# structural variants tracks
setwd("/Users/claire/Documents/adrneh/EDC/ashg2017/")
structvar <- fread("CNVindels-hg19Tables.txt", header=T, stringsAsFactors = F) # Obtained from UCSC table browser schema "DGV Struct Var / hg19.dgvMerged"
# edc structural variants track
structvar.plot <- ggplot(data=structvar) + geom_segment(aes(x=structvar[,chromStart], xend=structvar[,chromEnd], y=0, yend=0, color=factor(structvar[,varType])), size=50, alpha=0.15)+
  theme_bw() +  
  theme(text = element_text(size=15), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none", plot.margin = unit(c(-1.5,0,0.5,0), "lines"))+
  coord_cartesian(xlim=c(151973148, 153642014)) + ylim(0,1)+ scale_color_viridis(discrete = T)

structvar.il4r <- fread("IL4R-CNVindels-hg19Tables.txt", header=T, stringsAsFactors = F) # Obtained from UCSC table browser schema "DGV Struct Var / hg19.dgvMerged"
# structural variants track
structvar.il4r.plot <- ggplot(data=structvar.il4r) + geom_segment(aes(x=structvar.il4r[,chromStart], xend=structvar.il4r[,chromEnd], y=0, yend=0, color=factor(structvar.il4r[,varType])), size=50, alpha=0.15)+
  theme_bw() +
  theme(text = element_text(size=15), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none", plot.margin = unit(c(-1.5,0,0.5,0), "lines"))+
  coord_cartesian(xlim=c(26516211, 28185184)) + ylim(0,1)+ scale_color_viridis(discrete = T)

structvar.ifng <- fread("IFNG-CNVindels-hg19Tables.txt", header=T, stringsAsFactors = F)
structvar.ifng.plot <- ggplot(data=structvar.ifng) + geom_segment(aes(x=structvar.ifng[,chromStart], xend=structvar.ifng[,chromEnd], y=0, yend=0, color=factor(structvar.ifng[,varType])), size=50, alpha=0.15)+
  theme_bw() +
  theme(text = element_text(size=15), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none", plot.margin = unit(c(-1.5,0,0.5,0), "lines"))+
  coord_cartesian(xlim=c(67716475, 69385582)) + ylim(0,1)+ scale_color_viridis(discrete = T)

structvar.stat6 <- fread("STAT6-CNVindels-hg19Tables.txt", header=T, stringsAsFactors = F)
structvar.stat6.plot <- ggplot(data=structvar.stat6) + geom_segment(aes(x=structvar.stat6[,chromStart], xend=structvar.stat6[,chromEnd], y=0, yend=0, color=factor(structvar.stat6[,varType])), size=50, alpha=0.15)+
  theme_bw() +
  theme(text = element_text(size=15), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none", plot.margin = unit(c(-1.5,0,0.5,0), "lines"))+
  coord_cartesian(xlim=c(56662638, 58331753)) + ylim(0,1)+ scale_color_viridis(discrete = T)

# graph IL4R, IFNG, and STAT6 together---

cowplot::plot_grid(edc.plot, structvar.plot, il4r.plot, structvar.il4r.plot, ifng.plot, structvar.ifng.plot, stat6.plot,  structvar.stat6.plot, rel_heights=c(4,1,4,1,4,1,4,1), ncol = 1, align = "v")

# EDC variant overlap graph (venn diagram)-------
setwd("/Users/claire/Documents/adrneh/EDC/")
illumina <- fread("ADRN_761_plus_EH_43_chr01_EDC_SNPs_multiVCF_vcftools_ingotest.filled.QC.recode.vcf", header=T, skip=1740L)
illumina <- illumina[FILTER=="PASS",] # FILTER == PASS only, DP > 7, GQ > 20
platypus <- fread("ADRN_EDC_799_platypus_maxvar8_minreads5.vcf_PASS.TR7.GQ20.NR7.recode.vcf", header=T, skip=51L) # FILTER == PASS only, TR > 7, GQ > 20, NR > 7

# subset illumina to individuals in common
illumina.IDs <- names(illumina)[10:813]
platypus.IDs <- names(platypus)[10:808]

intersect.IDs <- intersect(illumina.IDs, platypus.IDs) #799 total in common
illumina <- subset(illumina, select=c(names(illumina)[1:9], intersect.IDs))
platypus <- subset(platypus, select=c(names(platypus)[1:9], intersect.IDs)) #unchanged

# separate by SNPs and indels
illumina[, Type:= ifelse( ((str_length(REF) > 1)| (str_length(ALT) > 1)), "INDEL", "SNP" )]
platypus[, Type:= ifelse( ((str_length(REF) > 1)| (str_length(ALT) > 1)), "INDEL", "SNP" )]

illumina.snp <- illumina[Type=="SNP",] #15500 snps
platypus.snp <- platypus[Type=="SNP",] #14098 snps

illumina.indel <- illumina[Type=="INDEL",] #68 indels
platypus.indel <- platypus[Type=="INDEL",] #1538 indels

variant.snp.intersect <- intersect(illumina.snp$POS, platypus.snp$POS)
euler <- euler(list(Platypus=c(platypus.snp$POS), Illumina=c(illumina.snp$POS)))
#Platypus=1900, Illumina=3302, "Platypus&Illumina"=12198
draw.pairwise.venn(15500, 14098, 12198, fill="#A9A9A9", col="transparent", alpha = rep(0.4, 2), cex = rep(3, 3), cat.cex=rep(3, 2), category = c("Illumina", "Platypus"), cat.dist=rep(0.06, 2), cat.pos=c(-45, 45), margin=0.1, fontfamily = "sans", cat.fontfamily="sans")
if(!is.null(dev.list())) dev.off()

#find snps not in platypus and whether they overlap with indels in platypus
platypus[,LengthRef:=nchar(REF)]
platypus.alts <- platypus[,c("POS", "ALT")]
platypus.alts[,c("A", "B", "C", "D", "E", "F", "G", "H"):= tstrsplit(ALT, ',', names=FALSE)]
mysubset <- c("A", "B", "C", "D", "E", "F", "G", "H")
platypus.alts[, (mysubset) := lapply(.SD, stri_length), .SDcols = mysubset]
platypus.alts[, LengthAlt := apply(.SD, 1, max, na.rm=T), .SDcols = mysubset ]

platypus <- merge(platypus, platypus.alts[,c("POS", "LengthAlt")], by="POS", all.x=T, all.y=T)
platypus[,MaxIndelLength := apply(.SD, 1, max, na.rm=T), .SDcols= c("LengthRef", "LengthAlt")]
platypus[,BPAffected := MaxIndelLength + POS]

illumina.exclusive.snps <- illumina.snp[!(illumina.snp$POS %in% platypus.snp$POS)]
illumina.exclusive.snps <- merge(illumina.exclusive.snps, platypus[,c("POS", "Type")], by="POS", all.x=T, all.y=F)
nrow(illumina.exclusive.snps[Type.y=="INDEL",]) #292 snps according to illumina are definitely indels according to platypus
nrow(illumina.exclusive.snps[Type.y=="SNP",]) #0

illumina.exclusive.snps.pos <- c(illumina.exclusive.snps$POS, use.names=F)

platypus.subset <- platypus[,c("POS", "BPAffected")]
illumina.subset <- illumina.exclusive.snps[,c("POS")]
illumina.subset[,end := POS]
names(illumina.subset)[1] <- "start"
names(platypus.subset) <- c("start", "end")
setkey(illumina.subset, start, end)
setkey(platypus.subset, start, end)
overlaps.table <- foverlaps(platypus.subset, illumina.subset, type="any", nomatch=0L) #573 "any" kind of overlaps
length(unique(overlaps.table$i.start)) #428
length(unique(overlaps.table$i.end)) #424

# next need to find what happened for the remaining 2729 illumina-exclusive snps
#platypus[((platypus$POS %in% overlaps.table$start.i)), c("POS", "REF", "ALT")]

remaining.3000 <- illumina.subset[!(illumina.subset$start %in% overlaps.table$i.start),]

#read in earlier qc step platypus data that is not PASS only
platypus.preQC <- fread("ADRN.platypus.EDC.filtered.recode.vcf", skip=48L, header=T)

remaining.2446 <- (platypus.preQC[POS %in% remaining.3000$start,]) #2446

remaining.352 <- illumina.subset[(!(illumina.subset$start %in% remaining.2446$POS)) & (!(illumina.subset$start %in% overlaps.table$start))] #352


#find snps not in illumina
platypus.exclusive.snps <- platypus.snp[!(platypus.snp$POS %in% illumina.snp$POS)] #1900

illumina.preQC <- fread("/Users/claire/Documents/adrneh/EDC/ashg2017/ADRN_chr1.EDC.Illumina.not.backfilled.vcf", header=T, skip=52L)
illumina.preQC.miss <- illumina.preQC[(illumina.preQC$POS %in% platypus.exclusive.snps$POS)] #1750

variant.indel.intersect <- intersect(illumina.indel$POS, platypus.indel$POS)
euler(list(Platypus=c(platypus.indel$POS), Illumina=c(illumina.indel$POS)))
#Platypus=1478, Illumina=8, "Platypus&Illumina"=60
draw.pairwise.venn(68, 1538, 60, fill="cornflowerblue", col="transparent", alpha = rep(0.5, 2), cex = rep(1.5, 3), cat.cex=rep(2, 2), category = c("Illumina", "Platypus"), cat.dist=rep(0.07, 2), margin=0.1)
if(!is.null(dev.list())) dev.off()

#

# FLG concordance graph------
setwd("/Users/claire/Documents/adrneh/EDC/ashg2017/")

concordance <- fread("concordance-longform.csv")
concordance <- concordance[source=="Platypus" | source=="Illumina",]

# Genotype counts for FLG variants by data source - organized by variant
ggplot(concordance, aes(source,count.updated2,fill=genotype))+
  geom_bar(position="stack",stat="identity")+
  facet_wrap(~variant, nrow=1)+
  scale_fill_viridis(discrete=T, direction=-1, name="Genotype")+
  theme_bw()+
  ylim(0,766)+
  theme(text = element_text(size=18), axis.title = element_blank())


