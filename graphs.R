# Claire Malley
# Graphs for ASHG2017 poster
# Compiled together on September 14

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

# EDC coverage graph -----
#edc data import
#setwd("/Users/claire/Documents/adrneh/EDC")
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
#setwd("/Users/claire/Documents/adrneh/EDC")
cpg <- fread("CpG_islands.bed", sep="\t", header=F, stringsAsFactors = F) # obtained from Genome Table Browser for hg19
names(cpg) <- c("chr", "start", "end", "cpginfo")
for (i in names(cpg[,c(4)])){
  cpg[[i]] <- sub('.*_', '', cpg[[i]])
  cpg <- cpg[, cpginfo:=as.numeric(cpginfo)]
}
cpg.edc <- cpg[chr=="chr1"]

#segdup data import
#setwd("/Users/claire/Documents/adrneh/Scripts_ADRN_from_start/Files_general")
segdup <- fread("hg19_segdup.bed", sep="\t", header=F, stringsAsFactors = F) # Obtained from 1000 Genomes phase 3 release
names(segdup) <- c("chr", "start", "end")
segdup.edc <- segdup[chr=="chr1" & start >= 151973148 & start <= 153642014,]
segdup.edc[,c("start") := as.numeric(start)]
segdup.edc[,c("end") := as.numeric(end)]

# structural variants data import
#setwd("/Users/claire/Documents/adrneh/EDC/ashg2017")
structvar <- fread("CNVindels-hg19Tables.txt", header=T, stringsAsFactors = F) # Obtained from UCSC table browser schema "DGV Struct Var / hg19.dgvMerged"
# edc structural variants track
structvar.plot <- ggplot(data=structvar) + geom_segment(aes(x=structvar[,chromStart], xend=structvar[,chromEnd], y=0, yend=0, color=factor(structvar[,varType])), size=50, alpha=0.15)+
  theme_bw() +  
  #annotate("text", label="Structural variants", x=151973148, hjust=0, vjust=1, y=1) +
  theme(text = element_text(size=15), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none", plot.margin = unit(c(-1.5,0,0.5,0), "lines"))+
  coord_cartesian(xlim=c(151973148, 153642014)) + ylim(0,1)+ scale_color_viridis(discrete = T)

#edc track:
edc.plot <- ggplot(data=edc) + geom_point(aes(x=as.numeric(edc$POS), y=as.numeric(edc$TC)), size=1, alpha=0.3) + theme_bw() + labs(y="Total reads per variant")+
  # window mean
  annotate(size=4, "text", hjust=1, vjust=1, label="Variants: 15636\nMean reads: 30961", x=153642014, y=50000, color="black") +
  ylim(0, 52000)+
  xlim(151973148, 153642014)+
  annotate("text", size=4, vjust=1, hjust=0, label="chr1:151973148-153642014", x=151973148, y=50000)+
  # low quality and FLG annotations
  #annotate(size=5, "text", label="Lower coverage spots", x=152235000, y=52000,vjust=1, hjust=0, color="skyblue4")+
 # annotate("rect", xmin=152555000, xmax=152590000, ymin=0, ymax=50000, color="skyblue4", fill="skyblue4", alpha=0.2)+
  #annotate("rect", xmin=152235000, xmax=152245000, ymin=0, ymax=50000, color="skyblue4", fill="skyblue4", alpha=0.2)+
  #annotate("rect", xmin=152759000, xmax=152770000, ymin=0, ymax=50000, color="skyblue4", fill="skyblue4", alpha=0.2)+
 # annotate("rect", xmin = 152274685, xmax = 152297517, ymin=0, ymax=50000, color="skyblue4", fill="skyblue4", alpha=0.2)+
  theme(text = element_text(size=15), legend.position="none", axis.title.x=element_blank(), axis.text.x = element_text(margin=margin(0,0,0,0,"pt"), hjust=0), axis.title.y=element_blank())

#segdup track:
# segdup.edc <- segdup[chr=="chr1" & start >= 151973148 & start <= 153642014,]
# segdup.edc[,c("start") := as.numeric(start)]
# segdup.edc[,c("end") := as.numeric(end)]
# segdup.plot <- ggplot(segdup.edc) + annotate("segment", x=segdup.edc[,start], xend=segdup.edc[,end], y=0, yend=0, size=50)+
#   theme_bw() +  
#   #annotate("text", label="Segmental duplications", x=151973148, hjust=0, vjust=1, y=1) +
#   theme(axis.ticks=element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
#   xlim(151973148, 153642014) + ylim(0,1)

#cpg islands track:
# cpg.plot <- ggplot(cpg.edc) + annotate("segment", x=cpg.edc[,start], xend=cpg.edc[,end], y=0, yend=0, size=50) + theme_bw()+
#   theme(axis.line.y = element_blank(), panel.grid.major.y=element_blank(), axis.title.y=element_blank(), text = element_text(size=15), axis.text.x = element_text(margin=margin(0,0,0,0,"pt"), hjust=0))+
#   annotate("text", label="CpG Islands", x=151973148, hjust=0, vjust=1, y=150)
# #trbl
#graph it:
#cowplot::plot_grid(edc.plot, structvar.plot, segdup.plot, rel_heights=c(6,2,1), ncol = 1, align = "v")

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
#setwd("/Users/claire/Documents/adrneh/EDC/")
stat6 <- fread("ADRN_stat6_799_platypus_maxvar8_minreads5.QC.vcf", header=T, skip=51L, stringsAsFactors = F)
stat6 <- splitInfoCol(stat6)
stat6.plot <- ggplot(data=stat6) + geom_point(aes(x=as.numeric(stat6$POS), y=as.numeric(stat6$TC)), size=1, alpha=0.3) + theme_bw() +
  # window mean
  annotate(size=4, "text", hjust=1, vjust=1, label="Variants: 12736\nMean reads: 28494", x=58331753, y=50000, color="black") +
  annotate("text", size=4, vjust=1, hjust=0, label="chr12:56662638-58331753", x=56662638, y=50000)+
  ylim(0, 50000)+xlim(56662638, 58331753)+
  theme(text = element_text(size=15), legend.position="none", axis.title.x=element_blank(), axis.text.x = element_text(margin=margin(0,0,0,0,"pt"), hjust=0), axis.title.y=element_blank())

# IFNG -----
#setwd("/Users/claire/Documents/adrneh/EDC/")
ifng <- fread("ADRN_ifng_799_platypus_maxvar8_minreads5.QC.vcf", header=T, skip=51L, stringsAsFactors = F)
ifng <- splitInfoCol(ifng)
ifng.plot <- ggplot(data=ifng) + geom_point(aes(x=as.numeric(ifng$POS), y=as.numeric(ifng$TC)), size=1, alpha=0.3) + theme_bw() +
  # window mean
  annotate(size=4, "text", hjust=1, vjust=1, label="Variants: 16761\nMean reads: 32307", x=69385582, y=50000, color="black") +
  annotate("text", size=4, vjust=1, hjust=0, label="chr16:67716475-69385582", x=67716475, y=50000)+
  ylim(0, 50000)+ xlim(67716475, 69385582)+
  theme(text = element_text(size=15), legend.position="none", axis.title.x=element_blank(), axis.text.x = element_text(margin=margin(0,0,0,0,"pt"), hjust=0), axis.title.y=element_blank())

# graph IL4R, IFNG, and STAT6 together---
#grid.arrange(stat6.plot, ifng.plot, il4r.plot, nrow=3, ncol=1)

#cowplot::plot_grid(il4r.plot, ifng.plot, stat6.plot, rel_heights=c(1,1,1), ncol = 1, align = "v")

cowplot::plot_grid(edc.plot, structvar.plot, il4r.plot, structvar.il4r.plot, ifng.plot, structvar.ifng.plot, stat6.plot,  structvar.stat6.plot, rel_heights=c(4,1,4,1,4,1,4,1), ncol = 1, align = "v")

#

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
#nrow(illumina.preQC.miss[!(illumina.preQC.miss$POS %in% platypus.exclusive.snps$POS)])


variant.indel.intersect <- intersect(illumina.indel$POS, platypus.indel$POS)
euler(list(Platypus=c(platypus.indel$POS), Illumina=c(illumina.indel$POS)))
#Platypus=1478, Illumina=8, "Platypus&Illumina"=60
draw.pairwise.venn(68, 1538, 60, fill="cornflowerblue", col="transparent", alpha = rep(0.5, 2), cex = rep(1.5, 3), cat.cex=rep(2, 2), category = c("Illumina", "Platypus"), cat.dist=rep(0.07, 2), margin=0.1)
if(!is.null(dev.list())) dev.off()

#

# FLG concordance graphs------
#setwd("/Users/claire/Documents/adrneh/EDC/ashg2017/")

concordance <- fread("concordance-longform.csv")

concordance <- concordance[source=="Platypus" | source=="Illumina",]

#concordance <- fread("concordance-updated.csv")


# Genotype counts for FLG variants by data source - organized by variant------
ggplot(concordance, aes(source,count.updated2,fill=genotype))+
  geom_bar(position="stack",stat="identity")+
  facet_wrap(~variant, nrow=1)+
  scale_fill_viridis(discrete=T, direction=-1, name="Genotype")+
  theme_bw()+
  ylim(0,766)+
  #labs(x="Source", y="Count, individuals")+
  theme(text = element_text(size=18), axis.title = element_blank())
  
# Genotype counts for FLG variants by data source - organized instead by source-------
g_legend<-function(a.gplot){
  #https://stackoverflow.com/a/12539820
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

conc.platypus <- concordance[source=="Platypus",]
platypus.flg <- ggplot(conc.platypus, aes(variant,count.updated2,fill=genotype))+
  geom_bar(position="stack",stat="identity")+
  scale_fill_viridis(discrete=T, direction=-1, name="Genotype")+
  theme_bw()+
  ylim(0,766)+
  labs(title="Platypus multi-sample VCF", y="Count, individuals")+
  theme(text = element_text(size=20), legend.position="none", axis.title.y = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle=45, hjust=1, margin=margin(-20,-20,0,-20,"pt")), axis.ticks.y = element_blank(), axis.ticks.x=element_line(color="white"),axis.line.x=element_line(color="white"), panel.border = element_blank())

conc.illumina <- concordance[source=="Illumina",]
illumina.flg <- ggplot(conc.illumina, aes(variant,count.updated2,fill=genotype, label=count.updated))+
  geom_bar(position="stack",stat="identity")+
  scale_fill_viridis(discrete=T, direction=-1, name="Genotype")+
  #scale_fill_brewer( type = "div" , palette = "Spectral", direction=1 )+
theme_bw()+
  ylim(0,766)+
  labs(title="Illumina Isaac single-sample VCFs")+
  theme(text = element_text(size=20), legend.position="none", axis.title.y = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle=45, hjust=1, margin=margin(-20,-20,0,-20,"pt")), axis.ticks.y = element_blank(), axis.ticks.x=element_line(color="white"),axis.line.x=element_line(color="white"), panel.border = element_blank())

#+  geom_text(size = 5, color="white", position = position_stack(vjust = 0.5))
# 
# conc.taqman <- concordance[source=="TaqMan",]
# taqman.flg <- ggplot(conc.taqman, aes(variant,count,fill=genotype))+
#   geom_bar(position="stack",stat="identity")+
#   scale_fill_viridis(discrete=T, direction=-1, name="Genotype")+
#   theme_hc()+
#   ylim(0,766)+
#   labs(title="TaqMan")+
#   theme(text = element_text(size=20), axis.title.y = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle=45, hjust=1, margin=margin(-20,-20,0,20,"pt")), axis.ticks = element_blank(), legend.position="none")

empty.plot <- ggplot()+theme_hc() #takes up second row, first column space in order to center the legend
#legend <- g_legend(platypus.flg)

#grid.arrange(platypus.flg, illumina.flg, taqman.flg, empty.plot, legend, ncol=3, nrow=2, heights=c(7/8, 1/8))
grid.arrange(platypus.flg, illumina.flg, empty.plot, ncol=2, nrow=2, heights=c(7/8, 1/8))

# EDC concordance, Illumina vs. Platypus for 6 individuals ----
#setwd("/Users/claire/Documents/adrneh/EDC/ashg2017/")
matches.data <- fread("matches-data.csv")
matches.plot <-ggplot(data=matches.data, aes(x=matches.data$`Type of match or mismatch`, y=matches.data$`Mean total`))+
  geom_bar(stat="identity")+
  theme_hc()+
  labs(x="Type of match or mismatch", y="Mean total")+
  theme(text = element_text(size=20), axis.text.x = element_text(angle=45, hjust=1, margin=margin(-20,-20,0,20,"pt")), axis.ticks = element_blank())+
  geom_text(x=2, y=8500, label="Mismatches/Matches = 4.8%", size=6, hjust=0)+
  geom_segment(x=2, xend=2.5, y=1000, yend=8000)+
  geom_segment(x=1.6, xend=2.5, y=6000, yend=8000)

matches.plot
