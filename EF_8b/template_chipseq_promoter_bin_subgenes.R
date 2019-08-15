## Author - Yue Zhang
## Date Aug 15, 2019
## Genentech, A Member of the Roche Group
library(GenomicRanges)
library(Rsamtools)
library(bamsignals)
library(dplyr)
library(foreach)
library(doMC)
library(reshape2)
library(ggplot2)
library(caTools)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(scales)
library(magrittr)
#library(xlsx)

# load list of genes that have Irf6 motif under an ATACseq peak near them
sub_gene <- read.table("ISRE_IRF_motif_fimo_intersect_ATACseq_intersect_genes.txt",sep="\t",stringsAsFactors = FALSE, header=FALSE)
sub_gene <- sub_gene[,c('V4'),drop=FALSE]
colnames(sub_gene) <-c("gname")
sub_gene <- sub_gene %>% group_by(gname) %>% filter(row_number() == 1)  # 930 -> only 869 present in the cluster



# colors for four clusters
cols <- hue_pal()(4)

# argument
conf_file <- "conf_all.txt"
genome <- "mm10"
outfile <- "EF_8b.pdf"
outfile_profile <- "EF_8b_profile.pdf"
bp_flank <- 2500           # upstream and downstram flanking sequence
binsize <- 50              # bin size = 50 bp

filterBycluster <- FALSE
SelectedCluster <- "3"

conf <- read.table(conf_file,sep="\t",header = T,stringsAsFactors = FALSE)

# require the number of cpu cores
ncpu <- nrow(conf)
registerDoMC(ncpu)

########## load genes annotation and sorted by cluter number
genome.coord <- read.table("GRCm38.IGIS4.0.genes.cut.tsv",head=T,stringsAsFactors = F,sep = "\t")
colnames(genome.coord) = c("chrom","start","end","strand","gname")
# remove chr at chromosome name
chrom <- gsub("chr","",genome.coord$chrom)
genome.coord$chrom <- chrom
genome.coord <- genome.coord %>% filter(chrom != "M")

# sort genes
genome.coord <- genome.coord %>% arrange(chrom, start)
# remove redundant transcripts/genes with same TSS and TTS
genome.coord.unique <- genome.coord %>% group_by(gname) %>% filter(row_number() == 1)

### sort gene by order in cluster and filter no expressed gene
gene.order <- read.table("count_in_RPM_sorted.txt",head=T,stringsAsFactors = F)
gene.order$gname <- rownames(gene.order)

if(filterBycluster) {
    gene.order <- gene.order %>% filter(cluster == SelectedCluster)
}

genome.coord.unique.sort <- left_join(gene.order,genome.coord.unique,by ="gname")
genome.coord.final <- genome.coord.unique.sort[,c(colnames(genome.coord.unique),"cluster")]


########## filter using subgenes, finally 815 genes are kept for analysis

genome.coord.final <- semi_join(genome.coord.final,sub_gene,by ="gname")

# select chrom, start, end, strand, gname columns from data base
# make a GRange object
gene <- GRanges(
    seqnames = genome.coord.final$chrom,
    ranges = IRanges(genome.coord.final$start, end = genome.coord.final$end, names = genome.coord.final$gname),
    strand = genome.coord.final$strand)


# parallel processing each bam files
d_all <- foreach(i = 1:nrow(conf), .combine = cbind) %dopar% {

   print(paste("I am processing the ",i, "th sample!",sep=""))
   # promoter
   proms <- GenomicRanges::promoters(gene, upstream=bp_flank, downstream=bp_flank)
   binnedSigs_up <- bamProfile(conf[i,"bam"], proms, binsize=binsize,verbose=FALSE)
   d_up <- t(alignSignals(binnedSigs_up))

   d_up <- d_up/conf[i,"total_reads"]*1000000
   d_up <- round(log2(d_up+1),3)

   return(d_up)
   print(paste("I finished the processing of ",i, "th sample!",sep=""))

}

colnames(d_all) <- 1:dim(d_all)[2]
rownames(d_all) <- genome.coord.final$gname

## Add RNAseq data as row annotation
rna = read.table('lfc_NGS2488_all.txt',sep="\t",stringsAsFactors = FALSE)
mat = data.frame(symbol = rownames(d_all), stringsAsFactors = FALSE)
mat = left_join(mat,rna,by='symbol')
mat = mat[,c('symbol','lfc','p')]
mat = data.frame(mat, row.names = 'symbol',stringsAsFactors = FALSE)


## To mark only significant up/down genes
lfc2 = rep(0, nrow(mat))
lfc2[which(mat$lfc > 1 & mat$p < 0.05)] = 1
lfc2[which(mat$lfc < -1 & mat$p < 0.05)] = -1

## To highlight only 6 genes that were validated
lfc3 = rep(0, nrow(mat))
lfc3[which(rownames(mat) %in% c('Rora','Grhl3','Ocln','Cers3','Sptlc3','Pnpla1'))] = 1

### heatmap
pdf(outfile)
ha_column = HeatmapAnnotation(df = data.frame(treatment = c(rep("WT1_H3K4me3", 100), rep("WT1_H3K27me3", 100),rep("WT1_H3K27Ac", 100), rep("WT1_ATAC", 100))),
                              col = list(treatment = c("WT1_H3K4me3" =  "red", "WT1_H3K27me3" = "darkred","WT1_H3K27Ac" =  "brown", "WT1_ATAC" = "orange")))

ha_row = rowAnnotation(df = data.frame(cluster= as.factor(as.character(genome.coord.final$cluster))),
                              col = list( cluster = c("1" =  cols[1], "2" = cols[2], "3" = cols[3], "4" = cols[4])))

## Uncomment below if you only want to plot significantly up/down genes
ha_row2 = rowAnnotation(df = data.frame("lfc"=lfc2), col=list("lfc"= colorRamp2(c(-1,0,1), c("blue","white", "red"))))


p1 <- Heatmap(d_all, bottom_annotation = ha_column, name = "log2RPM", col = colorRamp2(c(0,3), c("white", "red")),cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = F,row_names_gp = gpar(fontsize = 7),show_column_names = FALSE)
draw(p1 + ha_row + ha_row2)
dev.off()
