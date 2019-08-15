## Author - Yue Zhang
## Date - Aug 15, 2019

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

# the config file, which has four columns, bam, total_reads, label and color
# bam     total_reads     label   color
# test.bam   60108409     BAP1_KO_input   black

conf_file <- "conf.txt"
genome <- "mm10"
outfile <- "EF_8a.pdf"
bp_flank <- 2500        # upstream and downstram flanking sequence
conf <- read.table(conf_file,sep="\t",header = T,stringsAsFactors = FALSE)

# require the number of cpu cores
ncpu <- nrow(conf)
registerDoMC(ncpu)

# read the gene model
genome.coord <- read.table("GRCm38.IGIS4.0.genes.cut.tsv",head=T,stringsAsFactors = F,sep = "\t")
colnames(genome.coord) = c("chrom","start","end","strand","gname")
# remove chr at chromosome name
chrom <- gsub("chr","",genome.coord$chrom)
genome.coord$chrom <- chrom

genome.coord <- genome.coord %>% filter(chrom != "MT")
genome.coord <- genome.coord %>% filter(! str_detect(chrom,'CHR|GL|JH'))

# sort genes
genome.coord <- genome.coord %>% arrange(chrom, start)
# remove redundant transcripts/genes with same TSS and TTS
genome.coord.unique <- genome.coord %>% group_by(gname) %>% filter(row_number() == 1)

# select chrom, start, end, strand, gname columns from gene model and creat GRange object
data <- genome.coord.unique %>% select(chrom,start,end,strand,gname)
gene <- GRanges(
    seqnames = data$chrom,
    ranges = IRanges(data$start, end = data$end, names = data$gname),
    strand = data$strand)


# parallel processing each bam files
d_all <- foreach(i = 1:nrow(conf), .combine=rbind) %dopar% {

   print(paste("I am processing the ",i, "th sample!",sep=""))
   # promoter
   proms <- GenomicRanges::promoters(gene, upstream=bp_flank, downstream=bp_flank)
   binnedSigs_up <- bamCount(conf[i,"bam"], proms, verbose=FALSE)
   d_up <- round(binnedSigs_up/conf[i,"total_reads"]*1000000,3)

   print(paste("I finished the processing of ",i, "th sample!",sep=""))
   return(d_up)

}

d_all <- t(d_all)

colnames(d_all) <- conf$label
rownames(d_all) <- data$gname

# filter genes by chipseq signal
d_all <- as.data.frame(d_all)
d <- d_all[ rowSums(d_all) >= 3.2, ]

gene_name <- rownames(d)

d <- log2(d+1)
d <- d %>% mutate(K4 = rowMeans(d[,c("WT1_H3K4me3","WT2_H3K4me3")]))
d <- d %>% mutate(K27 = rowMeans(d[,c("WT1_H3K27me3","WT2_H3K27me3")]))

# set seed and do K-mean clustering
set.seed(10)

# set center for four clusters
centers=rbind(c(6,6,1.5,1.5),c(2,2,6,6),c(5,5,4.5,4.5),c(1.2,1.2,1.5,1.5))
clusters <- kmeans(d[,c(1:4)], centers=centers,iter.max = 15)

# very similar results
#clusters <- kmeans(d[,c(1:4)], 4)

d$cluster <- as.factor(clusters$cluster)

#### function of adding the distance of each point to its center of each cluster
distance <- function(points1, points2) {
  dis <- 0
  for(i in 1:length(points1)) {
    dis <- dis + (points1[i]-points2[i])^2
  }
  return(sqrt(dis))
}

dis <- list()
for(i in 1:dim(d)[1]) {
   dis[i] <- distance(d[i,c(1:4)],clusters$centers[as.numeric(d[i,"cluster"]),])
}

d$dis <- unlist(dis)
rownames(d) <- gene_name

#### sort by distance to center of each cluster, generate heatmap

d_sort <- d %>% rownames_to_column('gene') %>% arrange(cluster, dis) %>% dplyr::mutate_if(is.numeric, round, 4) %>% column_to_rownames('gene')
write.table(d_sort, "count_in_RPM_sorted.txt",sep="\t",quote = FALSE, row.names = T)

cols <- hue_pal()(4)

pdf(outfile)
ha_row = rowAnnotation(df = data.frame(cluster= as.factor(as.character(d_sort$cluster))),
                       col = list( cluster = c("1" =  cols[1], "2" = cols[2], "3" = cols[3], "4" = cols[4])))
p <- Heatmap(d_sort[,1:4], name = "Log2RPM", col = colorRamp2(c(0,4,8), c("blue","white", "red")),cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = FALSE)
draw(p + ha_row)
dev.off()

## Display number of genes per cluster
d_sort$cluster %>% table()

