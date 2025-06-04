#!/usr/bin/env Rscript

library("argparse")
library("dplyr")

parser <- ArgumentParser(description= 'Split eQTL files to smaller files based on genes.')
parser$add_argument('--headerfile', '-f', help= 'Header file. Example: /scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/eqtl_Exons/header_Exons')
parser$add_argument('--summ', '-s', help= 'Summary stat file based on chromosome. Example: /scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/eqtl_Exons/chr1.summ.bed')
parser$add_argument('--outdir', '-o', help= 'Output dir.')

xargs <- parser$parse_args()


library(data.table)
header <- read.table(xargs$headerfile, comment.char = "")
c <- fread(xargs$summ, header = F)
header[1,1] <- "#snp_chr"
colnames(c) <- as.character(header[1,])
library(dplyr)
a <- c %>% group_split(GeneID)

for (i in 1:length(a)) {
  write.table(a[[i]],
              paste0(xargs$outdir, "/",
                     a[[i]]$GeneID[1], "__InsPIRE_Islets_Exons__", gsub("chr", "", a[[i]][1, 1]), ":", a[[i]]$TSS[1], ".bed"),
              col.names = T, row.names = F, sep = "\t", quote = F)
}

