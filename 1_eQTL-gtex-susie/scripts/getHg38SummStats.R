#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tidyr")
library("data.table")

parser <- ArgumentParser(description= 'Adapt GTEx summary stats to custom reference genotype data.\n
                                       Some part of this script is hard-coding col names, so check carefully.')
parser$add_argument('--ukbb_hg38', '-u', help= 'VCF file on hg38 containing SNPs. Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/ukbb-rsid/ukbb-rsid.vcf.gz')
parser$add_argument('--summStat', '-t', help= 'Text file containing GTEx summary stats. Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/3_t1d-eQTL_GTEx-coloc/results/gtex_indexed/chr1.bed"')
parser$add_argument('--geneInfo', '-g', help= 'Text file containing gene name and coor in hg38.\n
                                              Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/0_data/Pancreas.v8.egenes.txt.gz"')
parser$add_argument('--outdir', '-o', help= 'Output dir.')

xargs <- parser$parse_args()


summStat <- fread(xargs$summStat, header = T)
#gene_id:                  GENCODE/Ensembl gene ID
#variant_id:               variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b38
#tss_distance:             distance between variant and transcription start site. Positive when variant is downstream of the TSS, negative otherwise
#ma_samples:               number of samples carrying the minor allele
#ma_count:                 total number of minor alleles across individuals
#pval_nominal:             nominal p-value
#slope:                    regression slope
#slope_se:                 standard error of the regression slope
#summStat[, c("snp_chr", "snp_stop", "ref_gtex", "alt_gtex", "gtex_code") := tstrsplit(variant_id, "_")]
#summStat <- as.data.frame(summStat)
#summStat$snp_chr <- as.character(summStat$snp_chr)
#summStat$snp_stop <- as.numeric(summStat$snp_stop)
#summStat$snp_start <- summStat$snp_stop-1

print("start merging summ stat file and vcf file")

ukbb_hg38 <- fread(xargs$ukbb_hg38, skip = "CHROM")
test <- inner_join(as.data.frame(ukbb_hg38[,1:5]), summStat, by = c("#CHROM" = "snp_chr", "POS" = "snp_stop")) #1:5 here means we will exclude QUAL, FILTER and INFO columns of ukbb file
#needed info - something like the following, depending on what GWAS file gave:
##snp_chrom	snp_start	snp_end	SNP	REF	ALT	minor_allele	maf	n_complete_samples	AC_hg19	ytx_hg19
##beta	se	tstat	pval	SNP_hg19	multiply

### check matching between ref/alt alleles in different ref
test$multiply <- ifelse(test$REF == test$alt_gtex & test$ALT == test$ref_gtex, -1,
                         ifelse(test$ALT == test$alt_gtex & test$REF == test$ref_gtex, 1, NA))
test <- test[!is.na(test$multiply),] #on lines where none of the cases in the above condition holds, drop these lines
test$slope <- test$slope * test$multiply
test <- test[, c("#CHROM", "snp_start", "POS", "ID", "REF", "ALT", "ref_gtex", "alt_gtex",
		 "gene_id", "gtex_code", "tss_distance", "ma_samples", "ma_count", "maf",
		 "pval_nominal", "slope", "slope_se", "multiply")]
colnames(test) <- c("#snp_chrom", "snp_start", "snp_end", "SNP", "REF", "ALT", "ref_gtex", "alt_gtex",
		    "gene_id", "gtex_code", "tss_distance", "ma_samples", "ma_count", "maf",
		    "pval_nominal", "slope", "slope_se", "multiply")

gene_info <- read.table(xargs$geneInfo, header = T, fill = T)
gene_info <- gene_info[,c("gene_id", "gene_name", "gene_chr", "gene_start")] #hardcode

test <- inner_join(test, gene_info, by = c("gene_id" = "gene_id", "#snp_chrom" = "gene_chr"))

a <- test %>% group_split(gene_id)

for (i in 1:length(a)) {
	sub <- a[[i]]
	print(sub[1, 1])
	colnames(sub) <- c("#snp_chrom", "snp_start", "snp_end", "SNP", "REF", "ALT", "ref_gtex", "alt_gtex",
			   "gene_id", "gtex_code", "tss_distance", "ma_samples", "ma_count", "maf",
			   "pval_nominal", "slope", "slope_se", "multiply", "gene_name", "TSS")
	write.table(sub,
              paste0(xargs$outdir, "/", sub$gene_id[1], "__GTEx_Pancreas_Gene__", gsub("chr", "", sub[1, 1]), ":", sub$TSS[1], ".bed"),
              col.names = T, row.names = F, sep = "\t", quote = F)
}

