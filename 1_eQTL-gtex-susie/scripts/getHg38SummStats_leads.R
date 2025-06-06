#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tidyr")

parser <- ArgumentParser(description= 'Adapt GTEx summary stats of lead gene-SNP pairs.\n
                                       Some part of this script is hard-coding col names, so check carefully.')
parser$add_argument('--ukbb_hg38', '-u', help= 'VCF file on hg38 containing SNPs. Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/ukbb-rsid/ukbb-rsid.vcf.gz')
parser$add_argument('--summStat', '-t', help= 'Text file containing eQTL summ stat. Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/0_data/GTEx_all/Pancreas.v8.egenes.txt.gz"')
parser$add_argument('--geneInfo', '-g', help= 'Text file containing gene name and coor in hg38.\n
                                              Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/0_data/GTEx_all/Pancreas.v8.egenes.txt.gz"')
parser$add_argument('--output', '-o', help= 'Output file.')

xargs <- parser$parse_args()

library(data.table)
library(dplyr)

options(scipen=999)

ukbb_hg38 <- fread(xargs$ukbb_hg38, skip = "CHROM")

summStat <- fread(xargs$summStat, header = T)
summStat <- summStat[summStat$qval <= 0.05,]
colnames(summStat)[1] <- "gene_id"
summStat <- summStat[, c("gene_id", "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df", "variant_id", "tss_distance", "maf", "ref_factor", "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta", "qval", "pval_nominal_threshold")]
summStat[, c("snp_chr", "snp_stop", "ref_gtex", "alt_gtex", "gtex_code") := tstrsplit(variant_id, "_")]
summStat <- as.data.frame(summStat)
summStat$snp_chr <- as.character(summStat$snp_chr)
summStat$snp_stop <- as.numeric(summStat$snp_stop)
summStat$snp_start <- summStat$snp_stop-1

ukbb_hg38 <- fread(xargs$ukbb_hg38, skip = "CHROM")
test <- inner_join(as.data.frame(ukbb_hg38[,1:5]), summStat, by = c("#CHROM" = "snp_chr", "POS" = "snp_stop")) #1:5 here means we will exclude QUAL, FILTER and INFO columns of ukbb file
#needed info - something like the following, depending on what GWAS file gave:
##snp_chrom     snp_start       snp_end SNP     REF     ALT     minor_allele    maf     n_complete_samples      AC_hg19 ytx_hg19
##beta  se      tstat   pval    SNP_hg19        multiply

### check matching between ref/alt alleles in different genome built
test$multiply <- ifelse(test$REF == test$alt_gtex & test$ALT == test$ref_gtex, -1,
                         ifelse(test$ALT == test$alt_gtex & test$REF == test$ref_gtex, 1, NA))
test <- test[!is.na(test$multiply),] #on lines where none of the cases in the above condition holds, drop these lines
test$slope <- test$slope * test$multiply
test <- test[, c("#CHROM", "snp_start", "POS", "ID", "REF", "ALT", "ref_gtex", "alt_gtex",
                 "gene_id", "gtex_code", "tss_distance", 
		 "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df", "maf", "ref_factor", "pval_nominal", 
		 "slope", "slope_se", "pval_perm", "pval_beta", "qval", "pval_nominal_threshold", "multiply")]
colnames(test) <- c("#snp_chrom", "snp_start", "snp_end", "SNP", "REF", "ALT", "ref_gtex", "alt_gtex",
		    "gene_id", "gtex_code", "tss_distance",
		    "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df", "maf", "ref_factor", "pval_nominal",
		    "slope", "slope_se", "pval_perm", "pval_beta", "qval", "pval_nominal_threshold", "multiply")

gene_info <- read.table(xargs$geneInfo, header = T, fill = T)
gene_info <- gene_info[,c("gene_id", "gene_name", "gene_chr", "gene_start")] #hardcode, only true for this GTEx file

test <- inner_join(test, gene_info, by = c("gene_id" = "gene_id", "#snp_chrom" = "gene_chr"))


write.table(test, xargs$output, col.names = T, row.names = F, quote = F, sep = "\t")

#to index this file, we may have to sort: head -n 1 eQTL_EUR_leads.txt > header; tail -n +2 eQTL_EUR_leads.txt | sort -k1,1 -k2,2n > tmp; cat header tmp > eQTL_EUR_leads.txt before bgzip and index
