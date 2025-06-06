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


complements <- list("A" = "T",
                    "T" = "A",
                    "G" = "C",
                    "C" = "G")

revcomp <- function(x) {
  out <- sapply(strsplit(x, "")[[1]], function(i) complements[[i]])
  return(paste(out, collapse = ""))
}

summStat <- fread(xargs$summStat, header = T)
#phenotype_id	variant_id	tss_distance	maf	ma_samples	ma_count	pval_nominal	slope	slope_se	snp_chr	snp_stop	ref_gtex	alt_gtex	gtex_code	pheno_chr	pheno_start	pheno_stop	pheno_clu	gene_id	snp_start

print("start merging summ stat file and vcf file")

ukbb_hg38 <- fread(xargs$ukbb_hg38, skip = "CHROM")
test <- inner_join(as.data.frame(ukbb_hg38[,1:5]), summStat, by = c("#CHROM" = "snp_chr", "POS" = "snp_stop")) #1:5 here means we will exclude QUAL, FILTER and INFO columns of ukbb file
#needed info - something like the following, depending on what GWAS file gave:
##snp_chrom	snp_start	snp_end	SNP	REF	ALT	minor_allele	maf	n_complete_samples	AC_hg19	ytx_hg19
##beta	se	tstat	pval	SNP_hg19	multiply

### check matching between ref/alt alleles in different genome built
# REF/ATL columns are ref/alt in reference genome, in this case hg38
# effect_allele/other_allele/effect_allele_frequency are in the original summ stat file, in this case hg19
#revEA <- unlist(lapply(test2$effect_allele, revcomp)) #obtain reverse compliments of the alleles
#revOA <- unlist(lapply(test2$other_allele, revcomp))

#if alt allele = reverse compliment of effect alleles & alt allele = reverse compliment of effect alleles, sign(beta) remains the same
test$multiply <- ifelse(test$REF == test$alt_gtex & test$ALT == test$ref_gtex, -1,
                         ifelse(test$ALT == test$alt_gtex & test$REF == test$ref_gtex, 1, NA))
                                #ifelse(test2$ALT == revEA & test2$REF == revOA, 1, NA))) 
test <- test[!is.na(test$multiply),] #on lines where none of the cases in the above condition holds, drop these lines
test$slope <- test$slope * test$multiply
test <- test[, c("#CHROM", "snp_start", "POS", "ID", "REF", "ALT", "ref_gtex", "alt_gtex",
		 "gene_id", "gtex_code", "tss_distance", "ma_samples", "ma_count", "maf",
		 "pval_nominal", "slope", "slope_se", "phenotype_id", "multiply")]
colnames(test) <- c("#snp_chrom", "snp_start", "snp_end", "SNP", "REF", "ALT", "ref_gtex", "alt_gtex",
		    "gene_id", "gtex_code", "tss_distance", "ma_samples", "ma_count", "maf",
		    "pval_nominal", "slope", "slope_se", "phenotype_id", "multiply")

gene_info <- read.table(xargs$geneInfo, header = T, fill = T)
gene_info <- gene_info[,c("gene_id", "gene_name", "gene_chr", "gene_start")] #hardcode

test <- inner_join(test, gene_info, by = c("gene_id" = "gene_id", "#snp_chrom" = "gene_chr"))

a <- test %>% group_split(phenotype_id)

for (i in 1:length(a)) {
	sub <- a[[i]]
	print(sub[1, 1])
	colnames(sub) <- c("#snp_chrom", "snp_start", "snp_end", "SNP", "REF", "ALT", "ref_gtex", "alt_gtex",
			   "gene_id", "gtex_code", "tss_distance", "ma_samples", "ma_count", "maf",
			   "pval_nominal", "slope", "slope_se", "phenotype_id", "multiply", "gene_name", "TSS")
	write.table(sub,
              paste0(xargs$outdir, "/", sub$gene_id[1], "__", sub$phenotype_id[1], "__GTEx_Pancreas_sGene__", gsub("chr", "", sub[1, 1]), ":", sub$TSS[1], ".bed"),
              col.names = T, row.names = F, sep = "\t", quote = F)
}

