#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tidyr")

parser <- ArgumentParser(description= 'Translate GWAS summary stats from hg19 to hg38.\n
                                       Some part of this script is hard-coding col names, so check carefully.')
parser$add_argument('--ukbb_hg38', '-u', help= 'VCF file on hg38 containing SNPs. Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/ukbb-rsid/ukbb-rsid.vcf.gz')
parser$add_argument('--summStat', '-t', help= 'Text file containing GWAS signals and summary stats in hg19. Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/temp-ukbb/chr{chrom}_Gene.hg19.bed"')
parser$add_argument('--bed_hg38', '-b', help= 'Bed file containing GWAS signal locations in hg38, obtained using LiftOver. Example: /scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/temp-ukbb/new.chr{chrom}_{analysis}.hg38')
parser$add_argument('--output', '-o', help= 'Output file.')

xargs <- parser$parse_args()

library(data.table)
library(dplyr)

ukbb_hg38 <- fread(xargs$ukbb_hg38, skip = "CHROM")

summStat <- fread(xargs$summStat, header = F)
tmp <- which(colSums(is.na(summStat))<nrow(summStat))
summStat <- as.data.frame(summStat)
summStat <- summStat[, names(tmp)]

if (length(grep("Exons", xargs$summStat)) > 0) {
	colnames(summStat) <- c("GeneID", "ExonsID", "GeneChr",	"TSS",	"DistanceToGene", "rsID", "SNPchr", "SNPposition", "other_allele", "effect_allele", "FreqREF", "FreqALT", "Pvalue", "Slope", "Lead", "T_stat", "SE")
} else {
	colnames(summStat) <- c("GeneID", "GeneChr", "TSS",  "DistanceToGene", "NumVarTested", "rsID", "SNPchr", "SNPposition", "other_allele", "effect_allele", "FreqREF", "FreqALT", "Pvalue", "Slope", "Lead", "T_stat", "SE")
}
summStat$GeneChr <- as.character(summStat$GeneChr)
summStat$GeneChr <- paste0("chr", as.character(summStat$GeneChr))

bed_hg38 <- read.table(xargs$bed_hg38, header = F)
colnames(bed_hg38) <- c("chr", "start", "end", "id")
bed_hg38 <- separate(bed_hg38, "id", c("snp", "GeneID"), sep = "_E")
bed_hg38$GeneID <- paste0("E", bed_hg38$GeneID)
bed_hg38 <- distinct(bed_hg38)
test <- inner_join(as.data.frame(ukbb_hg38[,1:5]), bed_hg38, by = c("#CHROM" = "chr", "POS" = "end")) #1:5 here means we will exclude QUAL, FILTER and INFO columns of ukbb file
test2 <- inner_join(test, summStat, by = c("#CHROM" = "GeneChr", "snp" = "rsID", "GeneID" = "GeneID"))
#needed info - something like the following, depending on what GWAS file gave:
##snp_chrom	snp_start	snp_end	SNP	REF	ALT	minor_allele	maf	n_complete_samples	AC_hg19	ytx_hg19
##beta	se	tstat	pval	SNP_hg19	multiply

#if reference allele is flipped, also flip the slope
test2$multiply <- ifelse(test2$REF == test2$effect_allele & test2$ALT == test2$other_allele, -1,
                         ifelse(test2$ALT == test2$effect_allele & test2$REF == test2$other_allele, 1, NA))
test2 <- test2[!is.na(test2$multiply),] #on lines where none of the cases in ifelse condition holds, drop these lines
test2$Slope <- test2$Slope * test2$multiply
if (length(grep("Exons", xargs$summStat)) > 0) {
  test2 <- test2[, c("#CHROM", "start", "POS", "ID", "REF", "ALT", "effect_allele", "other_allele",
                     "GeneID", "ExonsID", "TSS",  "DistanceToGene", "FreqREF", "FreqALT", "Pvalue", "Slope", "Lead", "T_stat", "SE", "multiply")]
  colnames(test2) <- c("snp_chrom", "snp_start", "snp_end", "SNP",  "REF",  "ALT",  "effect_allele", "other_allele",
                       "GeneID", "ExonsID", "TSS",  "DistanceToGene", "FreqREF", "FreqALT", "Pvalue", "Slope", "Lead", "T_stat", "SE", "multiply")
} else {
  test2 <- test2[, c("#CHROM", "start", "POS", "ID", "REF", "ALT", "effect_allele", "other_allele",
                     "GeneID", "TSS",  "DistanceToGene", "FreqREF", "FreqALT", "Pvalue", "Slope", "Lead", "T_stat", "SE", "multiply")]
  colnames(test2) <- c("snp_chrom", "snp_start", "snp_end", "SNP",  "REF",  "ALT",  "effect_allele", "other_allele",
                       "GeneID", "TSS",  "DistanceToGene", "FreqREF", "FreqALT", "Pvalue", "Slope", "Lead", "T_stat", "SE", "multiply")
}

write.table(test2, xargs$output, col.names = F, row.names = F, quote = F, sep = "\t")
