#!/usr/bin/env Rscript

library(susieR)
library(coloc)
library(glue)
library(tidyr)
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
library(optparse)
library(gtable)
library(grid)
library(gridExtra)
library(cowplot)
ggplot2::theme_set(theme_cowplot())
library(locuscomparer)

plot.set = function(setlist, trait, input, title, marker = NULL){
    print(marker)
    input$set = "none"
    for (i in seq_along(setlist)){
        input$set = ifelse(input$snp %in% names(setlist[[i]]), glue("set{i}"), input$set)
    }
    
    if (trait == "gwas"){
        p = ggplot(input, aes(x=snp_end, y=-log10(pval))) +
            geom_point(aes(color=set))
    } else {
        p = ggplot(input, aes(x=snp_end, y=-log10(p_nominal))) +
            geom_point(aes(color=set))
    }

    if (! is.null(marker)){
        p = p +
            geom_point(data = input[input$snp == marker,], shape=1, color="blue")
    }
    p = p +
        labs(x="position", title=title) +
        theme(legend.position="bottom", plot.title=element_text(size=7))
    return(p)
}

plot.result = function(result.vector, name, df, marker = NULL){
    new = list()
    new[[name]] = result.vector
    df = cbind(df, data.frame(new))
    p = ggplot(df, aes_string(x="snp_end", y=name)) +
        geom_point(shape=16, color="grey")
    if (! is.null(marker)){
        p = p +
            geom_point(data = df[df$snp == marker,], shape=1, color="blue")
    }
    
    return(p)
}

get_var = function(df, trait, betacol, secol, pcol){
    vcol = glue("{trait}_var")
    zcol = glue("{trait}_z")

    if (is.null(secol)) {
        df[[zcol]] = (df[[betacol]]/abs(df[[betacol]]))*qnorm(df[[pcol]]/2)
        df[[vcol]] = (df[[betacol]]/df[[zcol]])^2
    } else {
        df[[vcol]] = df[[secol]]^2
        df[[zcol]]=df[[betacol]]/df[[secol]]
    }
    return(df)
}

make_coloc_dataset = function(df, beta, var, type, position, snp, sdY=NULL, maf=NULL, N=NULL, LD=NULL) {
    D = list(beta = df[[beta]], varbeta = df[[var]],
             snp = df[[snp]], position = df[[position]], type = type, LD = LD)

    if (type == "quant"){
        if ( ! is.null(sdY)){
            D = c(D, list(sdY = sdY))
        } else {
            stopifnot((! is.null(maf)) & (! is.null(N)) )
            D = c(D, list(MAF = df[[maf]], N = df[[N]]))
        }
    }
    return(D)
}

susie_dropsets = function(res, drop){
    res$sets$cs = res$sets$cs[-drop]
    res$sets$cs_index = res$sets$cs_index[-drop]
    res$sets$coverage = res$sets$coverage[-drop]
    res$sets$purity = res$sets$purity[-drop,,drop=FALSE]
    return(res)
}

process_dosage = function(f, snplist){
    ld = read.csv(f, sep='\t', check.names = F)
    dups = ld[ (duplicated(ld$ID) | duplicated(ld$ID, fromLast = TRUE)),]
    print(glue("N duplicates = {nrow(dups)}"))
    ld = ld[! (duplicated(ld$ID) | duplicated(ld$ID, fromLast = TRUE)),]
    row.names(ld) = ld$ID
    ld$ID = NULL
    idlist = intersect(snplist, row.names(ld))
    ld = ld[idlist,]
    print(ld[1:5, 1:10])
    ld = cor(t(ld))
    return(ld)
}

option_list <- list(
    make_option(c("--trait1_ld"), type = "character", default=NULL, help = "dosage file to calculate correlation matrix for trait1 reference panel. [DEFAULT = NULL]"),
    make_option(c("--trait1"), type = "character", help = "[Required] trait2 dataframe"),
    make_option(c("--exon_id"), type = "character", default=NULL, help = "select this id from the exon_id column for trait2"),
    make_option(c("--trait1_name"), type = "character", help = "[Required] trait2 name"),
    make_option(c("--type"), type = "character", help = "[Required] trait2 type-  cc or quant"),
    make_option(c("--beta"), type = "character", help = "[Required] column name for trait2 beta"),
    make_option(c("--se"), type = "character", help = "column name for trait2 se. If not available, provide p value column in --trait1_p"),
    make_option(c("--effect"), type = "character", help = "column name for trait2 effect allele"),
    make_option(c("--non_effect"), type = "character", help = "column name for trait2 non effect allele"),
    make_option(c("--p"), type = "character", help = "column name for p value - will calculate z and se if se is not available for trait1"),
    make_option(c("--sdY"), type = "numeric", help = "If trait1 is quant, sdY (phenotype). 1 if the phenotypes were inverse normalized"),
    make_option(c("--maf"), type = "character", help = "if trait1 is quant and sdY is not known, provide the column name for maf"),
    make_option(c("--N"), type = "character", help = "if trait1 is quant and sdY is not known, provide the column name for N"),
    make_option(c("--coverage"), type = "numeric", default=0.95, help = "SuSiE coverage [DEFAULT = 0.95]"),
    make_option(c("--r2_prune"), type = "numeric", default=0.8, help = "SuSiE r2.prune [DEFAULT = 0.8]"),
    make_option(c("--maxit"), type = "numeric", default=10000, help = "SuSiE max iterations [DEFAULT = 10000]"),
    make_option(c("--min_abs_corr"), type = "numeric", default=0.1, help = "SuSiE min_abs_corr [DEFAULT = 0.1]"),
    make_option(c("--number_signals_default"), type = "numeric", default=10, help = "SuSie L [DEFAULT = 10]"),
    make_option(c("--s_threshold"), type = "numeric", default=0.3, help = "max s threshold above which change L [DEFAULT = 0.3]"),
    make_option(c("--number_signals_high_s"), type = "numeric", default=1, help = "SuSiE L if s > max s_threshold [DEFAULT = 1]"),
    make_option(c("--reduce_corr"), type = "numeric", default=NULL, help = "If no credible sets are found, sequentially reduce the min_abs_corr by 0.05 until this value is reached and see if sets can be found. [Default = NULL]"),
    make_option(c("--reduce_coverage"), type = "numeric", default=NULL, help = "If no credible sets are found, sequentially reduce the coverage by 0.05 until this value is reached and see if sets can be found. [Default = NULL]"),
    make_option(c("--marker"), type = "character", default=NULL, help = "highlight in red an inteded lead SNP by providing rsid"),
    make_option(c("--dropsets_if_not_contain"), type = "character", default=NULL, help = "drop the credible set if it doesn't contain this rsid"),
    make_option(c("--dropsets_if_no_proxy"), type = "numeric", default=NULL, help = "drop the credible set if it doesn't contain this rsid or at least 1 proxy with this r2 or greater. Flag --dropsets_if_not_contain if provided takes precedent over this."),
    make_option(c("--genome"), default = "hg19", type = "character", help = "genome. [DEFAULT = hg19]"),
    make_option(c("--prefix"), type = "character", help = "output prefix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

## read vars
prefix = opts$prefix
t2.beta = opts$beta
t2.se = opts$se
t2.p = opts$p
t2.type = opts$type
t2.sdY = opts$sdY
t2.maf = opts$maf
t2.N = opts$N
t2.EA = opts$effect
t2.NEA = opts$non_effect


## build LD matrix
print("start processing LD matrix")
## runsusie GWAS only technically needs beta and se since it is case control, but for calculating the s metric, we also need Z scores. i.e using p values column.
## read trait1
input1 = read.csv(opts$trait1, sep='\t', header=T, as.is=T, check.names=F)
snps = paste0(input1$SNP, "-", input1$REF, "-", input1$ALT)
print(head(input1))

ldf = process_dosage(opts$trait1_ld, snps)
snps = colnames(ldf)
print("length(snps)")
print(length(snps))
rslist_dup = duplicated(gsub("-.*", "", snps))
rslist = gsub("-.*", "", snps)
rslist <- setdiff(rslist, rslist[rslist_dup])
print(head(rslist))
print(length(rslist))

s = data.frame(snps)
colnames(s) =  "snp"
s = tidyr::separate(data = s, col = snp, into = c("rsid", "REF", "ALT"), sep = "-")
print("head(s) before filter")
print(head(s))
s <- s[s$rsid %in% rslist,]
print("dim(s)")
print(dim(s))


## read trait2 (call it trait 2 because we are analyzing eQTL summ stats, but the argument is still trait1
input2 = read.csv(opts$trait1, sep='\t', header=T, check.names=F)
if (! is.null(opts$exon_id)){
    if ("ExonsID" %in% colnames(input2)){
        input2 = input2[input2$ExonsID == opts$exon_id,]
    }
}

colnames(input2)[which(colnames(input2) %in% c("SNP", "Pvalue", "Slope", "SE") == TRUE)] <- c("rsid", "p_nominal", "beta", "se")

input2 <- distinct(input2)
print(dim(input2))
print(head(input2))
required.columns = c(t2.beta, t2.p, "rsid", t2.EA, t2.NEA)
for (i in required.columns){
    stopifnot(i %in% colnames(input2))
}
## check for snp and allele consistency
newt2.EA = "t2.EA"
newt2.NEA = "t2.NEA"

input2 = input2 %>% rename(!!newt2.EA := t2.EA)
input2 = input2 %>% rename(!!newt2.NEA := t2.NEA)

input2 = input2[input2$rsid %in% rslist,]
input2 = merge(input2, s, by="rsid")
print(head(input2))
input2 = input2[((input2$REF == input2[[newt2.NEA]]) & (input2$ALT == input2[[newt2.EA]])) | ((input2$REF == input2[[newt2.EA]]) & (input2$ALT == input2[[newt2.NEA]])),]
input2[[t2.beta]] = ifelse(input2$ALT == input2[[newt2.NEA]], -1*input2[[t2.beta]], input2[[t2.beta]])
input2$snp = glue("{input2$rsid}-{input2$REF}-{input2$ALT}")
duplicated.input = input2[duplicated(gsub("-.*", "", input2$snp), fromLast=TRUE), ]

if  (nrow(duplicated.input) > 0){
    print(glue("WARNING: {nrow(duplicated.input)} duplicates found in data. Removing for now. Check data."))
    print(duplicated.input)
    input2 = dplyr::setdiff(input2, duplicated.input)
}

print("using input2$snp as rownames")
rownames(input2) = input2$snp
input2 = get_var(input2, "trait2", betacol = t2.beta, secol = t2.se, pcol = t2.p)

## update input df - in case some variants were missing from the ld matrix.
snplist = intersect(snps, input2$snp)
input2 = input2[snplist,]
print(head(input2))
ldf = ldf[snplist, snplist]

print(nrow(input2))
print(nrow(ldf))

t2.ld = as.matrix(ldf)
gsub("-.*", "", rownames(ldf))[duplicated(gsub("-.*", "", rownames(ldf)))]
rownames(ldf) = gsub("-.*", "", rownames(ldf))
colnames(ldf) = gsub("-.*", "", colnames(ldf))
## make coloc dataset
position = "snp_end"
snp = "snp"
D2 = make_coloc_dataset(input2, beta = t2.beta, var = "trait2_var", type = t2.type,
                        position = position, snp = snp, sdY = t2.sdY, maf = t2.maf,
                        LD = t2.ld)

## Set default metrics
coverage = opts$coverage
r2.prune = opts$r2_prune
maxit = opts$maxit
min_abs_corr = opts$min_abs_corr

## Calculate diagnostic
diagnostic_s_gwas = estimate_s_rss(input2$trait2_z, t2.ld)
print(glue("s = {diagnostic_s_gwas}"))

## plot diagnostic
png(glue("{prefix}.s.png"), height=4, width=4, units="in", res=150)
condz_in = kriging_rss(input2$trait2_z, t2.ld)
print(condz_in$plot)
dev.off()

## if diagnostic s is very bad, set max signals to 1
if(diagnostic_s_gwas > opts$s_threshold){
    print("Setting numbersignals to 1")
    numbersignals = opts$number_signals_high_s
} else {
    numbersignals = opts$number_signals_default
}

## runsusie 
S2 = coloc::runsusie(D2, coverage = coverage, maxit = maxit,
                     min_abs_corr = min_abs_corr, repeat_until_convergence = F, L = numbersignals, check_prior=FALSE)
    
## if fail to converge, reduce max number of signals
while (S2$converged == F  & (numbersignals > 0)) {
    numbersignals = numbersignals - 1
    S2 = coloc::runsusie(D2, coverage = coverage, r2.prune = r2.prune, maxit = maxit,
                         min_abs_corr = min_abs_corr, repeat_until_convergence = F, L = numbersignals, check_prior=FALSE)
    if (S2$converged) {
        print(glue("Converged after considering max signals = {numbersignals}"))
    } else {
        print(glue("Didn't converge with max signals = {numbersignals}"))
    }
} # end max signals while loop
    
## if no signals found, reduce minimum abs correlation if that option was selected
if (!is.null(opts$reduce_corr)){
    while ((is.null(S2$sets$cs)) & (min_abs_corr >= opts$reduce_corr + 0.05) & (S2$converged == T)) {
        min_abs_corr = min_abs_corr - 0.05
        S2 = coloc::runsusie(D2, coverage = coverage, r2.prune = r2.prune, maxit = maxit,
                             min_abs_corr = min_abs_corr, repeat_until_convergence = F, L = numbersignals, check_prior=FALSE)
        if (! is.null(S2$sets$cs)){
            print(glue("Csets found after considering min_abs_corr = {min_abs_corr}"))
        } else {        
            print(glue("No csets found with min_abs_corr = {min_abs_corr}"))
        }
    } # end min abs correlation loop
}

print("S2$sets$cs")
print(S2$sets$cs)
## if still no signals found, reduce coverage sequentially if that option was selected
if (!is.null(opts$reduce_coverage)){
    step = 0.05
    while ((is.null(S2$sets$cs)) & (round(coverage - step, 2) >= opts$reduce_coverage) & (S2$converged == T)) {
        coverage = coverage - step
        print(glue("will try to find signals after reducing coverage to {opts$reduce_coverage}; current coverage = {coverage}"))
        S2 = coloc::runsusie(D2, coverage = coverage, r2.prune = r2.prune, maxit = maxit,
                             min_abs_corr = min_abs_corr, repeat_until_convergence = F, L = numbersignals, check_prior=FALSE)
        if (! is.null(S2$sets$cs)){
            print(glue("Csets found after considering coverage = {coverage}"))
        } else {
            print(glue("No csets found with coverage = {coverage}"))
        }
    }
}


## write susie results
sfname_gwas = glue("{prefix}.susie.Rda")
save(S2, file=sfname_gwas) 

results = c(glue("prefix\tdiagnostic_s_gwas\tmin_abs_corr\tnumbersignals\tcoverage\tnsets\tlensets\tn_selected_sets\tlen_selected_sets"))

if (is.null(S2$sets$cs)){
    nsets = 0
    lensets = 0
    n_selected_sets = 0
    len_selected_sets = 0
} else {
    nsets = length(S2$sets$cs)
    lensets = paste(lapply(S2$sets$cs, function(i){length(i)}), collapse = ',')
    n_selected_sets = nsets
    len_selected_sets = lensets
}


## color csets over GWAS P values
## intended lead SNP:
if (! is.null(opts$marker)){
    marker.id = input2[gsub("-.*", "", input2$snp) == opts$marker,]$snp
} else {
    marker.id = NULL
}
title = glue("{prefix}\ns = {diagnostic_s_gwas}, min_corr = {min_abs_corr}\n L = {numbersignals}, coverage = {coverage}\nnsets = {nsets}")
p1 = plot.set(S2$sets$cs, "eqtl",  input2, title, marker = marker.id)

## plot PIP
p2 = plot.result(S2$pip, "pip", input2, marker = marker.id)
## g2 = ggplotGrob(p2)

png(glue("{prefix}.png"), height=6, width=8, units="in", res=150)
cowplot::plot_grid(p1, p2, align = 'v', ncol = 1, rel_heights = c(0.6, 0.4)) ## display plotpip
dev.off()


results = c(results, glue("{prefix}\t{diagnostic_s_gwas}\t{min_abs_corr}\t{numbersignals}\t{coverage}\t{nsets}\t{lensets}\t{n_selected_sets}\t{len_selected_sets}"))
writeLines(results, glue("{prefix}.results.tsv"))
