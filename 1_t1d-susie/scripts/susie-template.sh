#!/bin/bash

################## running SuSiE for  ${susie_locus}:

## Susie 
/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/1_t1d-susie/scripts/susie-gwas.R --trait1 ${susie_locus}.gwas.tsv --prefix ${susie_locus} --ld_mat ${susie_locus}.ld.tsv  ${t1_params} ;
		
