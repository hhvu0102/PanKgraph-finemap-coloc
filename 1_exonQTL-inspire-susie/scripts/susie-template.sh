#!/bin/bash

################## running SuSiE for  ${susie_locus}:

## Susie 
/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/1_eQTL-inspire-susie/scripts/susie-eqtl.R --prefix ${susie_locus} ${t1_params} --trait1 ${eqtl_input} --trait1_ld ${susie_locus}.ukbb-dosages.tsv.gz ${extra_params}
		
