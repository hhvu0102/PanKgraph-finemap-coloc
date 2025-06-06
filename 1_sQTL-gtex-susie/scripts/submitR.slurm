#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=summry
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20Gb
#SBATCH --time=4:00:00
#SBATCH --account=scjp99
#SBATCH --partition=standard
#SBATCH --mail-user=vthihong@umich.edu
#SBATCH --output=../reports/summary.%j.out   # STDOUT output file
#SBATCH --error=../reports/summary.%j.err    # ERR output file

ml R/4.3.1
Rscript /scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/1_eQTL-inspire-susie/scripts/getSummary.R
