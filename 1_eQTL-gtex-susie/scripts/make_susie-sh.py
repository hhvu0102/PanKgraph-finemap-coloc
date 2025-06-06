#!/usr/bin/env python3
# coding: utf-8

# In[28]:


import pandas
import numpy
import sys
import os
import glob
import time
import argparse
from string import Template
import yaml
import configargparse

def getOpts():
    parser =  configargparse.ArgParser(config_file_parser_class=configargparse.YAMLConfigFileParser,
                                       description='From a list of GWAS lead SNPs, identify regions to run susie on and make a bash script for each region')
    parser.add('--config', required=True, is_config_file=True, type=yaml.safe_load, help='config file path')
    parser.add('--base', help="""base data dir if paths are not absolute""")
    parser.add_argument('--trait1-info', required=True, type=yaml.safe_load, help="""for each trait1, dictionary of trait type, column names for beta, se, maf etc.""")
    parser.add('--susie-template', required=True, help="""bash template file for susie run script""")
    parser.add('--make-ld-mat-template', required=False, help="""bash template file for make-ld-mat run script""")
    parser.add('--trait1-leads', required=True, help="""glob of headerless 4-column trait 1 bed files with lead SNPs. All should have the same --trait1-info etc.""")
    parser.add('--trait1-ref', required=True, help="""genotype dosage vcf for trait1""")
    parser.add('--trait1-ref-format', required=True, help="""Is the ref vcf chrom is of the format chr10 (chr) or 10 (int)""")
    parser.add('--susie-window', type=int, default=[250000], action="append", help="""Flank window on trait1 lead SNP for SuSiE""") #500kb when using this script for eQTL from Inspire
    parser.add('--dropsets-if-not-contain-lead', action='store_true', default=False, help="""fetch out the signal lead SNP id from the locus name, drop SuSiE sets that don't contain this snp if this parameter is set as True. Default False.""")
    parser.add('--prep-template', required=True, help="""bash template file for susie prep script """) 
    parser.add('--selected-stats', required=True, help="""glob of tsv file containing lead SNPs or TSS coordinates and where to find eqtl_input files. All should have the same --trait1-info etc.""")
    args = parser.parse_args()
    return args

#make-ld-mat-template.sh
def make_ld_mat_sh(x, make_ld_mat_template, sub_vars):
    
    ld_sh_filename = f"{x['susie_locus']}.ld.sh"

    for p in ["eqtl_input"]:
        assert os.path.exists(x[p])
        
    sub_dict = {k: x[k] for k in sub_vars}
    if "exon_id" in x:
        sub_dict['extra_params'] = f" --exon_id {x['exon_id']}"
    else:
        sub_dict['extra_params'] = ""
    with open(make_ld_mat_template, 'r') as f:
        txt = Template(f.read())
        cmds = txt.substitute(sub_dict)
        
    with open(ld_sh_filename, 'w+') as f:
        f.write(cmds)


def make_susie_sh(x, susie_template, sub_vars):

    susie_sh_filename = f"{x['susie_locus']}.susie.sh"

    for p in ["eqtl_input"]:
        assert os.path.exists(x[p])

    sub_dict = {k: x[k] for k in sub_vars}
    # variables to populate the susie.sh template
    #all_params = get_group_val(x, 't1_params')
    #print(all_params)
    #sub_dict.update({'t1_params': all_params})

    if "exon_id" in x:
        sub_dict['extra_params'] = f" --exon_id {x['exon_id']}"
    else:
        sub_dict['extra_params'] = ""
    
    with open(susie_template, 'r') as f:
        txt = Template(f.read())
        cmds = txt.substitute(sub_dict)
        
    with open(susie_sh_filename, 'w+') as f:
        f.write(cmds)


def checkpath(x):
    x = x.format(base = base)
    if os.path.exists(x):
        return(x)
    else:
        print(f"path {x} does not exist")
        sys.exit(1)


def susie_region(chrom, pos, cflank):
    start= pos - 1 - cflank
    start = 0 if start < 0 else start
    end = pos + cflank
    fetch = f"{chrom}:{start}-{end}"
    return fetch

def get_group_val(g, col):
    val = g.iloc[0][col]
    return val


if __name__ == '__main__':
    
    args = getOpts()

    base = args.base

    susie_window = args.susie_window

    selected = checkpath(args.selected_stats)

    trait1_info = args.trait1_info
   
    trait1_leads = args.trait1_leads.format(base = base)
    trait1_lead_files = glob.glob(trait1_leads)
    assert len(trait1_lead_files) > 0 #actually only working for one file for now
    trait1_replace = os.path.basename(trait1_leads).replace("*", "")

    for trait1 in trait1_lead_files:
        trait1_name = os.path.basename(trait1).replace(trait1_replace, "")

    susie_template = checkpath(args.susie_template)
    #make_ld_template = checkpath(args.make_ld_mat_template)
 
    colnames = ['chrom', 'start', 'snp_end', 'locus', 'gene_id', 'eqtl_input']
    d = pandas.read_csv(selected, sep='\t', header = 0, names=colnames,
                        dtype={'snp_end': int, 'start': int, 'window': int})

    #d['t1_params'] = " ".join([f'--{key} {trait1_info[key]}' for key in trait1_info.keys()])
    specialChars = ", !#$%^&*();[]"
    
    def replaceall(string, specialChars, replaceWith):
        for c in specialChars:
            string = string.replace(c, replaceWith)
        
        return string
    
    d['trait1_locus'] = d['locus'].map(lambda x: replaceall(x, specialChars, "_"))
    d['marker'] = d['trait1_locus'].map(lambda x: x.split('__')[1])
    d['trait1_name'] = trait1_name
    d['t1_params'] = d['trait1_name'].map(lambda x: " ".join([f'--{key} {trait1_info[x][key]}' for key in trait1_info[x].keys()]))
    d.drop(['locus'], axis=1, inplace=True)

    d['t1_params'] = d['t1_params'] + " --marker " + d['marker']
    
    if args.dropsets_if_not_contain_lead:
        d['t1_params'] = d['t1_params'] + " --dropsets_if_not_contain " + d['marker']

    if "window" in d.columns:
        d['fetch'] = d.apply(lambda x: susie_region(x['chrom'], x['snp_end'], x['window']), axis=1)
        d['susie_locus'] = d.apply(lambda x: f"{x['trait1_name']}__{x['trait1_locus']}__{x['fetch'].replace(':', '-')}__{int(x['window']/1000)}kb", axis=1)
        d.groupby(['fetch', 'susie_locus']).apply(make_susie_sh, prep_template, susie_template, args.trait1_ref_format)

    else:
        for window in susie_window:
            d['fetch'] = d.apply(lambda x: susie_region(x['chrom'], x['snp_end'], window), axis=1)
            d['susie_locus'] = d.apply(lambda x: f"{x['trait1_name']}__{x['trait1_locus']}__{x['fetch'].replace(':', '-')}__{int(window/1000)}kb", axis=1)
    
    print(d.head())
    print(d['t1_params'][0])
    
    # variables to populate the susieprep.sh template
    sub_vars = ['eqtl_input', 't1_params', 'trait1_locus', 'trait1_name', 'fetch', 'susie_locus', 'eqtl_input', 'marker']
       
    #d.apply(lambda x: make_ld_mat_sh(x, make_ld_template, sub_vars), axis=1) #skipping this and merge this step with Susie to save storage space
    d.apply(lambda x: make_susie_sh(x, susie_template, sub_vars), axis=1)
