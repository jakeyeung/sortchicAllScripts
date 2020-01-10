#!/bin/sh
# Jake Yeung
# run_snakemake.sh
#  
# 2019-11-26

maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all"

# edit config to look like: 
# {
#     "reference_file"  :  "/hpc/hub_oudenaarden/group_references/ensembl/98/danio_rerio/Danio_rerio.GRCz11.dna.primary_assembly.bgzipformat.fa.gz"
#     "counting_min_mq" : 40,
#     "counting_bin_sizes": [100_000],
#     "counting_sliding_window_increments": [20_000],
#     "mapper":"bwa"
# }

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; cd $maindir; submission.py "snakemake --cluster sge_wrapper.py  --jobs 20 --restart-times 3" -y -time 50
