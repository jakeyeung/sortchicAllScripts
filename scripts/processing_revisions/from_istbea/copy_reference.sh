#!/bin/sh
# Jake Yeung
# copy_reference.sh
#  
# 2022-01-22

inf="/hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/primary_assembly_97_129B6Masked_ERCC92.fa*"
outf="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/genomes"

scp t2:$inf $outf
