#!/bin/sh
# Jake Yeung
# 2-append_genome_with_lambdaphage.sh
#  
# 2020-07-15
# from 
# https://www.royfrancis.com/read-counts-of-rna-seq-spike-ins/

# inforig="/hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92.fa"
inforig="/hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/primary_assembly_NOMASK_ERCC92.fa"
infnew="/hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/primary_assembly_NOMASK_ERCC92.WithLambdaPhage.fa"
# infnew="/hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92_WithLambdaPhage.fa"
inflambda="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/genomes/LambdaPhageGenome.fa"

[[ ! -e $inforig ]] && echo "$inforig not found, exiting" && exit 1
[[ -e $infnew ]] && echo "$infnew found, exiting" && exit 1
[[ ! -e $inflambda ]] && echo "$inflambda not found, exiting" && exit 1

cp $inforig $infnew

cat $inflambda >> $infnew

