#!/bin/sh
# Jake Yeung
# 0-test_motevo.sh
#  
# 2020-11-07

jmem='8G'
jtime='1:00:00'

cmd="/hpc/hub_oudenaarden/jyeung/software/motevo_ver1.11/bin/motevo /hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged3_Peaks/H3K27me3/fastasplit/merged.H3K27me3.cutoff_analysis.merged.withchr2.annotated.fa.aa /hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged3_Peaks/H3K27me3/motevo_outputs/split/aa/param_files/Srf.pwm.param /hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm10_weight_matrices_v2_split/Srf.pwm"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged3_Peaks/H3K27me3/motevo_outputs/split/aa"
BNAME=$outdir/testmotevo.SRF.out
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=test --wrap "$cmd"
