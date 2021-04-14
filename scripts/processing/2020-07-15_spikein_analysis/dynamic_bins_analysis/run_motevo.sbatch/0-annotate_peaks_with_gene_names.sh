#!/bin/sh
# Jake Yeung
# 0-annotate_peaks_with_gene_names.sh
#  
# 2020-02-14

jmem='16G'
jtime='1:00:00'

inscript="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/annotate_bed_to_gene_and_distance.R"

[[ ! -e $inscript ]] && echo "$inscript not found, exiting" && exit 1

jmark="H3K27me3"
# dname="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_H3K9me3_dynamic_bins"
dname="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_top_6085_four_marks_dynamic_bins"

# inf="${dname}/coords_H3K9me3_dynamic_bins.noname.2021-01-28.bed"
inf="${dname}/top_6085.${jmark}.2021-02-17.bed"
outf="${dname}/top_6085.${jmark}.2021-02-17.annotated.bed"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ -e $outf ]] && echo "$outf found, continuing" && continue

BNAME=${dname}/AnnotateBeds.sbatch_out
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $inscript $inf $outf"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=AnnotateBeds --wrap "$cmd"
