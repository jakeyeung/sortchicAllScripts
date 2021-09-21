#!/bin/sh
# Jake Yeung
# 1-make_r1_read_only_K562_Zeller.sh
#  
# 2020-11-27

jmem='16G'
jtime='6:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/K562_merged/correlate_with_public_data/remove_r2_reads.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_G1filt_split_by_G1filt"
outdir="${indir}/r1only"

[[ ! -d $outdir ]] && mkdir $outdir

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
for jmark in $jmarks; do

    BNAME=${outdir}/${jmark}.r1only.sbatchout
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    inf="${indir}/K562_AllMerged_${jmark}.merged.sorted.tagged.G1filt.sorted.bam"
    outf="${outdir}/K562_AllMerged_${jmark}.merged.sorted.tagged.G1filt.r1only.sorted.bam"
   [[ -e $outf ]] && echo "$outf found, continuing" && continue
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
done
