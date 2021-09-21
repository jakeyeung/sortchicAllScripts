#!/bin/sh
# Jake Yeung
# 10-merge_bams_by_ctype.sh
# Merge plates together
# 2020-11-29

jmem='16G'
jtime='2:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/split_by_cluster.MAPQ_40.H3K27me3reseq/H3K27me3"

outdir="${indir}/merged_by_ctype"
[[ ! -d $outdir ]] && mkdir $outdir

cd $indir

ctypes=$(ls -d *.bam | cut -d"." -f4 | sort | uniq | tr "\n" " ")
echo $ctypes

for ctype in $ctypes; do
  echo $ctype
  outname="PZ-BM-rep3-H3K27me3-platemerged.${ctype}"
  outbam=${outdir}/${outname}.bam

  BNAME=${outdir}/${outname}.sbatchout
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

  # PZ-BM-rep3-H3K27me3-9.sorted.tagged.DCs.sorted.bam
  inbams=$(ls -d PZ-BM-rep3-H3K27me3-*.sorted.tagged.${ctype}.sorted.bam |  tr "\n" " ")
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outbam $inbams; samtools index $outbam"
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${ctype} --wrap "$cmd"
done

