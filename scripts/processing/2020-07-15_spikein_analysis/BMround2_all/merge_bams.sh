#!/bin/sh
# Jake Yeung
# merge_bams.sh
#  
# 2020-10-08

jmem='32G'
jtime='6:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/merged_bams"

cd $indir

jmarks="K4me1 K4me3 K27me3 K9me3"

# K4me1
# jmark="K4me1"

for jmark in $jmarks; do
  outf="${outdir}/PZ-BM-rep2rep3-${jmark}.merged.sorted.bam"
  infs=$(ls -d ${indir}/*.bam | tr "\n" " ")
  
  BNAME=$outdir/${jmark}_merge.log
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outf $infs; samtools index $outf"
  echo $cmd
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${outf}_%j.log --ntasks=1 --nodes=1 --job-name=merge_${jmark} --wrap "$cmd"
done



