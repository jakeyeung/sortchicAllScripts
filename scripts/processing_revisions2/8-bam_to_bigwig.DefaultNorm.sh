#!/bin/sh
# Jake Yeung
# 2-bam_to_bigwig.sh
#  
# 2020-01-10

# # WRAP UP
# while [[ `qstat | grep split | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='16G'
jtime='12:00:00'


# bs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.offset.DefaultNorm.sh"
bs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/old_new_bams_merged/split_by_celltype"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

# bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
bl="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/mm10.blacklist.copy.sorted.merged.nospikeins.nochromo.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

outmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/bigwigs_by_celltype"

bsizes="500"
for jmark in $jmarks; do
  for bsize in $bsizes; do
	  outdir=${outmain}/${jmark}
	  [[ ! -d $outdir ]] && mkdir $outdir
	  indir=${inmain}/${jmark}
      for inbam in `ls -d $indir/*.bam`; do
          echo $inbam
          bname=$(basename $inbam)
          bname=${bname%.*}  # strip extension 
          outbw=${outdir}/${bname}.${bsize}.CPM.bw
          [[ -e $outbw ]] && echo "$outbw found, continuing" && continue
  
          BNAME=$outdir/$bname.qsub
          DBASE=$(dirname "${BNAME}")
          [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
          cmd=". /nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh; conda activate scmo2022; bash $bs $inbam $outbw $bsize $bl"
		  # echo $cmd
          sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
      done
  done
done
