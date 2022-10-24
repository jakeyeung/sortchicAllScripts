#!/bin/sh
# Jake Yeung
# 5-make_count_tables_from_tagged_bams.sh
#  
# 2021-06-27


# # WRAP UP
# while [[ `squeue -u jyeung | wc -l` > 1 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='64G'
jtime='24:00:00'

# inbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/tagged_bams"
inbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/old_new_bams_merged"
[[ ! -d $inbase ]] && echo "$inbase not found, exiting" && exit 1

mapq=40
bl="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/mm10.blacklist.copy.sorted.merged.nospikeins.nochromo.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

binsizes="50000"

# outbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables"
outbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_bins_redo"
[[ ! -d $outbase ]] && mkdir $outbase

experi="BM"
# jmarks="k4me1-k9me3 k27me3 k9me3"
# jmarks="k4me1 k4me3 k27me3 k9me3"
jmarks="k4me1 k4me3 k27me3 k9me3"
# jmark="H3K36me3"

declare -A markarray

markarray[k4me1]=H3K4me1
markarray[k4me3]=H3K4me3
markarray[k27me3]=H3K27me3
markarray[k9me3]=H3K9me3


for jmark in $jmarks; do 
  echo $jmark
  for binsize in $binsizes; do
      stepsize=$binsize
  
      outmain="${outbase}/counts_tables_${binsize}"
      [[ ! -d $outmain ]] && mkdir $outmain
      outdir="${outmain}/${experi}_${jmark}"
      [[ ! -d $outdir ]] && mkdir $outdir

     # inmain="${inbase}/${experi}_${jmark}"
	 inmain=${inbase}
     [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
     # for inbam in `ls -d $inmain/BM_allmerged_${jmark}.bam`; do
	 # printf "%s is in %s\n" "$c" "${continent[$c]}"
	 # https://stackoverflow.com/questions/17403498/iterate-over-two-arrays-simultaneously-in-bash
	 
	 jmarkold=${markarray[${jmark}]}
	 echo $jmarkold
	 inbam="${inmain}/BM_allmerged_${jmarkold}.bam"
	 echo $inbam
	      [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
          inbambase=$(basename $inbam)
          inbambase=${inbambase%.*}
          bname=${inbambase}
          BNAME=$outdir/${bname}.sbatchout
          DBASE=$(dirname "${BNAME}")
          [[ ! -d $DBASE ]] && echo "dbase $DBASE not found, exiting" && exit 1
  
          outf1=$outdir/${bname}.countTable.binsize_${binsize}.csv
          [[ -e $outf1 ]] && echo "outf1 $outf1 found, continuing" && continue
  
          cmd=". /nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh; conda activate scmo2022; bamToCountTable.py $inbam -sliding $stepsize --filterXA -minMQ $mapq -o $outf1 -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
          sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
      done
  done
# done

