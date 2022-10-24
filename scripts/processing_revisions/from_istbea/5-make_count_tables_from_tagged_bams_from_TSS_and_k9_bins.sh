#!/bin/sh
# Jake Yeung
# 5-make_count_tables_from_tagged_bams_from_TSS_and_k9_bins.sh
#  
# 2022-07-18


# # WRAP UP
# while [[ `squeue -u jyeung | wc -l` > 1 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='16G'
jtime='12:00:00'

inbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/tagged_bams"
[[ ! -d $inbase ]] && echo "$inbase not found, exiting" && exit 1

mapq=40
bl="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/mm10.blacklist.copy.sorted.merged.nospikeins.nochromo.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

binsizes="10000"

outbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables_fromTSS_and_k9_bins"
[[ ! -d $outbase ]] && mkdir $outbase

experi="BM"
# jmarks="k4me1-k9me3 k27me3 k9me3"
jmarks="k4me1 k4me3 k27me3 k9me3"
# jmark="H3K36me3"

for jmark in $jmarks; do 
  echo $jmark
  for binsize in $binsizes; do
      stepsize=$binsize
	  # bedfile="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/refseq/MmRefseqTss.chromorenamed.${binsize}.again.nochromo.sorted.bed"
	  bedfile="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/refseq/MmRefseqTss.chromorenamed.${binsize}.with_k9_bins.2022-07-18.txt"
	  [[ ! -e $bedfile ]] && echo "$bedfile not found, exiting" && exit 1

      outmain="${outbase}/counts_tables_${binsize}"
      [[ ! -d $outmain ]] && mkdir $outmain
      outdir="${outmain}/${experi}_${jmark}"
      [[ ! -d $outdir ]] && mkdir $outdir

     inmain="${inbase}/${experi}_${jmark}"
     [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
     for inbam in `ls -d $inmain/*.bam`; do
          inbambase=$(basename $inbam)
          inbambase=${inbambase%.*}
          bname=${inbambase}
          BNAME=$outdir/${bname}.sbatchout
          DBASE=$(dirname "${BNAME}")
          [[ ! -d $DBASE ]] && echo "dbase $DBASE not found, exiting" && exit 1
  
          outf1=$outdir/${bname}.countTable.binsize_${binsize}.csv
          [[ -e $outf1 ]] && echo "outf1 $outf1 found, continuing" && continue
  
          cmd=". /nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh; conda activate scmo2022; bamToCountTable.py $inbam --filterXA -minMQ $mapq -o $outf1 -sampleTags SM -joinedFeatureTags reference_name -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels -bed ${bedfile}"
          sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
      done
  done


done

