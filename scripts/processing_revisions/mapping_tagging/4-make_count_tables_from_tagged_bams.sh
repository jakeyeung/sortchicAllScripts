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

jmem='16G'
jtime='24:00:00'

# inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/differentiation/round1/tagged_bams"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams"
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams"

mapq=40
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nospikeins.nochromo.bed"

binsizes="10000 50000"

# outbase="/hpc/hub_oudenaarden/jyeung/data/dblchic/differentiation/round1/count_tables"
outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/count_tables"
[[ ! -d $outbase ]] && mkdir $outbase

experi="BM"
jmarks="k4me1-k9me3 k27me3 k9me3"
# jmark="H3K36me3"

for jmark in $jmarks; do 
  echo $jmark
  for binsize in $binsizes; do
      stepsize=$binsize
  
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
  
          cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate SCMO_2021; bamToCountTable.py $inbam -sliding $stepsize --filterXA -minMQ $mapq -o $outf1 -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
          sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
      done
  done


done

