#!/bin/sh
# Jake Yeung
# 6-make_RZ_counts.sh
#  
# 2019-09-04

# sleep 3600

jmem='8G'
jtime='2:00:00'
# inbase="/hpc/hub_oudenaarden/jyeung/data/dblchic/differentiation/round1/tagged_bams"
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams"
[[ ! -d $inbase ]] && echo "$inbase not found, exiting" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/count_tables/LH_RZ_counts"
[[ ! -d $outdir ]] && mkdir $outdir

experi="BM"
jmarks="k4me1-k9me3 k27me3 k9me3"

for jmark in $jmarks; do
  inmain=${inbase}/${experi}_${jmark}
  echo $jmark
  for inbam in `ls -d $inmain/*${jmark}*.bam`; do
      bname=$(basename $inbam)
      bname=${bname%.*}
  
      outf=$outdir/${bname}.LH_counts.csv
      [[ -e $outf ]] && echo "$outf found, continuing" && continue
  
      BNAME=$outdir/${bname}.LHcounts.qsub
      DBASE=$(dirname "${BNAME}")
      [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
      cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate SCMO_2021; bamToCountTable.py $inbam -sampleTags SM -featureTags lh -o $outf --dedup --filterXA -minMQ 40 --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
      echo $cmd
      sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
  done
  # wait
done


