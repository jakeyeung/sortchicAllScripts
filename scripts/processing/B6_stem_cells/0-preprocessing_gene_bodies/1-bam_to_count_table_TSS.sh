#!/bin/sh
# Jake Yeung
# 1-bam_to_count_table_TSS.sh
#  
# 2019-06-23

mapq=60
# prefixs="H3K36me3 K36me3_K4me1 K36me3_K9me3 K9me3"
# inmain="/hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs"
# inmain="/hpc/hub_oudenaarden/hvinas/chic/merged_VAN3478_VAN3479_VAN3480"
# [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

inbed="/hpc/hub_oudenaarden/jyeung/data/databases/genebodies/nochr/gene_start_end.2019-08-17.pfilt_0.8.up_2000.down_2000.bed"
[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
suffix=$(basename $inbed)
suffix=${suffix%.*}

# outdir="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/count_tables_genebodies"
# [[ ! -d $outdir ]] && mkdir $outdir

inmain="/hpc/hub_oudenaarden/jyeung/raw_data_from_sequencer/AVO508/links_merged/raw_demultiplexed"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3700/countTables_genebodies"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='32G'
jtime='12:00:00'
# for prefix in $prefixs; do
for indir in $(ls -d ${inmain}/PZ*); do
    echo $indir
    fnamenoext=$(basename $indir)

    BNAME=$outdir/${fnamenoext}_gene_start_end
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1


    inbam="${indir}/tagged/bwaMapped.sorted.bam"
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    outf=$outdir/${fnamenoext}.filtered.mapq_${mapq}.${suffix}.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate singlecellmultiomicsenv2; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N $fnamenoext
    # echo "bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed"
done
