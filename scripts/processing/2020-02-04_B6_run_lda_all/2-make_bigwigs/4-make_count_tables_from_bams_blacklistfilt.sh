#!/bin/sh
# Jake Yeung
# 4-make_count_tables_from_bams_blacklistfilt.sh
#  
# 2020-05-28

jmem='16G'
jtime='5:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/SingleCellMultiOmics.ForDev/singlecellmultiomics/bamProcessing/bamToCountTable.WithBlacklist.py"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.count_tables_bl_filt"
mapq=40
stepsize=10000
bsize=10000
blfile="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"  # check chr names match bam



for inbam in `ls -d $indir/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}
    outf=${outdir}/${bname}.csv

    BNAME=$outdir/$bname.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $bsize -binTag DS --dedup -blacklist $blfile" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname.$bsize.countTable
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $bsize -binTag DS --dedup -blacklist $blfile"

done
