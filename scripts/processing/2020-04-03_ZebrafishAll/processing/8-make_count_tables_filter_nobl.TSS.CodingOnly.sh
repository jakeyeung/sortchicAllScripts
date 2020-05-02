#!/bin/sh
# Jake Yeung
# 6-make_count_and_RZ_tables.sh
#  
# 2019-12-21

# WRAP UP

# WRAP UP
wd="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks"
[[ ! -d $wd ]] && echo "$wd not found, exiting" && exit 1
cd $wd

# WRAP UP
while [[ `qstat | grep merge | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

inmain=$wd

jmem='32G'
jtime='6:00:00'

mapq=40
# binsize=10000  # try other binsizes?
binsize=10000  # try other binsizes?
stepsize=${binsize}
ps="/home/hub_oudenaarden/jyeung/projects/SingleCellMultiOmics.ForDev/singlecellmultiomics/bamProcessing/bamToCountTable.test.py"

inbed="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish/gene_tss.CodingOnly.winsize_${binsize}.species_drerio.nochr.asinteger.bed"

outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.TSS.CodingOnly.winsize_${binsize}"
[[ ! -d $outdir ]] && mkdir $outdir

for inbami in `ls -d $inmain/*.bam.bai`; do
    inbam=${inbami%.*}
    [[ ! -e $inbam ]] && echo "$inbam not found, continuing" && continue 
    bname=$(basename $inbam)
    bname=${bname%.*}

    BNAME=$outdir/${bname}.counttables.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    BNAME2=$outdir/${bname}.lhcounts.qsub
    DBASE2=$(dirname "${BNAME2}")
    [[ ! -d $DBASE2 ]] && echo "$DBASE2 not found, exiting" && exit 1

    outf=$outdir/${bname}.countTable.TSS.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam -sliding $stepsize --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname.countTable
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu  -N $bname.countTablesTSS
    # outf2=$outdir/${bname}.LHcounts.csv
    # [[ -e $outf2 ]] && echo "$outf2 found, continuing" && continue
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam -sampleTags SM -featureTags lh -o $outf2 --dedup --filterXA -minMQ $mapq" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME2}.out -e ${BNAME2}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname.RZcounts
done

