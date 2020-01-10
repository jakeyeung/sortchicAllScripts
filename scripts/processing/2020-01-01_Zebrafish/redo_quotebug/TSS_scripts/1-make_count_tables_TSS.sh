#!/bin/sh
# Jake Yeung
# 8-make_count_tables_TSS.sh
#  
# 2019-11-24

jmem='32G'
jtime='4:00:00'
mapq=40
jdists="20000 50000"

for jdist in $jdists; do
    inbed="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss.winsize_${jdist}.species_drerio.nochr.bed"
    [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
    suffix=$(basename $inbed)
    suffix=${suffix%.*}

    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all/raw_demultiplexed"
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all/countTables_geneTSS"
    [[ ! -d $outdir ]] && mkdir $outdir

    for indir in $(ls -d ${inmain}/PZ*); do
        echo $indir
        fnamenoext=$(basename $indir)

        BNAME=$outdir/${jdist}_${fnamenoext}_gene_start_end
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        inbam="${indir}/tagged/$fnamenoext.bwaMapped.tagged.bam"
        [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1

        outf=$outdir/${fnamenoext}.filtered.mapq_${mapq}.${suffix}.csv
        [[ -e $outf ]] && echo "$outf found, continuing" && continue

        echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N ${jdist}.$fnamenoext
    done
done
