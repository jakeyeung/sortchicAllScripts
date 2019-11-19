#!/bin/sh
# Jake Yeung
# 1-bam_to_count_table_TSS.sh
#  
# 2019-06-23

jmem='32G'
jtime='4:00:00'
mapq=40
# jdist="100000"
# oud="oud3910"
# oud="oud3909"
jdists="20000 50000 100000"
ouds="oud3909 oud3910"

for jdist in $jdists; do
    for oud in $ouds; do
        inbed="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss.winsize_${jdist}.species_drerio.nochr.bed"
        [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
        suffix=$(basename $inbed)
        suffix=${suffix%.*}

        inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/${oud}_HD_0/raw_demultiplexed"
        outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/${oud}_HD_0/countTables_geneTSS"
        [[ ! -d $outdir ]] && mkdir $outdir

        for indir in $(ls -d ${inmain}/PZ*); do
            echo $indir
            fnamenoext=$(basename $indir)

            BNAME=$outdir/${jdist}_${fnamenoext}_gene_start_end
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

            inbam="${indir}/tagged/bwaMapped.sorted.bam"
            [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
            outf=$outdir/${fnamenoext}.filtered.mapq_${mapq}.${suffix}.csv
            [[ -e $outf ]] && echo "$outf found, continuing" && continue
            echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N $fnamenoext
        done
    done
done

