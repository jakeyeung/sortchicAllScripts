#!/bin/sh
# Jake Yeung
# 5-run.calculate_cut_positions.sh
#  
# 2020-08-06

jmem='16G'
jtime='2:00:00'

radius=2000

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/5-measure_cut_distances.from_bedfile/calculate_cut_positions.py"
beddir="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/bams_split_by_clusters"

strands="neg pos"

for strand in $strands; do
    bedfile="${beddir}/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS.txt.${strand}.chromorenamed.bed.gz"
    # bedfile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.${strand}.bed.gz"

    # indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits"

    for indir in `ls -d $inmain/K562*`; do

        outdir="${indir}/cut_positions_TSS_rad_${radius}_${strand}.K562_spikeins"
        [[ ! -d $outdir ]] && mkdir $outdir

        for inf in `ls -d $indir/*.bam`; do
            bname=$(basename $inf)
            bname=${bname%.*}

            BNAME=${outdir}/${bname}.sbatchoutput
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

            outprefix=${outdir}/${bname}.TSS_cuts.radius_${radius}.${strand}

            # cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -outprefix $outprefix -bedfile $bedfile -radius $radius"
            cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps -infile $inf -outprefix $outprefix -bedfile $bedfile -radius $radius"
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
            # exit 0
        done

    done

done

