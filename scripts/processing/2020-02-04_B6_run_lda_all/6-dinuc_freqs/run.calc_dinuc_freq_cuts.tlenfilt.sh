#!/bin/sh
# Jake Yeung
# run.calc_dinuc_freq_cuts.sh
#  
# 2020-08-05

jmem='16G'
jtime='6:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/6-dinuc_freqs/calc_dinuc_freq_cuts.py"

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/bamFilteredByTlen"
refpath="/hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/primary_assembly_NOMASK_ERCC92.fa"

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/dinuc_freq/again2"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/dinuc_freq/tlen_147"
[[ ! -d $outdir ]] && mkdir $outdir

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits/H3K27me3-BM_AllMerged.HSPCs.sorted.cleaned.bam"
# [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

for inf in `ls -d $inmain/*.bam`; do
    bname=$(basename $inf)
    bname=${bname%.*}
    BNAME=$outdir/sbatchlog.${bname}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    outname=${bname}
    outf="${outdir}/${outname}"
    [[ -e $outf ]] && echo "$outf found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps -infile $inf -outprefix $outf -refpath $refpath -baseseq AT -mapqthres 40 -upstrm_extend 0 -downstrm_extend 200"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${bname}_%j.log --ntasks=1 --nodes=1 --job-name=$bname --wrap "$cmd"
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps -infile $inf -outprefix $outf -refpath $refpath -baseseq AT -mapqthres 40 -upstrm_extend 0 -downstrm_extend 200
    # exit 0
done
