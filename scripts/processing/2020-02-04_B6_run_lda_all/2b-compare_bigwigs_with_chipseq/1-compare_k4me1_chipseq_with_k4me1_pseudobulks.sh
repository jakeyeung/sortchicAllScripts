#!/bin/sh
# Jake Yeung
# 5-compare_bigwig_with_published.sh
#  
# 2020-02-27
# compare K4me1 with published data
# compare bigwigs for known celltype versus all k4me1 at TSS
# 

jmem='16G'
jtime='1:00:00'

inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science/renamed"
blist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
[[ ! -e $blist ]] && echo "$blist not found, exiting" && exit 1
# chicbwdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/bigwig_outputs/merged_bams.bigwigs_by_cluster"
chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster"
ref="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.BM_giladi.keeptop_1000.unsorted.bed"

mark="H3K4me1"
markalt="H3K4me1"

cd $inbase  # so ls doesn't give full name
ctypes=$(ls -d *.bw | sed 's/.bw//g' | cut -d"_" -f2 | tr '\n' ' ')

echo $ctypes

outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.compare_with_chipseq.${mark}"
[[ ! -d $outbase ]] && mkdir $outbase

n=0
maxjobs=2
for ctype in $ctypes; do
    echo $ctype
    infbase=$inbase/${mark}_${ctype}.bw
    [[ ! -e $infbase ]] && echo "$infbase not found, exiting" && exit 1
    infcompares=$(ls -d $chicbwdir/*${markalt}*.bw)
    infcompares=$(echo $infcompares | tr '\n' ' ')
    outf="$outbase/${mark}_${ctype}_comparison.npz"
    [[ -e $outf ]] && echo "$outf  found, continuing" && continue
    outftxt="$outbase/${mark}_${ctype}_comparison.txt"
    [[ -e $outftxt ]] && echo "$outftxt found, exiting" && exit 1

    BNAME=$outbase/${mark}_${ctype}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; multiBigwigSummary BED-file --BED $ref -b $infbase $infcompares -o $outf -bl $blist --outRawCounts $outftxt" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N bwsum_${ctype}
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; multiBigwigSummary BED-file --bed $ref -b $infbase $infcompares -o $outf -bl $blist --outRawCounts $outftxt"
    # exit 0

    # if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    #         # define maxjobs and n using maxjobsn skeleton
    #     wait # wait until all have finished (not optimal, but most times good enough)
    #     echo $n wait
    # fi
done

echo "Done"
