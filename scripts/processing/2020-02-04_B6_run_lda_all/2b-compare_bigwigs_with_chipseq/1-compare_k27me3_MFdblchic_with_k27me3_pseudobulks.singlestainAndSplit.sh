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

# jsuffix="unnorm"
# inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/public_chipseq_H3K27me3_celltypes/${jsuffix}"
inbase="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/bigwig_outputs/merged_bams.bigwigs_by_cluster"
[[ ! -d $inbase ]] && echo "$inbase not found, exiting" && exit 1
blist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
[[ ! -e $blist ]] && echo "$blist not found, exiting" && exit 1
chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster"
ref="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.BM_giladi.keeptop_1000.unsorted.bed"

# mark="H3K27me3"
mark="K27m3"
markalt="H3K27me3"

cd $inbase  # so ls doesn't give full name
# all_BM_K27m3_200119.topic11.sorted.bw => all_BM_K27m3_200119.topic11
ctypes=$(ls -d *.bw | sed 's/.bw//g' | cut -d"." -f1,2 | sort | uniq | tr '\n' ' ') 

outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.compare_with_MFdblchic.${mark}.again"
[[ ! -d $outbase ]] && mkdir $outbase

n=0
maxjobs=2
for ctype in $ctypes; do
    ctype=${ctype}.
    echo $ctype
    infbase=$inbase/${ctype}*.bw
    infcompares=$(ls -d $chicbwdir/*${markalt}*.bw)
    infcompares=$(echo $infcompares | tr '\n' ' ')
    outf="$outbase/${mark}_${ctype}_comparison.npz"
    [[ -e $outf ]] && echo "$outf  found, continuing" && continue
    outftxt="$outbase/${mark}_${ctype}_comparison.txt"
    [[ -e $outftxt ]] && echo "$outftxt found, exiting" && exit 1

    BNAME=$outbase/${mark}_${ctype}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo $infbase
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; multiBigwigSummary BED-file --BED $ref -b $infbase $infcompares -o $outf -bl $blist --outRawCounts $outftxt"
    # exit 0
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; multiBigwigSummary BED-file --BED $ref -b $infbase $infcompares -o $outf -bl $blist --outRawCounts $outftxt" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N bwsum_${ctype}
done

echo "Done"
