#!/bin/sh
# Jake Yeung
# 4d-filter_chromos.sh
# Filter chromos ERCC and stuff because it messes up bedgraph to bigwig 
# 2019-04-02
# 
# https://unix.stackexchange.com/questions/293684/basic-grep-awk-help-extracting-all-lines-containing-a-list-of-terms-from-one-f

csizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromlist.txt"
[[ ! -e $csizes ]] && echo "$csizes not found, exiting" && exit 1

suffix="build95_B6"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bedgraphs_${suffix}/log1p"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bedgraphs_${suffix}/log1p/chr_filt"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

[[ ! -d $outdir ]] && mkdir $outdir

for inf in `ls -d $indir/*.bedgraph`; do
    echo $inf
    fbase=$(basename $inf)
    outf=$outdir/$fbase
    [[ -e $outf ]] && echo "$outf not found, continuing" && continue
    grep -w -F -f $csizes $inf > $outf
done


