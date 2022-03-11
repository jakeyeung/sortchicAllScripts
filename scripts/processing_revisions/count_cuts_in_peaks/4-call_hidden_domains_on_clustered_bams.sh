#!/bin/sh
# Jake Yeung
# 1-call_hidden_domains_on_clustered_bams.sh 
# 2019-04-15
# 4-call_hidden_domains_on_clustered_bams.sh
# 2022-01-10

jmem='4G'
jtime='1:00:00'

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.hg38.filt.nochr.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

# inbase="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.imputevarfilt.lessstringent.mapq_40"jA
# inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_G1filt_split_by_G1filt_NoChrPrefix"
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams/K562/merged_bams/bams_split_by_good_cells_NoChrPrefix"

outmain="${inbase}/hiddendomains_outputs"
[[ ! -d $outmain ]] && mkdir $outmain

# marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
marks="k4me1 k4me3 k27me3 k9me3"

maxcounts=10
minlength=1000  # seems to work well at 1k?

for mark in $marks; do
    for inf in $(ls -d $inbase/*${mark}.GoodCells.sorted.bam); do
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        bname=$(basename $inf)
        bname=${bname%.*}

        bname=$bname.$minlength.cutoff
        [[ ! -d $outmain ]] && mkdir $outmain
        outdir=$outmain/$bname

        # [[ -d $outdir ]] && echo "$outdir found, skipping $mark" && continue

        [[ ! -d $outmain ]] && mkdir $outmain
        [[ ! -d $outdir ]] && mkdir $outdir
        BNAME=$outdir/$bname.log

        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname"
        sbatch --time=${jtime} --mem-per-cpu=${jmem} --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${mark}_${bname} --wrap "$cmd"
    done
done

