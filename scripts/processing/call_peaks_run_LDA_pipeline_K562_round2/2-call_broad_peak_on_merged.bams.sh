#!/bin/sh
# Jake Yeung
# 2-call_broad_peak_on_merged.bams.sh
# Call broad peak on merged bams 
# 2019-01-08: Run on K562 on across cell cycles

jmem='10G'
jtime='3:00:00'

marks="H3K27me3 H3K4me1"

genome="hs"

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged/BM_H3K4me1_merged.bam"

for mark in $marks; do
    # echo $mark
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_K562_round2/K562_${mark}_merged.bam"
    bname=$(basename $inf)
    bname=${bname%.*}

    cutoff=0.5
    minlength=1000  # seems to work well at 1k?
    bname=$bname.$cutoff.$minlength.cutoff
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/macs2_output_K562_round2"
    outdir=$outmain/$bname

    [[ -d $outdir ]] && echo "$outdir found, skipping $mark" && continue

    [[ ! -d $outmain ]] && mkdir $outmain
    [[ ! -d $outdir ]] && mkdir $outdir
    BNAME=$outdir/$bname.log

    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py2; macs2 callpeak --cutoff-analysis -t $inf --broad -g $genome --pvalue $cutoff --broad-cutoff $cutoff --min-length $minlength -n $bname --outdir $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
    # echo "Run macs2 for $mark"
done

