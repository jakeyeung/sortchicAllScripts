#!/bin/sh
# Jake Yeung
# link_all_tagged_bams.sh
#  
# 2020-10-08

indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams"
indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5230_BM/tagged_bams"
indir3="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM/mouse/tagged_bams"
indir4="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM/mouse/tagged_bams"
indir5="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5109_BM_twoplates/tagged_bams"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links"

for indir in ${indir1} ${indir2} ${indir3} ${indir4} ${indir5}; do
    echo $indir
    for f in `ls -d $indir/*.bam*`; do
        fbase=$(basename $f)
        outf=${outdir}/${fbase}
        [[ -e $outf ]] && echo "$outf found, continuing" && continue
        ln -s $f ${outf}
    done
done

