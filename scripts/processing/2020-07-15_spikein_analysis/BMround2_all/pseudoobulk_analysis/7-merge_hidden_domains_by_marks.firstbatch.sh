#!/bin/sh
# Jake Yeung
# 2-merge_hidden_domains_by_marks.sh
# Merge hidden domains outputs by marks  
# 2020-02-14

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

# inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output"

# inmain=${inbase}/hd_clusters.${mark}.minlength_${minlength}
# inmain="${inbase0}/hiddendomains_outputs.FromR" 
inbase0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.first_round/all_bams"
[[ ! -d $inbase0 ]] && echo "$inbase0 not found, exiting" && exit 1
inmain="${inbase0}/hiddendomains_outputs.FromR"

# inbase="${inbase0}/hiddendomains_outputs.FromR"
# # inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.first_round/all_bams/hiddendomains_outputs.FromR"

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
minlength=1000

for mark in $marks; do
    echo $mark
    outdir=${inbase0}/hd_merged.${mark}.minlength_${minlength}
    [[ ! -d $outdir ]] && mkdir $outdir
    outname="merged.${mark}.minlength_${minlength}.cutoff_analysis.bed"
    outnamemerged="merged.${mark}.minlength_${minlength}.cutoff_analysis.merged.bed"
    outnamemergednochr="merged.${mark}.minlength_${minlength}.cutoff_analysis.merged.nochr.bed"
    outnamemergedwithchr="merged.${mark}.minlength_${minlength}.cutoff_analysis.merged.withchr.bed"
    outf=${outdir}/${outname}
    outfmerged=${outdir}/${outnamemerged}
    outfmergednochr=${outdir}/${outnamemergednochr}
    outfmergedwithchr=${outdir}/${outnamemergedwithchr}
    [[ -e $outf ]] && echo "$outf found, skipping" && continue
    for indir in $(ls -d ${inmain}/${mark}*cutoff); do
        bdir=$(basename $indir)
        echo "Merging $bdir"
        inbed="${indir}/${bdir}_analysis.bed"
        [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
        cat $inbed | awk -v bdir="${bdir}" 'BEGIN {OFS="\t"}; $4=bdir"_"$4' >> $outf
    done
    # sort in place
    sort -k1,1 -k2,2n --output=$outf $outf
    echo "Merging $outf to $outfmerged"
    bedtools merge -i $outf -c 4 -o collapse > $outfmerged
    # remove chr
    sed 's/^chr//' $outfmerged | awk 'BEGIN {OFS="\t"}; $4="hd_"NR' > $outfmergednochr
    # with chr, remove fourth coolumn used as input for anntoating peaks to gene later, which adds a fourth coulmn
    cut -f1,2,3 $outfmerged > $outfmergedwithchr
    # sed 's/^chr//' $outfmerged | awk 'BEGIN {OFS="\t"}; $4="hd_"NR' > $outfmergednochr
done
