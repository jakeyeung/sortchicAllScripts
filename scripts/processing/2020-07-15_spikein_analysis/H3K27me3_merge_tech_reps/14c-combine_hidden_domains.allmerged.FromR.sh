#!/bin/sh
# Jake Yeung
# 6c-combine_hidden_domains.allmerged.FromR.sh
# Combine from different maxcount filters because it depends on the cluster 
# 2020-11-02

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/hiddendomains_outputs"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outdir="${inmain}/hiddendomains_outputs_minlength_500.mincount_-10.FromR.maxcount_10_40_60"
[[ ! -d $outdir ]] && mkdir $outdir

indir0="${inmain}/hiddendomains_outputs_minlength_500.FromR.maxcount_40_mincount_-10_minprob_0.6"

# indir1="${inmain}/hiddendomains_outputs_minlength_500.FromR.maxcount_60_mincount_-10_minprob_0.6"  # for Eryths handles missing chromo
indir2="${inmain}/hiddendomains_outputs_minlength_500.FromR.maxcount_10_mincount_-10_minprob_0.6"  # handles Basos and NKs missing chromosmoes 16

# baddir1a="PZ-BM-H3K27me3-Eryths-merged.500.cutoff"  # replace with maxcount 60
# baddir1ALL=${baddir1a}

baddir2a="PZ-BM-H3K27me3-Basophils-merged.500.cutoff"
baddir2b="PZ-BM-H3K27me3-NKs-merged.500.cutoff"
baddir2ALL="${baddir2a} ${baddir2b}"

# for bd1 in $baddir1ALL; do
#     cp -r ${indir1}/${bd1} ${outdir}/maxcount_60x500.${bd1}
# done

for bd2 in $baddir2ALL; do
    cp -r ${indir2}/${bd2} ${outdir}/maxcount_10x500.${bd2}
done

for indir in `ls -d $indir0/*cutoff`; do
    dname=$(basename $indir)
    echo $dname
    # [[ -e $outdir/maxcount_60x500.$dname ]] && echo "$outdir/maxcount_60x500.$dname found, continuing" && continue
    [[ -e $outdir/maxcount_10x500.$dname ]] && echo "$outdir/maxcount_10x500.$dname found, continuing" && continue
    echo $indir
    cp --no-clobber -r $indir $outdir/maxcount_40x500.${dname}
done

