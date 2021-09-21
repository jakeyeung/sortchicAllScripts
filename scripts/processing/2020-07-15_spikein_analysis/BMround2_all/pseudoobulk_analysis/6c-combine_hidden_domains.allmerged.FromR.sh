#!/bin/sh
# Jake Yeung
# 6c-combine_hidden_domains.allmerged.FromR.sh
# Combine from different maxcount filters because it depends on the cluster 
# 2020-11-02

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_60_80"
[[ ! -d $outdir ]] && mkdir $outdir

indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_60_mincount_0_minprob_0.6"


indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_mincount_0_minprob_0.6"  # in case you have Vitterbi error

indir3="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_80_mincount_0_minprob_0.6"  # incase you're missing chromosomes

# Eryths have Vitterbi error when using too high maxcount, need to lower it
baddir2a="BM_round1_round2_merged_H3K9me3_Eryth.2500.cutoff"  # indir1 is bad, replcae with indir2

# missing chromosomes, increase maxcount
# baddir3a="BM_round1_round2_merged_H3K9me3_Eryth.2500.cutoff"  # indir1 is bad, replcae with indir3
baddir3b="BM_round1_round2_merged_H3K27me3_Bcells.2500.cutoff"
baddir3c="BM_round1_round2_merged_H3K4me1_Granulocytes.2500.cutoff"
baddir3d="BM_round1_round2_merged_H3K4me1_HSPCs.2500.cutoff"
# baddir3e="BM_round1_round2_merged_H3K9me3_Eryth.2500.cutoff"

# baddir3ALL="$baddir3a $baddir3b $baddir3c $baddir3d $baddir3e"
baddir3ALL="$baddir3b $baddir3c $baddir3d"

# cp -r $indir1 $outdir

cp -r $indir2/$baddir2a ${outdir}/maxcount_40.${baddir2a} 

for bd3 in $baddir3ALL; do
    cp -r ${indir3}/${bd3} ${outdir}/maxcount_80.${bd3}
done

for indir in `ls -d $indir1/*cutoff`; do
    dname=$(basename $indir)
    echo $dname
    [[ -e $outdir/maxcount_40.$dname ]] && echo "$outdir/maxcount_40.$dname found, continuing" && continue
    [[ -e $outdir/maxcount_80.$dname ]] && echo "$outdir/maxcount_80.$dname found, continuing" && continue
    # for bd3 in $baddir3ALL; do
    #     [[ $dname == $bd3 ]] && echo "Matches not found, continuing for $dname" && continue
    # done
    echo $indir
    cp --no-clobber -r $indir $outdir/maxcount_60.${dname}
done


