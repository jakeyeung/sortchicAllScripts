#!/bin/sh
# Jake Yeung
# 6c-combine_hidden_domains.allmerged.FromR.sh
# Combine from different maxcount filters because it depends on the cluster 
# 2020-11-02

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outdir="${inmain}/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_60_R"
[[ ! -d $outdir ]] && mkdir $outdir

indir0="${inmain}/hiddendomains_outputs_minlength_2500"
indir1="${inmain}/hiddendomains_outputs_minlength_2500.FromR.maxcount_60_mincount_0_minprob_0.6"
indir2="${inmain}/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_mincount_0_minprob_0.6"  # in case you have Vitterbi error
indir3="${inmain}/hiddendomains_outputs_minlength_5000.FromR.maxcount_60_mincount_0_minprob_0.6"  # K4me1 basophils vitterbi error needs larger minlength?
indir4="${inmain}/hiddendomains_outputs_minlength_2500.FromR.maxcount_20_mincount_0_minprob_0.6"  # K4me3_pDCs missing chromsoome if maxcount > 40, maxcount=20 highest without missing chromosome. Too low may give Vitterbi error?

baddir0a="BM_round1_round2_merged_H3K27me3_DCs.2500.cutoff"
baddir0c="BM_round1_round2_merged_H3K4me1_Granulocytes.2500.cutoff"
baddir0d="BM_round1_round2_merged_H3K4me1_HSPCs.2500.cutoff"
baddir0e="BM_round1_round2_merged_H3K27me3_Bcells.2500.cutoff"
baddir0ALL="$baddir0a $baddir0b $baddir0c $baddir0d"

# Eryths have Vitterbi error when using too high maxcount, need to lower it
baddir2a="BM_round1_round2_merged_H3K9me3_Eryth.2500.cutoff"  # indir1 is bad, replcae with indir2

# missing chromosomes, increase maxcount
# baddir3a="BM_round1_round2_merged_H3K9me3_Eryth.2500.cutoff"  # indir1 is bad, replcae with indir3
# baddir3b="BM_round1_round2_merged_H3K27me3_Bcells.2500.cutoff"
# baddir3c="BM_round1_round2_merged_H3K4me1_Granulocytes.2500.cutoff"
# baddir3d="BM_round1_round2_merged_H3K4me1_HSPCs.2500.cutoff"
# baddir3e="BM_round1_round2_merged_H3K9me3_Eryth.2500.cutoff"


# baddir3ALL="$baddir3a $baddir3b $baddir3c $baddir3d $baddir3e"
# baddir3ALL="$baddir3b $baddir3c $baddir3d"
baddir3b="BM_round1_round2_merged_H3K4me1_Basophils.5000.cutoff"
baddir3ALL="$baddir3b"

baddir4a="BM_round1_round2_merged_H3K4me3_pDCs.2500.cutoff"

# cp -r $indir1 $outdir

for bd0 in $baddir0ALL; do
    cp -r ${indir0}/${bd0} ${outdir}/maxcount_Rx2500.${bd0}
done

cp -r $indir2/$baddir2a ${outdir}/maxcount_40x2500.${baddir2a} 

for bd3 in $baddir3ALL; do
    cp -r ${indir3}/${bd3} ${outdir}/maxcount_60x5000.${bd3}
done

cp -r $indir4/$baddir4a ${outdir}/maxcount_20x2500.${baddir4a} 

for indir in `ls -d $indir1/*cutoff`; do
    dname=$(basename $indir)
    echo $dname
    [[ -e $outdir/maxcount_Rx2500.$dname ]] && echo "$outdir/maxcount_Rx2500.$dname found, continuing" && continue
    [[ -e $outdir/maxcount_20x2500.$dname ]] && echo "$outdir/maxcount_40x2500.$dname found, continuing" && continue
    [[ -e $outdir/maxcount_40x2500.$dname ]] && echo "$outdir/maxcount_40x2500.$dname found, continuing" && continue
    [[ -e $outdir/maxcount_60x5000.$dname ]] && echo "$outdir/maxcount_60x5000.$dname found, continuing" && continue
    # for bd3 in $baddir3ALL; do
    #     [[ $dname == $bd3 ]] && echo "Matches not found, continuing for $dname" && continue
    # done
    echo $indir
    cp --no-clobber -r $indir $outdir/maxcount_60x2500.${dname}
done

echo "Manually delete:"
echo "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_60_R/maxcount_60x2500.BM_round1_round2_merged_H3K4me1_Basophils.2500.cutoff"
