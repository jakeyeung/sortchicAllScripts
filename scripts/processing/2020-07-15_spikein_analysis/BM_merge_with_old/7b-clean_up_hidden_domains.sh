#!/bin/sh
# Jake Yeung
# 2-merge_hidden_domains_by_marks.sh
# Merge hidden domains outputs by marks  
# 2020-02-14

# # WRAP UP
# while [[ `qstat | wc -l` > 1 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

inbase0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"
[[ ! -d $inbase0 ]] && echo "$inbase0 not found, exiting" && exit 1

# all our hiddendomain directories are bad



# jsuffix="FromR.maxcount_40_60_R"
# inmain="${inbase0}/hiddendomains_outputs_minlength_2500.${jsuffix}"  # 2500 mostly but 5000 for basophils probably because small cluster
# [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
# 
# marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# 
# for mark in $marks; do
#     echo $mark
# 
#     outdir="${inbase0}/hd_merged.${mark}.${jsuffix}"
#     [[ ! -d $outdir ]] && mkdir $outdir
#     outname="merged.${mark}.cutoff_analysis.bed"
#     outnamemerged="merged.${mark}.cutoff_analysis.merged.bed"
#     outnamemergednochr="merged.${mark}.cutoff_analysis.merged.nochr.bed"
#     outnamemergedwithchr="merged.${mark}.cutoff_analysis.merged.withchr.bed"
#     outf=${outdir}/${outname}
#     outfmerged=${outdir}/${outnamemerged}
#     outfmergednochr=${outdir}/${outnamemergednochr}
#     outfmergedwithchr=${outdir}/${outnamemergedwithchr}
#     [[ -e $outf ]] && echo "$outf found, skipping" && continue
# 
#     for indir in $(ls -d ${inmain}/*BM_round1_round2_merged_${mark}*cutoff); do
#         echo $indir
#         bname=$(basename $indir)
#         bdir=$(echo $bname | cut --complement -d"." -f1)
#         # indir: maxcount_60x2500.BM_round1_round2_merged_H3K4me3_HSPCs.2500.cutoff/
#         # minlength=$(echo $bdir | cut -d"." -f2)
#         # echo "Minlength: $minlength"
#         # [[ $minlength != [0-9]* ]] && echo "Must be integer: $minlength" && exit 1
# 
#         echo $mark
# 
#         # bdir=$(basename $indir)
#         inbed="${indir}/${bdir}_analysis.bed"
#         echo $inbed
#         ## exit 0
#         [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
#         cat $inbed | awk -v bdir="${bdir}" 'BEGIN {OFS="\t"}; $4=bdir"_"$4' >> $outf
#     done
#     # sort in place
#     sort -k1,1 -k2,2n --output=$outf $outf
#     echo "Merging $outf to $outfmerged"
#     bedtools merge -i $outf -c 4 -o collapse > $outfmerged
#     # remove chr
#     sed 's/^chr//' $outfmerged | awk 'BEGIN {OFS="\t"}; $4="hd_"NR' > $outfmergednochr
#     # with chr, remove fourth coolumn used as input for anntoating peaks to gene later, which adds a fourth coulmn
#     cut -f1,2,3 $outfmerged > $outfmergedwithchr
#     # sed 's/^chr//' $outfmerged | awk 'BEGIN {OFS="\t"}; $4="hd_"NR' > $outfmergednochr
# done
