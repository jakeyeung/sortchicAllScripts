#!/bin/sh
# Jake Yeung
# 0-annotate_peaks_with_gene_names.sh
#  
# 2020-02-14

jmem='32G'
jtime='6:00:00'

# inscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/annotate_bed_to_gene_and_distance.R"
inscript="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/annotate_bed_to_gene_and_distance.R"  # new location, loads scchicFuncs

[[ ! -e $inscript ]] && echo "$inscript not found, exiting" && exit 1
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
jlength=1000

for mark in $marks; do
    indir="${inmain}/hd_merged.${mark}.minlength_${jlength}"
    outdir=${indir}
    jbase="merged.${mark}.minlength_${jlength}.cutoff_analysis"
    inname=${jbase}.bed
    outname=${jbase}.annotated.bed

    inf=${indir}/${inname}
    outf=${outdir}/${outname}

    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=${outdir}/${mark}.qsub_annot
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $inscript $inf $outf"
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $inscript $inf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N Anno_${mark}
done

