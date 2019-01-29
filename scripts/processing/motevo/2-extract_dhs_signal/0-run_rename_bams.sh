#!/bin/sh
# Jake Yeung
# 0-run_rename_bams.sh
# Run rename bams but dont merge 
# 2016-08-05

rscript="/Home/jyeung/projects/sleep_deprivation/atacseq_analysis/2-extract_dhs_signal/unmerge_rename_bams_to_long_format.R"

inbam="/archive/epfl/upnae/jyeung/sleep_deprivation/data_from_charlotte.atacseq/bams_renamed_matching_rnaseq"
outdir="/archive/epfl/upnae/jyeung/sleep_deprivation/data_from_charlotte.atacseq/bams_unmerged.atacseq"

Rscript $rscript $inbam $outdir
