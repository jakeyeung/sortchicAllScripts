#!/bin/sh
# Jake Yeung
# 5-run.split_bam_by_barcodes.sh
# Run split bam 
# 2018-12-15

jmem='4G'
jtime='0:30:00'

jscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/split_bam_by_barcodes.py"
# bname="PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14"
inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95"
# bcmain="/home/hub_oudenaarden/jyeung/projects/scChiC/outputs_R/barcode_summaries"
bcmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/barcode_summaries_BM_build95/summarized"

[[ ! -e $jscript ]] && echo "$jscript not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_build95.withchr"
[[ ! -d $outmain ]] && mkdir $outmain
for bcf in $(ls -d $bcmain/barcode_summary.PZ*.thres.0.chip.*.txt); do
    # echo $bcf
    bname=$(echo $bcf | awk '{split($0, a, "."); print(a[2])}')
    inbam=$inmain/$bname/$bname.filtered.sorted.bam
    outdir="$outmain/$bname"
    [[ -d $outdir ]] && echo "$outdir exists, skipping..." && continue
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    [[ ! -e $bcf ]] && echo "$bcf not found, exiting" && exit 1
    [[ ! -d $outdir ]] && mkdir $outdir
    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $jscript --add_chr_prefix $inbam $bcf $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $jscript --add_chr_prefix $inbam $bcf $outdir
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $jscript --add_chr_prefix $inbam $bcf $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
    # exit 0
    # echo "python $jscript $inbam $bcf $outdir"
   #  | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
done

# bcf="/home/hub_oudenaarden/jyeung/projects/scChiC/outputs_R/barcode_summaries/barcode_summary.PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.thres.10000.chip.H3K27me3.txt"
