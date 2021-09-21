#!/bin/sh
# Jake Yeung
# 3-run_plot_topics_after_LDA_WKM_and_BM.sh
#  
# 2020-06-23

# for WKM only 

jmem='16G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/plot_LDA_topics.R"
Rscript $rs $inf $outf
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisTSS_BM_WKM_dists"  
refbed="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish.SelectedTSSfromWKM/gene_tss.SelectedFromWKM.winsize_50000.species_drerio.bed"

for indir in `ls -d $inmain/*Zebrafish*2020-06-22*`; do
    inf=$(ls -d $indir/*.Robj)
    echo $inf
    bname=$(basename $inf)
    bname=${bname%.*}

    BNAME=${indir}/${bname}.sbatch.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${indir}/${bname}.topics_plot.pdf
    # cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf -keeptop 150 -inftss $refbed"
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf -keeptop 150 --TopicsOnly"
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=$bname --wrap "$cmd"
done
