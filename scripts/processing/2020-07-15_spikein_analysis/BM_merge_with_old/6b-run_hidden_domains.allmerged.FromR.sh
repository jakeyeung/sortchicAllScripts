#!/bin/sh
# Jake Yeung
# 1-call_hidden_domains_on_clustered_bams.sh 
# 2019-04-15

# # WRAP UP
# while [[ `squeue -u jyeung | grep "BM_round" wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='2G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/hiddendomains_analysis/run_hidden_domains.R"

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.filt.nochr.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

# inmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.first_round/all_bams"
# inmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.second_round/merged_bams"

# maxcount="40"
minprob="0.6"
# minlengths="2500 5000 10000"
minlength="5000"
# maxcounts="100 120 140"
# maxcounts="160 180 200 220"
# maxcounts="60 80 100"
maxcounts="10 20 40"
mincount="0"

chromonames="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y"
# for minlength in $minlengths; do
for maxcount in $maxcounts; do

    # inmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
    inmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"
    inmain="${inmain0}/hiddendomains_outputs_minlength_${minlength}"
    [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

    outmain="${inmain0}/hiddendomains_outputs_minlength_${minlength}.FromR.maxcount_${maxcount}_mincount_${mincount}_minprob_${minprob}"
    [[ ! -d $outmain ]] && mkdir $outmain

    # marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"


    for indir in `ls -d $inmain/*cutoff`; do
        dname=$(basename $indir)
        outdir=${outmain}/${dname}
        [[ ! -d $outdir ]] && mkdir $outdir

        inf="${indir}/${dname}_treatment_bins.txt"
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        outf="${outdir}/${dname}_domains.txt"
        [[ -e $outf ]] && echo "$outf found, continuing" && continue
        BNAME=$outdir/$dname.log

        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -maxreadcount $maxcount -minreadcount $mincount -minprob $minprob -chromonames ${chromonames}" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N "R_HD_$dname"
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -maxreadcount $maxcount -minreadcount $mincount -minprob $minprob -chromonames ${chromonames}"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=HD_downstream_${bname} --wrap "$cmd"
    done


done
# minlength="2500"

