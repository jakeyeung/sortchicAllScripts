#!/bin/sh
# Jake Yeung
# 1c-run.domainsToBed.sh
#  
# 2020-04-24

# WRAP UP
while [[ `qstat | grep R_HD_BM | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='8G'
jtime='2:00:00'

jscript="/hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/domainsToBed.pl"
jscript2="/hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/domainsMergeToBed.pl"

[[ ! -e $jscript ]] && echo "$jscript not found, exiting" && exit 1
[[ ! -e $jscript2 ]] && echo "$jscript2 not found, exiting" && exit 1


chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.filt.nochr.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1
# bwidth="1000"
bwidth="2500"
[[ $bwidth != [0-9]* ]] && echo "Must be integer: $bwidth" && exit 1
minpost="0.6"
# maxcount=60
# mincount=0

# inmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.first_round/all_bams"
# inmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.second_round/merged_bams"
inmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
# inmain="${inmain0}/hiddendomains_outputs/hiddendomains_outputs_minlength_${bwidth}.FromR.maxcount_${maxcount}_mincount_${mincount}_minprob_${minpost}"
# inmain="${inmain0}/hiddendomains_outputs_minlength_2500.FromR.maxcount_60_maxcount_40"
inmain="${inmain0}/hiddendomains_outputs_minlength_2500.FromR.maxcount_100_mincount_0_minprob_0.6"
# [[ ! -d $outmain ]] && mkdir $outmain

for indir in `ls -d $inmain/*cutoff`; do
    # bname=$(basename $indir)
    # dname=$(echo $bname | cut --complement -d"." -f1)
    dname=$(basename $indir)
    inf="${indir}/${dname}_domains.txt"
    outf="${indir}/${dname}_vis.bed"
    outf2="${indir}/${dname}_analysis.bed"

    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1


    # [[ -e $ouf ]] && echo "$ouf found, continuing" && continue
    # [[ -e $outf2 ]] && echo "$outf2 found, continuing" && continue

    BNAME=$indir/HD_after_qsub_$dname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "$jscript -t -b $bwidth -g $chromsizes $inf > $outf; $jscript2 -b $bwidth -g $chromsizes -p $minpost $inf > $outf2" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N HD_after_$dname
    # echo "Debug exiting after one"
    # exit 0
done

