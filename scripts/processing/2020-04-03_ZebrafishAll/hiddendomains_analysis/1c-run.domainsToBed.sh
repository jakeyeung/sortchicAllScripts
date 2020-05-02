#!/bin/sh
# Jake Yeung
# 1c-run.domainsToBed.sh
#  
# 2020-04-24

jmem='4G'
jtime='1:00:00'

jscript="/hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/domainsToBed.pl"
jscript2="/hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/domainsMergeToBed.pl"

[[ ! -e $jscript ]] && echo "$jscript not found, exiting" && exit 1
[[ ! -e $jscript2 ]] && echo "$jscript2 not found, exiting" && exit 1


chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/danRer11.chrom.sizes"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1
bwidth="1000"
[[ $bwidth != [0-9]* ]] && echo "Must be integer: $bwidth" && exit 1
minpost="0.6"

inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/hiddendomains_outputs.FromR"
[[ ! -d $inmain ]] && mkdir $inmain

for indir in `ls -d $inmain/PZ*`; do
    dname=$(basename $indir)
    inf="${indir}/${dname}_domains.txt"
    outf="${indir}/${dname}_vis.bed"
    outf2="${indir}/${dname}_analysis.bed"

    # [[ -e $ouf ]] && echo "$ouf found, continuing" && continue
    # [[ -e $outf2 ]] && echo "$outf2 found, continuing" && continue

    BNAME=$indir/HD_after_qsub_$dname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "$jscript -t -b $bwidth -g $chromsizes $inf > $outf; $jscript2 -b $bwidth -g $chromsizes -p $minpost $inf > $outf2" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N HD_after_$dname
    # echo "Debug exiting after one"
    # exit 0
done

