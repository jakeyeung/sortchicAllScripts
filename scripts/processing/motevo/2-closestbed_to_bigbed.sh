#!/bin/sh
# Jake Yeung
# 7-closestbed_to_bigbed.sh
# Convert closestbed to bigbed for visualization on genome browser 
# 2016-01-06
# 2019-02-01  # redo as postdoc

# qsub params
jmem='8G'
jtime='1:00:00'

bedtobigbed="/hpc/hub_oudenaarden/jyeung/software/ucsc_utils/bedToBigBed"

# maindir="/scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs_cleaned"
maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs"

closestbeddir="$maindir/bed"  # or bed_stranded
outbeddir="$maindir/bigbeds"
chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.txt"

# tmpdir="/scratch/el/daily/jyeung/tmp"
# tmpdir="/tmp/jake"
tmpdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/tmp"
mkdir -p $tmpdir
nohupmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/nohups_closestbed_to_bed"

[[ ! -e $bedtobigbed ]] && echo "$bedtobigbed not found, exiting" && exit 1
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1
[[ ! -d $closestbeddir ]] && echo "$closestbeddir not found, exiting" && exit 1
[[ ! -d $outbeddir ]] && mkdir $outbeddir
[[ ! -d $nohupmain ]] && mkdir $nohupmain

for b in `ls -d $closestbeddir/*.bed`; do
	bbase=$(basename $b)
	bbase_notext=${bbase%%.*}
	outf=$outbeddir/$bbase_notext.bb
	# round label to 2 decimals, print only rows with scores greater than 0.5
	jcmd="awk -F $'\t' 'BEGIN {OFS = FS} {lab=sprintf(\"%.2f\", \$5); if(\$5 >= 0.5) print \$1, \$2, \$3, lab, int(\$5*1000)}' $b | LC_COLLATE=C sort -k1,1 -k2,2n > $tmpdir/$bbase; $bedtobigbed $tmpdir/$bbase $chromsizes $outf"
	# echo $jcmd
	# break

    BNAME=$nohupmain/$bbase_notext
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # maybe faster without qsub
    echo "$jcmd" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
	# bsub -o $nohupmain/$bbase_notext.out -e $nohupmain/$bbase_notext.err "$jcmd"
	# [[ -e $outf ]] && echo "Deleting tmp" && rm $tmpdir/$bbase
done

