#!/bin/sh
# Jake Yeung
# 1-copy_bed_split_extract_cat_beds.sh
# Combine first 3 steps of /home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/motevo_dhs_scripts_clean/2-extract_dhs_signal
# Copy to VITALIT and run jobs.
# 2016-07-22

nohupmain="/scratch/el/monthly/jyeung/nohups/extract_regions.unmerged.longqueue"
[[ ! -d $nohupmain ]] && mkdir $nohupmain

bedf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions/regions_to_extract.merged.windows500.distfilt.50000.bed"
bedbname=$(basename $bedf)

bamdir="/archive/epfl/upnae/jyeung/sleep_deprivation/data_from_charlotte.atacseq/bams_unmerged.atacseq"
# covdir="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/covbed.atacseq"
covdir="/archive/epfl/upnae/jyeung/sleep_deprivation/atacseq_signal/covbed.unmerged.atacseq"
[[ ! -d $covdir ]] && mkdir $covdir  # final output directory

# create scratch folders
scratchdailyORweekly="/scratch/el/weekly"
[[ ! -d $scratchdailyORweekly ]] && echo "$scratchdaily must be daily or weekly" && exit 1
scratchmain=$scratchdailyORweekly/jyeung
[[ ! -d $scratchmain ]] && mkdir $scratchmain
inbname=$(basename $bamdir)
outbname=$(basename $covdir)
bamdirtmpbed=$scratchmain/bedregions
bamdirtmp=$scratchmain/$inbname
covdirtmp=$scratchmain/$outbname
[[ ! -d $covdirtmp ]] && mkdir -p $covdirtmp
[[ ! -d $bamdirtmp ]] && mkdir -p $bamdirtmp
[[ ! -d $bamdirtmpbed ]] && mkdir $bamdirtmpbed

# copy gene regions
inbedtmp=$bamdirtmpbed/$bedbname
cp --no-clobber $bedf $inbedtmp

# split into about 100 files
[[ ! -e $inbedtmp ]] && echo "$inbedtmp not found, exiting" && exit 1
splitbeddirtmp=$bamdirtmpbed/"split_inbed"
[[ ! -d $splitbeddirtmp ]] && mkdir $splitbeddirtmp
[[ ! -d $splitbeddirtmp ]] && echo "$splitbeddirtmp not found, exiting" && exit 1
outsplitprefix=$bedbname"."
# split --lines=46854 $inbedtmp $splitbeddirtmp/$outsplitprefix
split --lines=110000 $inbedtmp $splitbeddirtmp/$outsplitprefix  # enough files to launch all into "long" queue simultaneously 
# split --lines=20000 $inbedtmp $splitbeddirtmp/$outsplitprefix  # lots of files so that you finish in less than 24 hrs

# # copy bam to vitalit and get bams_str
[[ ! -d $bamdirtmp ]] && echo "$bamdirtmp not found, exiting" && exit 1
for bamin in `ls -d $bamdir/*.bam*`; do
	cp --no-clobber $bamin $bamdirtmp
done

# get string for list of bam files
bamstr=`find $bamdirtmp -type f -name "*.bam" | sort`
# bamstrordered=`find $bamdirtmp -type f -name "*.bam" | sort`
# echo $bamstrordered
# write to output dir to link columns to colnames in output file
echo $bamstr | tr " " "\n" >  $covdir/metadata.txt

bamstr=`echo $bamstr | tr "\n" " "`
# bams_str=`cat $bams | tr '\n' ' '`

# 
outprefix="atacseq_signal_windows500.mat"
for bedin in `ls -d $splitbeddirtmp/*.bed*`; do
        base=$(basename $bedin)
        ext="${base##*.}"
        outname=$outprefix.$ext
        outfile=$covdirtmp/$outname
		# run only if outfile does not exist
		[[ ! -e $outfile ]] && bsub -q long -n 1 -o $nohupmain/$ext.out -e $nohupmain/$ext.err "bedtools2 multicov -bams $bamstr -bed $bedin > $outfile"
		# echo "bedtools2 multicov -bams $bamstr -bed $bedin > $outfile"
done


# WRAP UP
while [[ `bjobs | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done
echo "Copying back to archive"

# cat beds
[[ ! -e $covdir/$outprefix ]] && cat $covdirtmp/*.mat* >> $covdir/$outprefix
# rsync -vaP $covdirtmp/ $covdir
