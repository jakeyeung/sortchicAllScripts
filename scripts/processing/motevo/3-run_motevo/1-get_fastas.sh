#!/bin/sh
# Jake Yeung
# 1-get_fastas.sh
# Get fastas onto scratch directory 
# process and run motevo 
# 2016-07-24

## BEGIN GET FASTA ## 
dailyORweekly="daily"
scratchmain="/scratch/el/$dailyORweekly/jyeung"
windowsf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions/regions_to_extract.merged.windows500.distfilt.50000.bed"
windowsbname=$(basename $windowsf)
windowsnoext=${windowsbname%.*}
windowsdirtmp="$scratchmain/bedregions_formotevo"
windowsftmp=$windowsdirtmp/$windowsbname
[[ ! -d $windowsdirtmp ]] && mkdir $windowsdirtmp
[[ ! -d $windowsdirtmp ]] && echo "$windowsdirtmp not found, exiting" && exit 1

# copy to scratch
cp --no-clobber $windowsf $windowsftmp

fastadirtmp=$scratchmain/fasta
[[ ! -d $fastadirtmp ]] && mkdir $fastadirtmp
[[ ! -d $fastadirtmp ]] && echo "$fastadirtmp not found, exiting" && exit 1
fastaref="/archive/epfl/upnae/jyeung/annotations/fasta_reference_ucsc/mm10.fa"
[[ ! -e $fastaref ]] && echo "$fastaref not found, exiting" && exit 1

fname=${windowed_name%%.*}
fastaftmp=$fastadirtmp/$windowsnoext.fa
[[ ! -e $fastaftmp ]] && bedtools2 getfasta -fi $fastaref -bed $windowsftmp -fo $fastaftmp && ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: getfasta failed" && exit 1
## END GET FASTA ## 

## BEGIN SPLIT BEDS ## 
[[ ! -e $fastaftmp ]] && echo "$fastaftmp not found, exiting" && exit 1

fastasplitdirtmp="$scratchmain/fastasplit"
[[ ! -d $fastasplitdirtmp ]] && mkdir $fastasplitdirtmp
[[ ! -d $fastasplitdirtmp ]] && echo "$fastasplitdirtmp not found, exiting" && exit 1

# use same number as in dhs_merged_tissue
n=60000  # must be EVEN number because fasta
rem=$(( $n % 2 ))
if [ $rem -eq 0 ]
then
  echo "Even number check: $n is OK!"
else
  echo "Even number check: $n is NOT OK... exiting!"
  exit 1
fi

# make into motevo format also
fastabase=$(basename $fastaftmp)
sed 's/>/>>mm10_/' $fastaftmp | split --lines=$n - $fastasplitdirtmp/$fastabase.

## END SPLIT BEDS ##

## BEGIN MOTEVO ## 

# wmdir="/archive/epfl/upnae/jyeung/databases/WMs/SwissRegulon"  # motevo runs in vitalit dont use archive
wmdir="/scratch/el/monthly/jyeung/databases/WMs/SwissRegulon"
[[ ! -d $wmdir ]] && echo "$wmdir not found, exiting" && exit 1

sitecountscript="/Home/jyeung/projects/shared_tissue_specificity/calculate_sitecount.py"
[[ ! -e $sitecountscript ]] && echo "$sitecountscript not found, exiting" && exit 1

# check split fastas exist
[[ ! -d $fastasplitdirtmp ]] && echo "$fastasplitdirtmp not found, exiting" && exit 1

paramdirname="param_files"

motevodirtmpbase="$scratchmain/motevo_outputs"
[[ ! -d $motevodirtmpbase ]] && mkdir $motevodirtmpbase
motevodirtmpsplit="$motevodirtmpbase/split"
[[ ! -d $motevodirtmpsplit ]] && mkdir $motevodirtmpsplit
[[ ! -d $motevodirtmpsplit ]] && echo "$motevodirtmpsplit not found, exiting" && exit 1

for infasta in `ls -d $fastasplitdirtmp/*.fa.*`
do
        bname=$(basename $infasta)
        chunk="${bname##*.}"
        motevodirtmpsplitchunk=$motevodirtmpsplit/$chunk
        paramsdir=$motevodirtmpsplitchunk/$paramdirname   

        [[ -d $motevodirtmpsplitchunk ]] && echo "chunk: $chunk found Skipping" && continue
        [[ ! -d $motevodirtmpsplitchunk ]] && mkdir $motevodirtmpsplitchunk
        [[ ! -d $paramsdir ]] && mkdir $paramsdir
        cd $paramsdir
        python $sitecountscript -w $wmdir -f $infasta -c 0 -o $motevodirtmpsplitchunk
		ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
done

## END MOTEVO ## 

# WAIT FOR MOTEVO JOBS TO FINISH # 
while [[ `bjobs | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# DOWNSTREAM HANDLING OF FILES # 

mergescript="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/mara_dhs/3-run_motevo_on_regions/merge_motevo_output.py"
[[ ! -e $mergescript ]] && echo "$mergescript not found, exiting" && exit 1

# check required directories
[[ ! -d $wmdir ]] && echo "$wmdir not found, exiting" && exit 1
[[ ! -d $motevodirtmpsplit ]] && echo "$motevodirtmpsplit not found, exiting" && exit 1

# Run mergescript
motevodirtmpmerged="$motevodirtmpbase/merged"
wmnames="/archive/epfl/upnae/jyeung/annotations/swissregulon/WMs.list"  # for checking we got all the motifs
if [[ -d $motevodirtmpmerged ]]
then
	echo "$motevodirtmpmerged exists, skipping"
else
	echo "Running merge script"
	mkdir $motevodirtmpmerged
	python $mergescript --wmnamesfile $wmnames $motevodirtmpsplit $motevodirtmpmerged
	ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
fi

# make bed
convertscript="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/mara_dhs/3-run_motevo_on_regions/merged_sites_to_bed.py"
[[ ! -e $convertscript ]] && echo "$convertscript not found, exiting" && exit 1

[[ ! -d $motevodirtmpmerged ]] && echo "$motevodirtmpmerged not found, exiting" && exit 1
motevodirtmpbed="$motevodirtmpbase/bed"

if [[ -d $motevodirtmpbed ]]
then
	echo "$motevodirtmpbed exists, skipping"
else
	echo "Running convert script"
	mkdir $motevodirtmpbed
	python $convertscript $motevodirtmpmerged $motevodirtmpbed --get_exact_region
fi


# Assign gene to bed
assignscript="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/merged_dhs/assign_nearest_gene_bed.py"
# output dir
motevodirtmpbedclosest=$motevodirtmpbase/closestbed_multiple_genes
# check input dirs and files
[[ ! -d $motevodirtmpbed ]] && echo "$motevodirtmpbed not found, exiting" && exit 1  # output dir
[[ ! -e $windowsftmp ]] && echo "$windowsftmp not found, exiting" && exit 1  # our ref bed
nohupdir="/scratch/el/monthly/jyeung/nohups/motevo_assign_nearestmulti"

[[ ! -e $assignscript ]] && echo "$assignscript not found, exiting" && exit 1
[[ ! -d $motevodirtmpbed ]] && echo "$motevodirtmpbed not found, exiting" && exit 1
[[ ! -d $motevodirtmpbedclosest ]] && mkdir $motevodirtmpbedclosest
[[ ! -d $nohupdir ]] && mkdir $nohupdir

for inbed in `ls -d $motevodirtmpbed/*.bed`; do
	base=$(basename $inbed)	
	bedout=$motevodirtmpbedclosest/${base%%.*}.closestmulti.bed
	[[ -e $bedout ]] && echo "$bedout found, continuing" && continue
	bsub -o $nohupdir/$base.out -e $nohupdir/$base.err -M 10000000 "python $assignscript $inbed $windowsftmp $bedout --has_motevo_id --save_pickle"
done

# wait for jobs
while [[ `bjobs | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# assign to each line
decompress="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/merged_dhs/convert_compressed_bed_to_long.py"
[[ ! -d $motevodirtmpbedclosest ]] && echo "$motevodirtmpbedclosest not found, exiting" && exit 1
motevodirtmpbedclosestlong=$motevodirtmpbedclosest/long_format
[[ ! -d $motevodirtmpbedclosestlong ]] && mkdir $motevodirtmpbedclosestlong
nohupdir="/scratch/el/monthly/jyeung/nohups/closestmulti_to_long_motifpeak"

[[ ! -e $decompress ]] && echo "$decompress not found, exiting" && exit 1
[[ ! -d $motevodirtmpbedclosest ]] && echo "$motevodirtmpbedclosest not found, exiting" && exit 1
[[ ! -d $motevodirtmpbedclosestlong ]] && mkdir $motevodirtmpbedclosestlong
[[ ! -d $nohupdir ]] && mkdir $nohupdir

for b in `ls -d $motevodirtmpbedclosest/*.bed`; do
        base=$(basename $b)
        basenoext=${base%.*}.long.bed
        bout=$motevodirtmpbedclosestlong/$basenoext
        [[ -e $bout ]] && echo "$bout found, continuing" && continue
        errf=$nohupdir/$basenoext.err
        outf=$nohupdir/$basenoext.out
        bsub -e $errf -o $outf -M 4000000 "python $decompress $b $bout --has_dist --gene_col_i 5"
        # python $decompress $b $bout --has_dist --gene_col_i 5
done

# wait for jobs
while [[ `bjobs | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# cat beds
merged_motevodirtmpbed="$motevodirtmpbed/merged_bed_closestbed_long"
[[ ! -d $merged_motevodirtmpbed ]] && mkdir $merged_motevodirtmpbed
[[ ! -d $merged_motevodirtmpbed ]] && echo "$merged_motevodirtmpbed not found, exiting" && exit 1

mergedname="motevo_merged.closest.long.bed"
echo "Running cat beds"
if [[ -e $merged_motevodirtmpbed/$mergedname ]]
then
	echo "Catt'd already"
else
	echo "Catting beds"
	# cat $motevodirtmpbedclosestlong/*.bed >> $merged_motevodirtmpbed/$mergedname
	# reorganize
	cat $motevodirtmpbedclosestlong/*.bed | sed  's/\;/\t/g' | sed 's/mm10_//g' | awk -F"\t" -v OFS="\t" ' { t = $5; $5 = $6; $6 = $7; $7 = $8; $8 = t; print; } ' | tr -d $"\r" >> $merged_motevodirtmpbed/$mergedname
fi


# save back to archive
merged_motevodir="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/motevo_merged_bed_long"
[[ ! -d $merged_motevodir ]] && mkdir $merged_motevodir
[[ ! -d $merged_motevodir ]] && echo "$merged_motevodir not found, exiting" && exit 1
cp --no-clobber $merged_motevodirtmpbed/$mergedname $merged_motevodir/$mergedname

