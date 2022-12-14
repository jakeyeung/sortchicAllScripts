#!/bin/sh
# Jake Yeung
# 1-get_fastas.sh
# Get fastas onto scratch directory 
# process and run motevo 
# 2019-03-22

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py2

jmark="H3K9me3"
scratchmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/${jmark}"

## BEGIN GET FASTA ## 

windowsf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output/BM_${jmark}_merged.1000.cutoff/BM_${jmark}_merged.1000.cutoff_analysis.blacklistfilt.annot.bed"
[[ ! -e $windowsf ]] && echo "$windowsf not found, exiting" && exit 1
windowsbname=$(basename $windowsf)
windowsnoext=${windowsbname%.*}

fastaref="/hpc/hub_oudenaarden/jyeung/data/databases/fasta/mm10.fa"
[[ ! -e $fastaref ]] && echo "$fastaref not found, exiting" && exit 1

fastadirtmp=$scratchmain/fasta
[[ ! -d $fastadirtmp ]] && mkdir $fastadirtmp
fname=${windowsbname%%.*}
fastaftmp=$fastadirtmp/$windowsnoext.fa
# go to UPPER
[[ ! -e $fastaftmp ]] && bedtools2 getfasta -fi $fastaref -bed $windowsf | awk 'NR % 2 { print } !(NR % 2) {print toupper($0)}' > $fastaftmp && ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: getfasta failed" && exit 1
# echo "bedtools2 getfasta -fi $fastaref -bed $windowsf -fo $fastaftmp"
# [[ ! -e $fastaftmp ]] && bedtools2 getfasta -fi $fastaref -bed $windowsf -fo $fastaftmp && ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: getfasta failed" && exit 1

## END GET FASTA ## 

## BEGIN SPLIT BEDS ## 
[[ ! -e $fastaftmp ]] && echo "$fastaftmp not found, exiting" && exit 1

fastasplitdirtmp="$scratchmain/fastasplit"
[[ ! -d $fastasplitdirtmp ]] && mkdir $fastasplitdirtmp
[[ ! -d $fastasplitdirtmp ]] && echo "$fastasplitdirtmp not found, exiting" && exit 1

# use same number as in dhs_merged_tissue
# n=60000 # works for H3K4me1 and H3K4me3
n=10000  # must be EVEN number because fasta
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
wmdir="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm10_weight_matrices_v2_split"
[[ ! -d $wmdir ]] && echo "$wmdir not found, exiting" && exit 1

# sitecountscript="/Home/jyeung/projects/shared_tissue_specificity/calculate_sitecount.py"
sitecountscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/calculate_sitecount.py"
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
        python $sitecountscript -w $wmdir -f $infasta -c 0 -o $motevodirtmpsplitchunk -d $paramsdir
		ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
done

# echo "Created motevo inputs, but not yet run"

## END MOTEVO ## 

# WAIT FOR MOTEVO JOBS TO FINISH # 
while [[ `qstat | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# DOWNSTREAM HANDLING OF FILES # 

mergescript="/home/hub_oudenaarden/jyeung/projects/from_PhD/tissue-specificity-shellscripts/mara_dhs/3-run_motevo_on_regions/merge_motevo_output.py"
[[ ! -e $mergescript ]] && echo "$mergescript not found, exiting" && exit 1

# check required directories
[[ ! -d $wmdir ]] && echo "$wmdir not found, exiting" && exit 1
[[ ! -d $motevodirtmpsplit ]] && echo "$motevodirtmpsplit not found, exiting" && exit 1

# Run mergescript
motevodirtmpmerged="$motevodirtmpbase/merged"
# wmnames="/archive/epfl/upnae/jyeung/annotations/swissregulon/WMs.list"  # for checking we got all the motifs
wmnames="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm10_v2_WMs.list"
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
# convertscript="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/mara_dhs/3-run_motevo_on_regions/merged_sites_to_bed.py"
convertscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/merged_sites_to_bed.py"
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
# assignscript="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/merged_dhs/assign_nearest_gene_bed.py"
assignscript="/home/hub_oudenaarden/jyeung/projects/from_PhD/tissue-specificity-shellscripts/merged_dhs/assign_nearest_gene_bed.py"
# output dir
motevodirtmpbedclosest=$motevodirtmpbase/closestbed_multiple_genes
# check input dirs and files
[[ ! -d $motevodirtmpbed ]] && echo "$motevodirtmpbed not found, exiting" && exit 1  # output dir
[[ ! -e $windowsf ]] && echo "$windowsf not found, exiting" && exit 1  # our ref bed
# nohupdir="/scratch/el/monthly/jyeung/nohups/motevo_assign_nearestmulti"
nohupdir="$scratchmain/nohups_motevo_assign_nearestmulti"

[[ ! -e $assignscript ]] && echo "$assignscript not found, exiting" && exit 1
[[ ! -d $motevodirtmpbed ]] && echo "$motevodirtmpbed not found, exiting" && exit 1
[[ ! -d $motevodirtmpbedclosest ]] && mkdir $motevodirtmpbedclosest
[[ ! -d $nohupdir ]] && mkdir $nohupdir

jmem="10G"
jtime="8:00:00"
for inbed in `ls -d $motevodirtmpbed/*.bed`; do
	base=$(basename $inbed)	
	bedout=$motevodirtmpbedclosest/${base%%.*}.closestmulti.bed
	[[ -e $bedout ]] && echo "$bedout found, continuing" && continue
	# bsub -o $nohupdir/$base.out -e $nohupdir/$base.err -M 10000000 "python $assignscript $inbed $windowsf $bedout --has_motevo_id --save_pickle"
    BNAME=$nohupdir/$base
	echo "python $assignscript $inbed $windowsf $bedout --has_motevo_id --save_pickle" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
done

# wait for jobs
while [[ `qstat | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# assign to each line
# decompress="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/merged_dhs/convert_compressed_bed_to_long.py"
decompress="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/convert_compressed_bed_to_long.py"
[[ ! -d $motevodirtmpbedclosest ]] && echo "$motevodirtmpbedclosest not found, exiting" && exit 1
motevodirtmpbedclosestlong=$motevodirtmpbedclosest/long_format
[[ ! -d $motevodirtmpbedclosestlong ]] && mkdir $motevodirtmpbedclosestlong
# nohupdir="/scratch/el/monthly/jyeung/nohups/closestmulti_to_long_motifpeak"
nohupdir="$scratchmain/nohups_closestmulti_to_long_motifpeak"

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
        # bsub -e $errf -o $outf -M 4000000 "python $decompress $b $bout --has_dist --gene_col_i 5"
        BNAME=$nohupdir/$base
        echo "python $decompress $b $bout --has_dist --gene_col_i 5" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
done

# wait for jobs
while [[ `qstat | wc -l` > 1 ]]; do
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

