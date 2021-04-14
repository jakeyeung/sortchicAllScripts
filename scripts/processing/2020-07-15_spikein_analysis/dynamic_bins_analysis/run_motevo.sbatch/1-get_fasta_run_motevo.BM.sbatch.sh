#!/bin/sh
# Jake Yeung
# 1-get_fastas.sh
# Get fastas onto scratch directory 
# process and run motevo 
# 2019-03-22

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py2

# redefined below
jmem="8G"
jtime="2:00:00"
# jmem="10G"
# jtime="8:00:00"

# jmark="H3K9me3"
# jmark="H3K9me3"
jmark="H3K27me3"

prefix="cluster"
suffix="BM_dynamic"
suffix2="bins"
minlength=1000
scratchmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_${prefix}_${suffix}_${suffix2}_${minlength}/${jmark}"

# comment uot this if scratchmain already exists and you want to rerun because script failed halfway
# getfasta, split, and running motevo will be skipped if outptput dirs are found
[[ -d $scratchmain ]] && echo "$scratchmain found, exiting" && exit 1 
[[ ! -d $scratchmain ]] && mkdir -p $scratchmain

## BEGIN GET FASTA ## 
# windowsf="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/hiddendomains_outputs.FromR/hd_merged.${jmark}.minlength_1000/merged.${jmark}.minlength_${minlength}.cutoff_analysis.merged.withchr.annotated.bed"
# windowsf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_H3K9me3_dynamic_bins/coords_H3K9me3_dynamic_bins.bed"
# windowsf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_H3K9me3_dynamic_bins/coords_H3K9me3_dynamic_bins.noname.2021-01-28.annotated.bed"
windowsf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_top_6085_four_marks_dynamic_bins/top_6085.${jmark}.2021-02-17.annotated.bed"

# windowsf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/hd_merged.H3K4me3.minlength_1000/merged.H3K4me3.minlength_1000.cutoff_analysis.merged.withchr.annotated.test.bed"
[[ ! -e $windowsf ]] && echo "$windowsf not found, exiting" && exit 1
windowsbname=$(basename $windowsf)
windowsnoext=${windowsbname%.*}

fastaref="/hpc/hub_oudenaarden/jyeung/data/databases/fasta/mm10.fa"
# fastaref="/hpc/hub_oudenaarden/jyeung/data/databases/fasta/Danio_rerio.GRCz11.dna.primary_assembly.withchr.fa"
[[ ! -e $fastaref ]] && echo "$fastaref not found, exiting" && exit 1

fastadirtmp=$scratchmain/fasta
[[ ! -d $fastadirtmp ]] && mkdir $fastadirtmp
fname=${windowsbname%%.*}
fastaftmp=$fastadirtmp/$windowsnoext.fa
# go to UPPER
# [[ ! -e $fastaftmp ]] && bedtools2 getfasta -fi $fastaref -bed $windowsf | awk 'NR % 2 { print } !(NR % 2) {print toupper($0)}' > $fastaftmp && ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: getfasta failed" && exit 1

[[ ! -e $fastaftmp ]] && bedtools getfasta -fi $fastaref -bed $windowsf | awk 'NR % 2 { print } !(NR % 2) {print toupper($0)}' > $fastaftmp && ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: getfasta failed" && exit 1

# exit 0
## END GET FASTA ## 

## BEGIN SPLIT BEDS ## 
[[ ! -e $fastaftmp ]] && echo "$fastaftmp not found, exiting" && exit 1

fastasplitdirtmp="$scratchmain/fastasplit"

if [ ! -d $fastasplitdirtmp ]
then
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
else
	echo "fastasplitdirtmp $fastasplitdirtmp exists, skipping the split"
fi


## END SPLIT BEDS ##

## BEGIN MOTEVO ## 

wmdir="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm10_weight_matrices_v2_split"
[[ ! -d $wmdir ]] && echo "$wmdir not found, exiting" && exit 1

# sitecountscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/calculate_sitecount.py"
sitecountscript="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/calculate_sitecount.sbatch.py"
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

        [[ -d $paramsdir ]] && echo "$paramsdir found, continuing" && continue

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
while [[ `squeue -u jyeung | grep MOTEVO | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# DOWNSTREAM HANDLING OF FILES # 

# mergescript="/home/hub_oudenaarden/jyeung/projects/from_PhD/tissue-specificity-shellscripts/mara_dhs/3-run_motevo_on_regions/merge_motevo_output.py"
mergescript="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/merge_motevo_output.py"
[[ ! -e $mergescript ]] && echo "$mergescript not found, exiting" && exit 1

# check required directories
[[ ! -d $wmdir ]] && echo "$wmdir not found, exiting" && exit 1
[[ ! -d $motevodirtmpsplit ]] && echo "$motevodirtmpsplit not found, exiting" && exit 1

# Run mergescript
motevodirtmpmerged="$motevodirtmpbase/merged"
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
# convertscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/merged_sites_to_bed.py"
convertscript="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/merged_sites_to_bed.py"
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
# assignscript="/home/hub_oudenaarden/jyeung/projects/from_PhD/tissue-specificity-shellscripts/merged_dhs/assign_nearest_gene_bed.py"
assignscript="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/assign_nearest_gene_bed.py"
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
jtime="2:00:00"

for inbed in `ls -d $motevodirtmpbed/*.bed`; do
	base=$(basename $inbed)	
	bedout=$motevodirtmpbedclosest/${base%%.*}.closestmulti.bed
	[[ -e $bedout ]] && echo "$bedout found, continuing" && continue
	# bsub -o $nohupdir/$base.out -e $nohupdir/$base.err -M 10000000 "python $assignscript $inbed $windowsf $bedout --has_motevo_id --save_pickle"
    BNAME=$nohupdir/$base
	# echo "python $assignscript $inbed $windowsf $bedout --has_motevo_id --save_pickle" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -N MOTEVO
	cmd="python $assignscript $inbed $windowsf $bedout --has_motevo_id --save_pickle"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=MOTEVO --wrap "$cmd" 
done

# wait for jobs
while [[ `squeue -u jyeung | grep MOTEVO |  wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# assign to each line
# decompress="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/merged_dhs/convert_compressed_bed_to_long.py"
# decompress="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/convert_compressed_bed_to_long.py"
decompress="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/convert_compressed_bed_to_long.py"
[[ ! -d $motevodirtmpbedclosest ]] && echo "$motevodirtmpbedclosest not found, exiting" && exit 1
motevodirtmpbedclosestlong=$motevodirtmpbedclosest/long_format
[[ ! -d $motevodirtmpbedclosestlong ]] && mkdir $motevodirtmpbedclosestlong
# nohupdir="/scratch/el/monthly/jyeung/nohups/closestmulti_to_long_motifpeak"
nohupdir="$scratchmain/nohups_closestmulti_to_long_motifpeak"

[[ ! -e $decompress ]] && echo "$decompress not found, exiting" && exit 1
[[ ! -d $motevodirtmpbedclosest ]] && echo "$motevodirtmpbedclosest not found, exiting" && exit 1
[[ ! -d $motevodirtmpbedclosestlong ]] && mkdir $motevodirtmpbedclosestlong
[[ ! -d $nohupdir ]] && mkdir $nohupdir

jmem="1G"
jtime="1:00:00"
for b in `ls -d $motevodirtmpbedclosest/*.bed`; do
        base=$(basename $b)
        basenoext=${base%.*}.long.bed
        bout=$motevodirtmpbedclosestlong/$basenoext
        [[ -e $bout ]] && echo "$bout found, continuing" && continue
        errf=$nohupdir/$basenoext.err
        outf=$nohupdir/$basenoext.out
        # bsub -e $errf -o $outf -M 4000000 "python $decompress $b $bout --has_dist --gene_col_i 5"
        BNAME=$nohupdir/$base
        # echo "python $decompress $b $bout --has_dist --gene_col_i 5" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -N MOTEVO
        cmd="python $decompress $b $bout --has_dist --gene_col_i 5"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=MOTEVO --wrap "$cmd"
done

# wait for jobs
while [[ `squeue -u jyeung | grep MOTEVO | wc -l` > 0 ]]; do
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

echo "Done catting beds... script finished"
