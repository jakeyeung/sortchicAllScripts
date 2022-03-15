#!/bin/sh
# Jake Yeung
# 4b-sort_lone_bam_file.sh
#  
# 2022-01-24

jmem='64G'
jtime='24:00:00'
ncores=8

inbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data"
inmain="${inbase}/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
tmpdir=$inmain/tmpdirmulti
[[ ! -d $tmpdir ]] && mkdir $tmpdir

wd="${inbase}/workdir"
[[ ! -d $wd ]] && mkdir $wd
cd $wd

for indir in `ls -d ${inmain}/PZ-sortChIC-BM-SL4-k27me3-2`; do
    bname=$(basename $indir)
    inbam=$indir/bwaMapped.bam
	outdir=$indir
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    BNAME=$outdir/$bname.sort_only.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    sortedbam=$indir/bwaMapped.sorted.bam
    outbamtagged=$outdir/${bname}.sorted.tagged.bam
    [[ -e $outbamtagged ]] && echo "$outbamtagged found, continuing" && continue
    cmd=". /nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh; conda activate scmo2022; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam"
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

