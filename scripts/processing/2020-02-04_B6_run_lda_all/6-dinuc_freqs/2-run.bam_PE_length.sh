#!/bin/sh
# Jake Yeung
# 2-run.bam_PE_length.sh
#  
# 2020-07-02

jmem='16G'
jtime='4:00:00'

# inmain="/hpc/hub_oudenaarden/jyeung/data/chimeseq/bams/from_AvO/OUD4506_tagged"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17"
outdir="${inmain}/bamPEFragmentSize"
[[ ! -d $outdir ]] && mkdir $outdir

BNAME=${outdir}/PEfragsize_sbatchlog
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
for inbam in `ls -d $inmain/*.bam`; do 

    bname=$(basename $inbam)
    bname=${bname%.*}

    outpng=${outdir}/${bname}_PEfragsize.png
    outtxt=${outdir}/${bname}_PEfragsize.txt

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamPEFragmentSize --bamfiles ${inbam} --histogram ${outpng} --plotFileFormat png --samplesLabel ${bname} --blackListFileName $bl --outRawFragmentLengths $outtxt"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=PEfragsize --wrap "$cmd"

done

