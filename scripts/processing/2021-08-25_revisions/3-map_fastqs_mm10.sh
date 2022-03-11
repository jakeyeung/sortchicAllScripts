#!/bin/sh
# Jake Yeung
# 3-map_fastqs_mm10.sh
#  
# 2021-08-25

bwabin="bwa"
indxbase="/hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/primary_assembly_NOMASK_ERCC92.fa"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/fastqs_from_bam"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/bam_remapped_mm10"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='32G'
jtime="24:00:00"
ncores=8

for f in `ls -d $indir/*.fastq`; do
  bname=$(basename $f)
  bname=${bname%.*}
  BNAME=${outdir}/${bname}.sbatch_log
  outf=${outdir}/${bname}.remapped.bam
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; $bwabin mem -t $ncores $indxbase $f | samtools view -Sb - > $outf"
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=$ncores --job-name=${bname} --wrap "$cmd"
done
