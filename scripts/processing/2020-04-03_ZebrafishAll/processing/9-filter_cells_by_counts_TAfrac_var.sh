#!/bin/sh
# Jake Yeung
# 0c-filter_cells_by_counts_TAfrac_var.sh
#  
# 2020-04-03

jmem='16G'
jtime='1:00:00'

winsize=50000

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/filter_cells_by_totalcuts_TAfrac_intrachromovar.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/bonemarrow/RZdat.all.deeper_EtOH_definitive_NoTcells"
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_${winsize}"
indirmat=${indir}
# indirmat="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/bonemarrow/countTables.all.deeper_EtOH_definitive_NoTcells.50000_50000"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $indirmat ]] && echo "$indirmat not found, exiting" && exit 1

m1="H3K4me1"
m2="H3K4me3"
m3="H3K27me3"
m4="H3K9me3"

infmat1=${indir}/"PZ-ChIC-ZF_${m1}_2020-04-07.countTable.csv"
infmat2=${indir}/"PZ-ChIC-ZF_${m2}_2020-04-07.countTable.csv"
infmat3=${indir}/"PZ-ChIC-ZF_${m3}_2020-04-07.countTable.csv"
infmat4=${indir}/"PZ-ChIC-ZF_${m4}_2020-04-07.countTable.csv"

inf1=${indir}/"PZ-ChIC-ZF_${m1}_2020-04-07.LHcounts.csv"
inf2=${indir}/"PZ-ChIC-ZF_${m2}_2020-04-07.LHcounts.csv"
inf3=${indir}/"PZ-ChIC-ZF_${m3}_2020-04-07.LHcounts.csv"
inf4=${indir}/"PZ-ChIC-ZF_${m4}_2020-04-07.LHcounts.csv"

chromos="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25"

ccutoff="1000 500 1000 1000"
tcutoff="0.5"
countsmax="30000 10000 50000 50000"
varmin="0.1 0.1 0.1 0.075"

# stringent
# countsmax="1000000 1000000 1000000 1000000"
# varmin="0.12 0.2 0.15 0.1"

varmax="1 1 1 1"
outdir="${indir}.chromo_filtered_by_counts_TAfrac_var.lessstringent.2020-04-14"
[[ ! -d $outdir ]] && mkdir $outdir

BNAME=$outdir/filter_cells_qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -countcutoffmin $ccutoff -countcutoffmax $countsmax -varcutoffmin $varmin -varcutoffmax $varmax -TAcutoff $tcutoff  -infilerz $inf1 $inf2 $inf3 $inf4 -infilecounts $infmat1 $infmat2 $infmat3 $infmat4 -names $m1 $m2 $m3 $m4 -chromoskeep $chromos -outdir $outdir" --overwrite

# echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -countcutoffmin $ccutoff -countcutoffmax $countsmax -varcutoffmin $varmin -varcutoffmax $varmax -TAcutoff $tcutoff  -infilerz $inf1 $inf2 $inf3 $inf4 -infilecounts $infmat1 $infmat2 $infmat3 $infmat4 -names $m1 $m2 $m3 $m4 -chromoskeep $chromos -outdir $outdir" --overwrite | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N filter_cells_ZF
