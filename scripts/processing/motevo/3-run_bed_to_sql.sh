#!/bin/sh
# Jake Yeung
# 2-run_bed_to_sql.sh
# convert sh to sql 
# 2016-07-30

# sqlscript="/home/hub_oudenaarden/jyeung/projects/from_PhD/make_sql_scripts/bed_to_sql.R"
sqlscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/bed_to_sql.R"
inftmp="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"

outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/sql"
outname="motevo_merged.closest.long.peakindx.sqlite3"

[[ ! -d $outbase ]] && mkdir $outbase

outftmp=$outbase/$outname

cnames="chromo,start,end,motif,sitecount,gene,dist,peak"
classes="character,integer,integer,factor,numeric,character,integer,character"
tblname="closestbed_multigene_swissregulon_v2"
indx="chromo,gene,dist,motif,peak"
nohupdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/nohups_sql"

[[ ! -e $sqlscript ]] && echo "$sqlscript not found, exiting" && exit 1
[[ ! -e $inftmp ]] && echo "$inftmp not found, exiting" && exit 1

jmem='96G'
jtime='1:00:00'
BNAME=$nohupdir/bed_to_sql
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

# echo "Rscript $sqlscript $inftmp $outftmp $cnames $classes $tblname $indx" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
echo "Rscript $sqlscript $inftmp $outftmp $cnames $classes $tblname $indx"

