#!/bin/sh
# Jake Yeung
# 2-run_bed_to_sql.sh
# convert sh to sql 
# 2016-07-30

sqlscript="/Home/jyeung/projects/make_sql_scripts/bed_to_sql.R"
inf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/motevo_merged_bed_long/motevo_merged.closest.long.reorg.bed"
infbase=$(basename $inf)
inftmp="/scratch/el/daily/jyeung/$infbase"
outf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/motevo_merged_bed_long/motevo_merged.closest.long.motifindex.sqlite3"
outfbase=$(basename $outf)
outftmp="/scratch/el/daily/jyeung/$outfbase"
cnames="chromo,start,end,motif,sitecount,gene,dist,peak"
classes="character,integer,integer,factor,numeric,character,integer,character"
tblname="closestbed_multigene_swissregulon"
indx="chromo,start,end,gene,dist,motif"
nohupdir="/scratch/el/monthly/jyeung/nohups2"

cp --no-clobber $inf $inftmp
bsub -o $nohupdir/bed_to_sql.out -e $nohupdir/bed_o_sql.err -M 100000000 "Rscript $sqlscript $inftmp $outftmp $cnames $classes $tblname $indx"

# WRAP UP
while [[ `bjobs | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done
echo "Copying back to archive"
cp $outftmp $outf
