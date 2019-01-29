#!/bin/sh
# Jake Yeung
# 4-run_matlong_to_sql.sh
# Run matlong to sql. Use cleaner code 
# 2016-08-03

sqlscript="/Home/jyeung/projects/make_sql_scripts/bed_to_sql.R"
# inf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/covbed.atacseq/dhs_signal_windows500.long.longforsql.mat"
inf="/archive/epfl/upnae/jyeung/sleep_deprivation/atacseq_signal/covbed.unmerged.atacseq/atacseq_signal_windows500.long.longforsql.mat"
infbase=$(basename $inf)
inftmp="/scratch/el/daily/jyeung/$infbase"
# outf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/covbed.atacseq/dhs_signal_windows500.long.motifindexed.sqlite"
outf="/archive/epfl/upnae/jyeung/sleep_deprivation/atacseq_signal/covbed.unmerged.atacseq/atacseq_signal_windows500.long.motifindexed.sqlite3"
outfbase=$(basename $outf)
outftmp="/scratch/el/daily/jyeung/$outfbase"
cnames="NA"
classes="character,integer,integer,character,integer,character,integer"
tblname="sleep_deprivation_atacseq_covbed_unmerged"
indx="chromo,start,end,gene,dist,sample"
nohupdir="/scratch/el/monthly/jyeung/nohups2"

cp --no-clobber $inf $inftmp
bsub -o $nohupdir/bed_to_sql.unmerged.out -e $nohupdir/bed_o_sql.unmerged.err -M 100000000 "Rscript $sqlscript $inftmp $outftmp $cnames $classes $tblname $indx"

# WRAP UP
while [[ `bjobs | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done
echo "Copying back to archive"
cp $outftmp $outf
