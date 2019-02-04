#!/bin/sh
# Jake Yeung
# 2019-02-02_send_trackhub_to_websync.sh
# Send trackhub to websync 
# 2019-02-02

indir="/Users/yeung/data/trackhubs/H3K4me1_peaks"
outdir="websync@upnaesrv1.epfl.ch:/data/web/sites/motevo_from_peaks/"

rsync -avrL --copy-links $indir $outdir
