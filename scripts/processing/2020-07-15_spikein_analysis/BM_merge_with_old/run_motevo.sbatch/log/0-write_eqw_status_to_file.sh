#!/bin/sh
# Jake Yeung
# 0-write_eqw_status_to_file.sh
#  
# 2020-02-18

qstat | awk 'NR > 2' | for J in `cut -d" " -f 1`; do qstat -j $J; done > /home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/4-run_motevo/log/eqw_jobs_need_rerun.txt
