#!/bin/sh
# Jake Yeung
# 1-liftover_dr7_to_dr10.sh
#  
# 2020-04-15

export PATH=/hpc/hub_oudenaarden/jyeung/software/ucsc_utils:$PATH

cf="/hpc/hub_oudenaarden/jyeung/data/databases/chainfiles/danRer7ToDanRer10.over.chain.gz"
inbed="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer7_CNEs_PEGASUS.forliftover.bed"
outbed="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer10_CNEs_PEGASUS.forliftover.bed"
unmapped="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer10_CNEs_PEGASUS.forliftover.unmapped"

liftOver $inbed $cf $outbed $unmapped


