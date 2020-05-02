#!/bin/sh
# Jake Yeung
# 2-liftover_dr10_to_dr11.sh
#  
# 2020-04-15

export PATH=/hpc/hub_oudenaarden/jyeung/software/ucsc_utils:$PATH

cf="/hpc/hub_oudenaarden/jyeung/data/databases/chainfiles/danRer10ToDanRer11.over.chain.gz"
[[ ! -e $cf ]] && echo "$cf not found, exiting" && exit 1

inbed="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer10_CNEs_PEGASUS.forliftover.bed"
outbed="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer11_CNEs_PEGASUS.forliftover.bed"
unmapped="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer11_CNEs_PEGASUS.forliftover.unmapped"

liftOver $inbed $cf $outbed $unmapped

