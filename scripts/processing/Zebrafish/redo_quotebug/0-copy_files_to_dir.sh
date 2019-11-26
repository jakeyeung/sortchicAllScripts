#!/bin/sh
# Jake Yeung
# 0-copy_files_to_dir.sh
# Copy ZF files from archive to working directory 
# 2019-11-22

ouds="OUD3913 OUD3909 OUD3910"

inmain="/hpc/archive/hub_oudenaarden/seqdata"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all"

find . -name "*PZ*ChIC*ZF*fastq.gz" 2>/dev/null | xargs -I {} cp {} $outdir
