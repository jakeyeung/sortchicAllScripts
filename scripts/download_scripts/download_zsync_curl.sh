#!/bin/sh
# Jake Yeung
# download_zsync_curl.sh
#  
# 2019-01-14

outdir="/hpc/hub_oudenaarden/jyeung/software"

cd $outdir

# Download (with wget or curl):
# wget https://resources.aertslab.org/cistarget/zsync_curl
curl -O https://resources.aertslab.org/cistarget/zsync_curl

# Make executable:
chmod a+x zsync_curl

# Display full path to zsync_curl.
ZSYNC_CURL="${PWD}/zsync_curl"
echo "${ZSYNC_CURL}"
