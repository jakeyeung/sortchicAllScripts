#!/bin/sh
# Jake Yeung
# 2-liftover_coordinates.sh
# Liftover coorindates 
# 2016-07-22

downloaddir="/archive/epfl/upnae/jyeung/databases/chainfiles"
[[ ! -d $downloaddir ]] && mkdir $downloaddir
# download chainfiles
chainurl="http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz"
# liftOver copied from /home/yeung/bin/liftOver RSTUDIOSERVER

cd $downloaddir
mapchain=$downloaddir/"mm10Tomm9.over.chain.gz"
[[ ! -e $mapchain ]] && echo "Download chain files" && curl -o $mapchain $chainurl

inf="/archive/epfl/upnae/jyeung/MARA_promoters/mara_promoters_gene_name_association.cut.mm9.bed"  # mm9
outf="/archive/epfl/upnae/jyeung/MARA_promoters/mara_promoters_gene_name_association.mm10.bed"  # mm10
dir=$(dirname $outf)
bname=$(basename $outf)
bname=${bname%%.*}
unmapped=$dir/$bname.mm9tomm10.unmapped

liftOver $inf $mapchain $outf $unmapped
