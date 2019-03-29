#!/bin/sh
# Jake Yeung
# bedgraph_liftover_bigwig.sh
# wig_liftover_bigwig.sh
# Neutrophil data is bedGraph mm9. Convert to mm10 then bigwig
# Convert to mm10 then make bigwig files 
# https://www.biostars.org/p/81185/
# liftover, bedgraph to bigwig

inf=$1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
# chain=$2
# [[ ! -e $chain ]] && echo "$chain not found, exiting" && exit 1
outbase=$2  # should have an extension like .bw
[[ ! -d $outbase ]] && mkdir $outbase

chain="/Users/yeung/data/databases/liftover_chains/mm9ToMm10.over.chain.gz"
mm9size="/Users/yeung/data/databases/chrom_sizes/mm9.chrom.sizes"
mm10size="/Users/yeung/data/databases/chrom_sizes/mm10.chrom.sizes"

[[ ! -e $mm9size ]] && echo "$mm9size not found, exiting" && exit 1
[[ ! -e $mm10size ]] && echo "$mm10size not found, exiting" && exit 1
[[ ! -e $chain ]] && echo "$chain not found, exiting" && exit 1

# put things in tmp file
bname=$(basename $inf)
bname=${bname%%.*}

tmpdir="/tmp"

bedoutmm9=$tmpdir/$bname.mm9.bedgraph
bedoutmm10=$tmpdir/$bname.mm10.bedgraph
bedoutmm10mergedsorted=$tmpdir/$bname.mm10.mergedsorted.bedgraph
bwoutmm10=$outbase/$bname.mm10.bw

[[ -e $bwoutmm10 ]] && echo "$bwoutmm10  found, exit quietly" && exit 0

# first line is trackhub=, which should be removed or commented out https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/_kDGyPgWOmM
gzcat $inf | awk '{if (NR!=1) {print}}' > $bedoutmm9

echo "liftOver $bedoutmm9 $chain $bedoutmm10 $outbase/$bname.unmapped"
liftOver $bedoutmm9 $chain $bedoutmm10 $outbase/$bname.unmapped
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1

# errors with overlapping regions?
# https://groups.google.com/forum/#!topic/bedtools-discuss/6CJ1436E2uY
# d -1 to not merge book-end reads, which are almost all the lines
sort -k1,1 -k2,2n $bedoutmm10 | bedtools merge -d -1 -c 4 -o mean > $bedoutmm10mergedsorted

echo "bedGraphToBigWig $bedoutmm10mergedsorted $mm10size $outbase/$bname.bw"
bedGraphToBigWig $bedoutmm10mergedsorted $mm10size $outbase/$bname.bw
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1

