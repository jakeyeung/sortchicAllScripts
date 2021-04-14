######################################################################
### How to generate a refseq TSS collection:


# python2.7 refseq2txt.py $inf $outf
ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/measure_cut_distances.TSS/refseq2txt.py"
# inf="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_genomic.gff"
inf="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_genomic.gbff"
outbase="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_genomic.parsed"
outf=${outbase}.parsed.init.txt
outftss=${outbase}.parsed.TSS.txt
outftes=${outbase}.parsed.TES.txt
outfcds=${outbase}.parsed.CDS.txt
outf="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_genomic.parsed.txt"
echo "Running python..."
# hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf > $outf
python $ps $inf > $outf
echo "Running python... done"

# 1. Tss
awk 'BEGIN{FS=OFS="\t"} $5 == "mRNA" && $7 ~ "NM_" {if($4 == "+"){tss = $2}else{tss = $3}; print $1, "TSS", tss, $4, 1, $7 ".." $6}' $outf | sort -k1,1 -k3,3n -k4,4 > ${outftss}

# 2. Transcription end sites
awk 'BEGIN{FS=OFS="\t"} $5 == "mRNA" && $7 ~ "NM_" {if($4 == "+"){tss = $3}else{tss = $2}; print $1, "TES", tss, $4, 1, $7 ".." $6}' $outf | sort -k1,1 -k3,3n -k4,4  > ${outftes}

# 1. CDS
awk 'BEGIN{FS=OFS="\t"} $5 == "CDS" && $7 ~ "NP_" {if($4 == "+"){tss = $2}else{tss = $3}; print $1, "CDS", tss, $4, 1, $7 ".." $6}' $outf | sort -k1,1 -k3,3n -k4,4 > ${outfcds}
