#!/bin/sh
# Jake Yeung
# download_erythryoblast_bigwigs.sh
# Download erythryblast bigwigs 
# 2019-03-20

# inurl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM946549"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_GenomeRes_2014"
# inurl="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM946549&format=file&file=GSM946549%5Fmm9%5FwgEncodePsuHistoneErythroblH3k09me3BE14halfCd1InputPk%2EbroadPeak%2Egz"

url1="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946549/suppl/GSM946549_mm9_wgEncodePsuHistoneErythroblH3k09me3BE14halfCd1InputRepSignalRep1.bigWig"
url2="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946549/suppl/GSM946549_mm9_wgEncodePsuHistoneErythroblH3k09me3BE14halfCd1InputRepSignalRep2.bigWig"
url3="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946549/suppl/GSM946549_mm9_wgEncodePsuHistoneErythroblH3k09me3BE14halfCd1InputSig.bigWig"

url4="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946547/suppl/GSM946547_mm9_wgEncodePsuHistoneErythroblH3k27me3BE14halfCd1InputRepSignalRep1.bigWig"
url5="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946547/suppl/GSM946547_mm9_wgEncodePsuHistoneErythroblH3k27me3BE14halfCd1InputRepSignalRep2.bigWig"
url6="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946547/suppl/GSM946547_mm9_wgEncodePsuHistoneErythroblH3k27me3BE14halfCd1InputSig.bigWig"

url7="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946524/suppl/GSM946524_mm9_wgEncodePsuHistoneErythroblH3k04me3BE14halfCd1InputRepSignalRep1.bigWig"
url8="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946524/suppl/GSM946524_mm9_wgEncodePsuHistoneErythroblH3k04me3BE14halfCd1InputRepSignalRep2.bigWig"
url9="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946524/suppl/GSM946524_mm9_wgEncodePsuHistoneErythroblH3k04me3BE14halfCd1InputSig.bigWig"

url10="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946536/suppl/GSM946536_mm9_wgEncodePsuHistoneErythroblH3k04me1BE14halfCd1InputRepSignalRep1.bigWig"
url11="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946536/suppl/GSM946536_mm9_wgEncodePsuHistoneErythroblH3k04me1BE14halfCd1InputRepSignalRep2.bigWig"
url12="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946536/suppl/GSM946536_mm9_wgEncodePsuHistoneErythroblH3k04me1BE14halfCd1InputSig.bigWig"

url13="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946523/suppl/GSM946523_mm9_wgEncodePsuHistoneMegakaryoH3k27me3BE14halfCd1InputRepSignalRep1.bigWig"
url14="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946523/suppl/GSM946523_mm9_wgEncodePsuHistoneMegakaryoH3k27me3BE14halfCd1InputRepSignalRep2.bigWig"
url15="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946523/suppl/GSM946523_mm9_wgEncodePsuHistoneMegakaryoH3k27me3BE14halfCd1InputSig.bigWig"

url16="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946525/suppl/GSM946525_mm9_wgEncodePsuHistoneMegakaryoH3k04me1BE14halfCd1InputRepSignalRep1.bigWig"
url17="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946525/suppl/GSM946525_mm9_wgEncodePsuHistoneMegakaryoH3k04me1BE14halfCd1InputRepSignalRep2.bigWig"
url18="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946525/suppl/GSM946525_mm9_wgEncodePsuHistoneMegakaryoH3k04me1BE14halfCd1InputSig.bigWig"

url19="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946527/suppl/GSM946527_mm9_wgEncodePsuHistoneMegakaryoH3k04me3BE14halfCd1InputRepSignalRep1.bigWig"
url20="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946527/suppl/GSM946527_mm9_wgEncodePsuHistoneMegakaryoH3k04me3BE14halfCd1InputRepSignalRep2.bigWig"
url21="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM946nnn/GSM946527/suppl/GSM946527_mm9_wgEncodePsuHistoneMegakaryoH3k04me3BE14halfCd1InputSig.bigWig"

cd $outdir

for i in $(echo {1..21}); do
    # echo $i
    inurl_var=url${i}
    inurl=$(eval echo "\$$inurl_var")
    # echo $inurl
    curl -O $inurl
done
