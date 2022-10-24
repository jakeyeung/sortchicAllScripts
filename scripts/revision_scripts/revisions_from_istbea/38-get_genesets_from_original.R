# Jake Yeung
# Date of Creation: 2022-05-21
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/38-get_genesets_from_original.R
# 


# Load gene sets ----------------------------------------------------------


inf.genes.k4me1 <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/geneset_H3K4me1_TSS_topics2.filt.colorcoded2.txt"
dat.genes.k4me1 <- fread(inf.genes.k4me1)

inf.genes <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/celltype_specific_genes_defined_by_K4me3_TSS.txt"
dat.genes <- fread(inf.genes)


# Process genes -----------------------------------------------------------



dat.genes.all <- bind_rows(dat.genes, subset(dat.genes.k4me1, select = c(-topic)))

dat.genes.sub <- subset(dat.genes, !(jset == "Basophils" & rnk > 100))
dat.genes.sub <- subset(dat.genes.sub, !(jset == "pDCs" & rnk > 100))

dat.genes.sub <- subset(dat.genes, rnk < 150)

# swap basophils
dat.genes.sub <- subset(dat.genes.sub, jset != "Basophils" & jset != "pDCs")

basos.k4me3 <- subset(dat.genes, jset %in% c("Basophils"))$gene
basos.k4me1 <- subset(dat.genes.k4me1, jset %in% c("Basophils"))$gene

pdcs.k4me3 <- subset(dat.genes, jset %in% c("pDCs"))$gene
pdcs.k4me1 <- subset(dat.genes.k4me1, jset %in% c("pDCs"))$gene

basos.manual <- c("Il4", "Il6", "Cpa3", "Il1r1")
basos.rname <- subset(dat.genes.all, symbol %in% basos.manual)$gene
basos.both <- c(intersect(basos.k4me1, basos.k4me3), basos.rname)
pdcs.both <- intersect(pdcs.k4me1, pdcs.k4me3)

dat.to.add <- subset(dat.genes.all, gene %in% c(basos.both, pdcs.both))

dat.genes.sub.join <- bind_rows(dat.genes.sub, dat.to.add)

dat.genes.sub.join$jset <- factor(dat.genes.sub.join$jset, levels = c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs"))

dat.genes.sub.join <- dat.genes.sub.join %>%
  arrange(jset, rnk)



# Saev output -------------------------------------------------------------
outf <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/dat_genes_sub_join_from_orig.rds"
saveRDS(dat.genes.sub.join, file = outf)
