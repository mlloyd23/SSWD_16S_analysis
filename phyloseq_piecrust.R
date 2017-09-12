library(phyloseq)
library(DESeq2)
library(ggplot2)
theme_set(theme_bw())
library('biom')

x = read_biom("predicted_metagenome_L3.biom")
otumat = as(biom_data(x), "matrix")
OTU = otu_table(otumat, taxa_are_rows=TRUE)


mapping <- import_qiime_sample_data(mapfilename = 'map.txt')

phylo <- merge_phyloseq(OTU, mapping)
phylo

#########################################################
###Compare by phenotype
##########################################################
phylo_subset = subset_samples(phylo, Phenotype != "Dead")
sample_data(phylo_subset)$individual<-factor(sample_data(phylo_subset)$individual)
head(sample_data(phylo_subset)$Phenotype)

pheno = phyloseq_to_deseq2(phylo_subset, ~ individual + Phenotype)
pheno_results = DESeq(pheno, test="Wald")
pheno_res = results(pheno_results)
summary(pheno_res)
head(pheno_res)

alpha = 0.05
pheno_sigtab = pheno_res[which(pheno_res$padj < alpha), ]
head(pheno_sigtab)
write.table(pheno_sigtab, "pheno_L3.txt", sep="\t")



