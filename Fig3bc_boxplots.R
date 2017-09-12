library(phyloseq)
library(DESeq2)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(tidyr)

otutable <- import_biom(BIOMfilename = 'biom_for_phyloseq.biom', 
                        treefilename = 'rep_set_no_chimeras.tre', 
                        parseFunction = parse_taxonomy_greengenes)
#warnings()

mapping <- import_qiime_sample_data(mapfilename = 'map.txt')

phylo <- merge_phyloseq(otutable, mapping)
phylo

####BENEFICIAL BOXPLOT
sample_data(phylo)$individual<-factor(sample_data(phylo)$individual)
sample_data(phylo)$pn<-factor(sample_data(phylo)$pn)

pheno_num <- phyloseq_to_deseq2(phylo, ~ individual + pn)
pheno_num_results <- DESeq(pheno_num, test="Wald")
pheno_num_res_0_1<- results(pheno_num_results, contrast=c("pn","0","1"),tidy=TRUE)
head(pheno_num_res_0_1)
poslog<-filter(pheno_num_res_0_1, log2FoldChange  > 0) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
print(tbl_df(poslog), n=23)

goi <- poslog$row[1:6]
goi


tcounts <- t(log2((counts(pheno_num_results[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(pheno_num_results), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

tcounts$Day<-as.factor(tcounts$Day)
tcounts$pn<-as.factor(tcounts$pn)

tcounts %>% 
  select(Row.names, Trajectory,Day, pn,Phenotype, gene, expression) %>% 
  head %>% 
  knitr::kable()

class(tcounts$gene)

palette<-c("#ffffcc", "#7fcdbb", "#2c7fb8")



labels <- c("532994"= c("532994\nPseudoalteromonas"), "540617"= "540617\nPseudoalteromonas",
		"789777"= "789777\nPseudoalteromonas", "808031"= "808031\nPseudoalteromonas",
		"830290"= "830290\nPseudoalteromonas", "New.ReferenceOTU4597"= "N.R.OTU4597\nPseudoalteromonas")

pdf("Beneficial_boxplot.pdf")
ggplot(tcounts, aes(pn, expression, fill=pn)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y",  labeller=labeller(gene = labels)) + 
  labs(x=NULL, 
       y="Expression (log normalized counts)", 
       fill="Symptom Stage")+
	theme(axis.title = element_text(color="#666666", face="bold", size=16)) +
	theme(legend.text=element_text(size=12)) + 
	theme(legend.title = element_text(colour="#666666", size=12, face="bold"))+
	scale_fill_manual(values= palette, labels=c("Healthy", "Early Stage", "Late Stage"))+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black"))+
	theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

###############Pathogenic
neglog<-filter(pheno_num_res_0_1, log2FoldChange  < 0) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
print(tbl_df(neglog), n=23)

neggoi <- neglog$row[1:6]
neggoi


tcounts <- t(log2((counts(pheno_num_results[neggoi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(pheno_num_results), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

tcounts$Day<-as.factor(tcounts$Day)
tcounts$pn<-as.factor(tcounts$pn)

tcounts %>% 
  select(Row.names, Trajectory,Day, pn,Phenotype, gene, expression) %>% 
  head %>% 
  knitr::kable()

class(tcounts$gene)

palette<-c("#ffffcc", "#7fcdbb", "#2c7fb8")



labels <- c("748706"= c("748706\nTenacibaculum"), "New.ReferenceOTU6124"= "N.R.OTU6124\nFlavobacteriaceae",
		"New.CleanUp.ReferenceOTU1371945"= "N.CU.R.OTU1371945\nTenacibaculum", "317388"= "317388\nFlavobacteriaceae",
		"79936"= "79936\nColwelliaceae", "New.ReferenceOTU4268"= "N.R.OTU4268\nNA")

pdf("Pathogenic_boxplot.pdf")
ggplot(tcounts, aes(pn, expression, fill=pn)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y",  labeller=labeller(gene = labels)) + 
  labs(x=NULL, 
       y="Expression (log normalized counts)", 
       fill="Symptom Stage")+
	theme(axis.title = element_text(color="#666666", face="bold", size=16)) +
	theme(legend.text=element_text(size=12)) + 
	theme(legend.title = element_text(colour="#666666", size=12, face="bold"))+
	scale_fill_manual(values= palette, labels=c("Healthy", "Early Stage", "Late Stage"))+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black"))+
	theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

