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
goi<-c("830290","532994", "789777", "540617", "808031", "New.ReferenceOTU4597")

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

labels <- c("532994"= "532994\nPseudoalteromonas\npadj=1.44e-5", 
            "540617"= "540617\nPseudoalteromonas\npadj=0.0002",
	        	"789777"= "789777\nPseudoalteromonas\npadj=1.44e-5", 
	        	"808031"= "808031\nPseudoalteromonas\npadj=0.0015",
		        "830290"= "830290\nPseudoalteromonas\npadj=1.44e-5", 
	        	"New.ReferenceOTU4597"= "N.R.OTU4597\nPseudoalteromonas\npadj=0.0008")

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
neggoi<-c("748706", "New.ReferenceOTU3763", "81397", "575317", 
          "New.ReferenceOTU256","541771")

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

labels <- c("748706"= "748706\nTenacibaculum\npadj=1.78e-5",
            "New.ReferenceOTU3763"= "N.R.OTU6124\nPhaeobacter\npadj=0.0001",
		        "81397"= "81397\nPolaribacter\npadj=0.0002", 
		        "575317"= "575317\nPhaeobacter\npadj=0.0003",
	          "New.ReferenceOTU256"= "N.R.OTU256\nPhaeobacter\npadj=0.0003", 
		        "541771"= "541771\nPhaeobacter\npadj=0.0007")

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

