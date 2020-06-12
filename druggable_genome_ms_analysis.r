library(readr)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggrepel)

# read in data, making sure cols are in right format
setwd("/data/Wolfson-UKBB-Dobson/SMR/results")
druggable_genome = read_csv("druggable_genome.csv")
eqtlgen = read_table2("overall_ms_smr_results.txt",col_types=cols(
  probeID = col_character(),
  ProbeChr = col_double(),
  Gene = col_character(),
  Probe_bp = col_double(),
  topSNP = col_character(),
  topSNP_chr = col_double(),
  topSNP_bp = col_double(),
  A1 = col_character(),
  A2 = col_character(),
  Freq = col_double(),
  b_GWAS = col_double(),
  se_GWAS = col_double(),
  p_GWAS = col_double(),
  b_eQTL = col_double(),
  se_eQTL = col_double(),
  p_eQTL = col_double(),
  b_SMR = col_double(),
  se_SMR = col_double(),
  p_SMR = col_double(),
  p_HEIDI = col_double(),
  nsnp_HEIDI = col_double()))
cage = read_table2("overall_ms_cage_smr_results.txt",col_types=cols(
  probeID = col_character(),
  ProbeChr = col_double(),
  Gene = col_character(),
  Probe_bp = col_double(),
  topSNP = col_character(),
  topSNP_chr = col_double(),
  topSNP_bp = col_double(),
  A1 = col_character(),
  A2 = col_character(),
  Freq = col_double(),
  b_GWAS = col_double(),
  se_GWAS = col_double(),
  p_GWAS = col_double(),
  b_eQTL = col_double(),
  se_eQTL = col_double(),
  p_eQTL = col_double(),
  b_SMR = col_double(),
  se_SMR = col_double(),
  p_SMR = col_double(),
  p_HEIDI = col_double(),
  nsnp_HEIDI = col_double()))

meth_ms = read_table2("overall_meth_ms_smr_results.txt",col_types=cols(
  probeID = col_character(),
  ProbeChr = col_double(),
  Gene = col_character(),
  Probe_bp = col_double(),
  topSNP = col_character(),
  topSNP_chr = col_double(),
  topSNP_bp = col_double(),
  A1 = col_character(),
  A2 = col_character(),
  Freq = col_double(),
  b_GWAS = col_double(),
  se_GWAS = col_double(),
  p_GWAS = col_double(),
  b_eQTL = col_double(),
  se_eQTL = col_double(),
  p_eQTL = col_double(),
  b_SMR = col_double(),
  se_SMR = col_double(),
  p_SMR = col_double(),
  p_HEIDI = col_double(),
  nsnp_HEIDI = col_double()))


meth_expr = read_table2("overall_meth_expr_smr_results.txt",col_types=cols(
  Expo_ID = col_character(),
  Expo_Chr = col_double(),
  Expo_Gene = col_character(),
  Expo_bp = col_double(),
  Outco_ID = col_character(),
  Outco_Chr = col_double(),
  Outco_Gene = col_character(),
  Outco_bp = col_double(),
  topSNP = col_character(),
  topSNP_chr = col_double(),
  topSNP_bp = col_double(),
  A1 = col_character(),
  A2 = col_character(),
  Freq = col_double(),
  b_Outco = col_double(),
  se_Outco = col_double(),
  p_Outco = col_double(),
  b_Expo = col_double(),
  se_Expo = col_double(),
  p_Expo = col_double(),
  b_SMR = col_double(),
  se_SMR = col_double(),
  p_SMR = col_double(),
  p_HEIDI = col_double(),
  nsnp_HEIDI = col_double()
))


geuvadis = read_table2("overall_ms_geuvadis_smr_results.txt",col_types=cols(
  probeID = col_character(),
  ProbeChr = col_double(),
  Gene = col_character(),
  Probe_bp = col_double(),
  topSNP = col_character(),
  topSNP_chr = col_double(),
  topSNP_bp = col_double(),
  A1 = col_character(),
  A2 = col_character(),
  Freq = col_double(),
  b_GWAS = col_double(),
  se_GWAS = col_double(),
  p_GWAS = col_double(),
  b_eQTL = col_double(),
  se_eQTL = col_double(),
  p_eQTL = col_double(),
  b_SMR = col_double(),
  se_SMR = col_double(),
  p_SMR = col_double(),
  p_HEIDI = col_double(),
  nsnp_HEIDI = col_double()))


########################
# remove random bits of txt & get vars in right format
########################


eqtlgen = eqtlgen %>% mutate("ProbeChr" =  as.numeric(ProbeChr)) %>% filter(!is.na(ProbeChr))
cage = cage %>% mutate("ProbeChr" =  as.numeric(ProbeChr)) %>% filter(!is.na(ProbeChr))
meth_ms = meth_ms %>% mutate("ProbeChr" =  as.numeric(ProbeChr)) %>% filter(!is.na(ProbeChr))
meth_expr = meth_expr %>% mutate("Expo_Chr" =  as.numeric(Expo_Chr)) %>% filter(!is.na(Expo_Chr))

##########################
# eqtlgen
##########################


# merge with druggable genes
# nb this is merging on ensembl gene id
colnames(eqtlgen)[1]="ensembl_gene_id"
combo = eqtlgen %>% filter(ensembl_gene_id %in% druggable_genome$ensembl_gene_id) %>% 
  left_join(druggable_genome %>% select(hgnc_names,ensembl_gene_id,druggability_tier,start_b37,end_b37),by="ensembl_gene_id") %>% 
  filter(!is.na(b_SMR))

# get vars in numeric format

# exclude se mhc (hg 19 coords 6:25,000,000 - 35,000,000)
combo$topSNP_bp = as.numeric(as.character(combo$topSNP_bp))
combo$start_b37 = as.numeric(as.character(combo$start_b37))
combo$end_b37 = as.numeric(as.character(combo$end_b37))
combo$p_SMR = as.numeric(as.character(combo$p_SMR))

paste0("The number of tested genes before MHC filter is ",nrow(combo))

combo = combo %>% filter(!((ProbeChr == 6 & start_b37 >25000000 & start_b37 <35000000) |  (ProbeChr == 6 & end_b37 >25000000 & end_b37 <35000000)| (ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)))

# fdr adjusted p values for remaining probes
combo = combo %>% mutate(Q=p.adjust(p_SMR,method="fdr"))
paste0("The number of tested genes in total is ",nrow(combo))
combo = combo %>% filter(Q<0.05)
paste0("The number of genes passing Q<0.05 is ",nrow(combo))
combo = combo %>% filter(p_HEIDI>0.01)
paste0("The number of genes passing HEIDI > 0.05 is ",nrow(combo))
combo_tier1 = combo %>% filter(druggability_tier=="Tier 1")
paste0("The number of genes that are tier 1 druggable is ",nrow(combo_tier1))
write_csv(combo,"eqtl_sig_genes.csv")
table(combo$druggability_tier)
indiv_genes = combo %>% distinct(hgnc_names) 
paste0("The number of distinct genes is ",nrow(indiv_genes))

tbl = data.frame(table(combo$druggability_tier))
p=ggplot(tbl,aes(Var1,Freq,fill=Var1))+geom_col(colour="black")+
  theme_classic()+
  labs(x="Druggability Tier",y="Number of SMR-prioritised druggable genes")+
  scale_fill_manual(values = wes_palette("GrandBudapest1",4))+
  theme(legend.position = "none")

png("druggable_genes.png",height=8,width=8,units="in",res=300)
p
dev.off()

# plot direction of effect
plotting_genes = combo %>% select(hgnc_names, b_GWAS,b_eQTL,topSNP)

plotting_genes$b_eQTL = as.numeric(as.character(plotting_genes$b_eQTL))
plotting_genes$b_GWAS = as.numeric(as.character(plotting_genes$b_GWAS))
p=ggplot(plotting_genes,aes(b_eQTL,b_GWAS))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  geom_point()+geom_label_repel(aes(label=hgnc_names))+
  theme_classic()+
  labs(x="Effect on gene expression (Beta from eQTL data)",y="Effect on MS risk (Beta/log(OR) from MS GWAS)")

png("effect_of_genes_gwas_eqtl.png",height=8,width=8,units="in",res=300)
p
dev.off()

p2=ggplot(plotting_genes,aes(b_eQTL,b_GWAS))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  geom_point()+geom_label_repel(aes(label=topSNP))+
  theme_classic()+
  labs(x="Effect on gene expression (Beta from eQTL data)",y="Effect on MS risk (Beta/log(OR) from MS GWAS)")

png("effect_of_genes_gwas_eqtl_snplabels.png",height=8,width=8,units="in",res=300)
p2
dev.off()

##########################
# cage
##########################

# merge with druggable genes
# nb cage doesn't have ensembl ids, just illumina probe names
# this joins on hugo gene name
colnames(cage)[3]="hgnc_names"

combo_cage = cage %>% filter(hgnc_names %in% druggable_genome$hgnc_names) %>% 
  left_join(druggable_genome %>% select(hgnc_names,ensembl_gene_id,druggability_tier,start_b37,end_b37),by="hgnc_names") %>% 
  filter(!is.na(b_SMR))


# exclude mhc
combo_cage$topSNP_bp = as.numeric(as.character(combo_cage$topSNP_bp))
combo_cage$start_b37 = as.numeric(as.character(combo_cage$start_b37))
combo_cage$end_b37 = as.numeric(as.character(combo_cage$end_b37))
paste0("The number of tested genes before MHC filter is ",nrow(combo_cage))

combo_cage = combo_cage %>% filter(!((ProbeChr == 6 & start_b37 >25000000 & start_b37 <35000000) |  (ProbeChr == 6 & end_b37 >25000000 & end_b37 <35000000)| (ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)))

combo_cage = combo_cage %>% mutate(Q=p.adjust(p_SMR,method="fdr"))
paste0("The number of tested genes in total is ",nrow(combo_cage))
combo_cage = combo_cage %>% filter(Q<0.05)
paste0("The number of genes passing Q<0.05 is ",nrow(combo_cage))
combo_cage = combo_cage %>% filter(p_HEIDI>0.01)
paste0("The number of genes passing HEIDI > 0.05 is ",nrow(combo_cage))
combo_cage_tier1 = combo_cage %>% filter(druggability_tier=="Tier 1")
paste0("The number of genes that are tier 1 druggable is ",nrow(combo_cage_tier1))
write_csv(combo_cage,"cage_sig_genes.csv")
table(combo_cage$druggability_tier)


# cage plots
tbl = data.frame(table(combo_cage$druggability_tier))
p=ggplot(tbl,aes(Var1,Freq,fill=Var1))+geom_col(colour="black")+
  theme_classic()+
  labs(x="Druggability Tier",y="Number of SMR-prioritised druggable genes")+
  scale_fill_manual(values = wes_palette("GrandBudapest1",4))+
  theme(legend.position = "none")

png("cage_druggable_genes.png",height=8,width=8,units="in",res=300)
p
dev.off()

# plot direction of effect
plotting_genes = combo_cage %>% select(hgnc_names, b_GWAS,b_eQTL,topSNP)

plotting_genes$b_eQTL = as.numeric(as.character(plotting_genes$b_eQTL))
plotting_genes$b_GWAS = as.numeric(as.character(plotting_genes$b_GWAS))
p=ggplot(plotting_genes,aes(b_eQTL,b_GWAS))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  geom_point()+geom_label_repel(aes(label=hgnc_names))+
  theme_classic()+
  labs(x="Effect on gene expression (Beta from eQTL data)",y="Effect on MS risk (Beta/log(OR) from MS GWAS)")

png("cage_effect_of_genes_gwas_eqtl.png",height=8,width=8,units="in",res=300)
p
dev.off()

p2=ggplot(plotting_genes,aes(b_eQTL,b_GWAS))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  geom_point()+geom_label_repel(aes(label=topSNP))+
  theme_classic()+
  labs(x="Effect on gene expression (Beta from eQTL data)",y="Effect on MS risk (Beta/log(OR) from MS GWAS)")

png("cage_effect_of_genes_gwas_eqtl_snplabels.png",height=8,width=8,units="in",res=300)
p2
dev.off()
indiv_genes = combo_cage %>% distinct(hgnc_names) 
paste0("The number of distinct genes is ",nrow(indiv_genes))


# combine cage and eqtlgen
# this is to compare effect estimates from these two datasets

eqtl_summ = combo %>% select(ensembl_gene_id,druggability_tier,hgnc_names,b_SMR,p_SMR,p_HEIDI)
cage_summ = combo_cage %>% select(ensembl_gene_id,druggability_tier,hgnc_names,b_SMR,p_SMR,p_HEIDI) %>% arrange(p_SMR) %>% distinct(hgnc_names,.keep_all=TRUE)
cage_eqtl_summ = eqtl_summ %>% left_join(cage_summ,by=c("ensembl_gene_id","druggability_tier","hgnc_names")) %>% filter(!is.na(b_SMR.y)) %>% arrange(p_SMR.x)
write_csv(cage_eqtl_summ,"replicated_genes.csv")
indiv_genes = cage_eqtl_summ %>% distinct(hgnc_names) 
paste0("The number of distinct genes is ",nrow(indiv_genes))


# plot
cage_eqtl_summ$b_SMR.x = as.numeric(as.character(cage_eqtl_summ$b_SMR.x))
cage_eqtl_summ$b_SMR.y = as.numeric(as.character(cage_eqtl_summ$b_SMR.y))
cage_eqtl_summ$p_SMR.x = as.numeric(as.character(cage_eqtl_summ$p_SMR.x))
p=ggplot(cage_eqtl_summ,aes(b_SMR.x,b_SMR.y))+geom_abline(slope=1,intercept=0)+geom_point()+
  geom_label_repel(aes(label=hgnc_names),max.iter=5000)+
  theme_classic()+labs(x="Beta(SMR) - eQTLgen",y="Beta(SMR) - CAGE")
png(file="cage_eqtlgen_replicated.png",width=8,height=8,units="in",res=300)
p
dev.off()

cor.test(cage_eqtl_summ$b_SMR.x,cage_eqtl_summ$b_SMR.y)

##########################
# full smr eqtlgen
##########################

# this is to run SMR for all genes (not just druggable)

eqtlgen$topSNP_bp = as.numeric((eqtlgen$topSNP_bp))
eqtlgen$p_SMR = as.numeric((eqtlgen$p_SMR))

full_smr = eqtlgen %>% filter(!(ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)) %>% mutate(Q=p.adjust(p_SMR,method="fdr")) %>% filter(Q<0.05) %>% filter(p_HEIDI>0.01)

# this bit combines the SMR results with a gene list from USCS browser
# all coords are hg19

glist = read_table2("gene_list_hg19")
glist = glist %>% select(hg19.ensGtp.gene,hg19.ensemblToGeneName.value) %>% distinct(hg19.ensGtp.gene,.keep_all=TRUE)
colnames(glist)=c("ensembl_gene_id","gene_name")
full_smr %>% left_join(glist,by="ensembl_gene_id")
full_smr = full_smr %>% left_join(glist,by="ensembl_gene_id") %>% arrange(Q)

write_csv(full_smr,"full_smr.csv")

##########################
# multiomic smr
##########################

# M --> E
colnames(meth_expr)[7]="ensembl_gene_id"

# filter M2E pairs to druggable genes only
# exclude mhc
# fdr adjust
# filter for Q < 0.05
# filter p HEIDI > 0.01

meth_expr=meth_expr %>% 
  filter(ensembl_gene_id %in% druggable_genome$ensembl_gene_id) %>%
  filter(!(Expo_Chr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)) %>% 
  mutate(Q=p.adjust(p_SMR,method="fdr")) %>% 
  filter(Q<0.05) %>% 
  filter(p_HEIDI>0.01) 
write_csv(meth_expr,"methylation_gene_pairs.csv")

# E --> D
# filter E2D associations to druggable genes only
# exclude mhc
# fdr adjust
# filter for Q < 0.05
# filter p HEIDI > 0.01

expr_ms = eqtlgen %>% filter(ensembl_gene_id %in% druggable_genome$ensembl_gene_id) %>%
  filter(!(ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)) %>% 
  mutate(Q=p.adjust(p_SMR,method="fdr")) %>% 
  filter(Q<0.05) %>% 
  filter(p_HEIDI>0.01)

# M --> D
# filter M2D associations to druggable genes only
# exclude mhc
# fdr adjust
# filter for Q < 0.05
# filter p HEIDI > 0.01

meth_ms = meth_ms %>% filter(!is.na(b_SMR)) 
print(paste0("There are ",nrow(meth_ms)," probes before MHC filtering"))
meth_ms = meth_ms %>% 
  filter(!(ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)) %>% 
  mutate(Q=p.adjust(p_SMR,method="fdr")) %>% 
  filter(Q<0.05) %>% 
  filter(p_HEIDI>0.01)
write_csv(meth_ms,"methylation_assoc_with_ms.csv")

# filter M2E pairs to only those for which there is an E2D association (i.e. M2E2D)
m_e_pairs = meth_expr %>% select(Expo_ID,ensembl_gene_id)
e_d_genes = expr_ms %>% select(ensembl_gene_id)
e_d_genes_with_m2e_pair = m_e_pairs %>% filter(ensembl_gene_id %in% e_d_genes$ensembl_gene_id)
m_d_probes = meth_ms %>% select(probeID)
# then filter this M2E2D list again to only those associations where there is a significant 
# M2D association
e_d_genes_with_m2e_pair_with_sig_m2d = e_d_genes_with_m2e_pair %>% filter(Expo_ID %in% m_d_probes$probeID)


meth_ms_druggable = meth_ms %>% filter(probeID %in% meth_expr$Expo_ID) %>% rename("Expo_ID" = probeID) %>% left_join(meth_expr,by="Expo_ID")
overall_omics = meth_ms_druggable %>% filter(ensembl_gene_id %in% expr_ms$ensembl_gene_id) %>% left_join(druggable_genome %>% select(hgnc_names,ensembl_gene_id,druggability_tier),by="ensembl_gene_id")
write_csv(overall_omics,"multiomic_prioritised_ms_genes.csv")




##########################
# multiomic smr with CAGE
##########################



expr_ms = cage %>% filter(hgnc_names %in% druggable_genome$hgnc_names) %>%
  filter(!(ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)) %>% 
  mutate(Q=p.adjust(p_SMR,method="fdr")) %>% 
  filter(Q<0.05) %>% 
  filter(p_HEIDI>0.01)

meth_ms = meth_ms %>% filter(!is.na(b_SMR)) %>% 
  filter(!(ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)) %>% 
  mutate(Q=p.adjust(p_SMR,method="fdr")) %>% 
  filter(Q<0.05) %>% 
  filter(p_HEIDI>0.01)
colnames(meth_expr)[7]="ensembl_gene_id"

meth_expr=meth_expr %>% 
  filter(ensembl_gene_id %in% druggable_genome$ensembl_gene_id) %>%
  filter(!(Expo_Chr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)) %>% 
  mutate(Q=p.adjust(p_SMR,method="fdr")) %>% 
  filter(Q<0.05) %>% 
  filter(p_HEIDI>0.01) 

meth_ms_druggable = meth_ms %>% filter(probeID %in% meth_expr$Expo_ID) %>% rename("Expo_ID" = probeID) %>% left_join(meth_expr,by="Expo_ID")
expr_ms = expr_ms %>% left_join(druggable_genome %>% select(hgnc_names,ensembl_gene_id,druggability_tier),by="hgnc_names")
overall_omics_cage = meth_ms_druggable %>% filter(ensembl_gene_id %in% expr_ms$ensembl_gene_id) %>% left_join(expr_ms,by="ensembl_gene_id")
write_csv(overall_omics_cage,"multiomic_prioritised_ms_genes_cage.csv")

replicated_multiomic_genes = overall_omics %>% distinct(ensembl_gene_id,.keep_all=TRUE) %>% filter(hgnc_names %in% overall_omics_cage$hgnc_names) %>% 
  filter(druggability_tier =="Tier 1") %>% 
  select(ensembl_gene_id,ProbeChr)
write_tsv(replicated_multiomic_genes,"replicated_multiomic_probelist",col_names=FALSE)




##########################
# geuvadis lcl
##########################
colnames(geuvadis)[1]="ensembl_gene_id"
colnames(geuvadis)[3]="hgnc_names"

ids = c()
for (i in 1:nrow(geuvadis)){
  ids <<- c(ids,str_split(geuvadis$ensembl_gene_id,pattern=fixed("."),2)[[i]][1])
}
geuvadis$ensembl_gene_id = ids

# merge with druggable genes

combo = geuvadis %>% filter(ensembl_gene_id %in% druggable_genome$ensembl_gene_id) %>% 
  left_join(druggable_genome %>% select(hgnc_names,ensembl_gene_id,druggability_tier,start_b37,end_b37),by="ensembl_gene_id") %>% 
  filter(!is.na(b_SMR))

# get vars in numeric format

# exclude mhc
combo$topSNP_bp = as.numeric(as.character(combo$topSNP_bp))
combo$start_b37 = as.numeric(as.character(combo$start_b37))
combo$end_b37 = as.numeric(as.character(combo$end_b37))
combo$p_SMR = as.numeric(as.character(combo$p_SMR))

paste0("The number of tested genes before MHC filter is ",nrow(combo))

combo = combo %>% filter(!((ProbeChr == 6 & start_b37 >25000000 & start_b37 <35000000) |  (ProbeChr == 6 & end_b37 >25000000 & end_b37 <35000000)| (ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)))

combo = combo %>% mutate(Q=p.adjust(p_SMR,method="fdr"))
paste0("The number of tested genes in total is ",nrow(combo))
combo = combo %>% filter(Q<0.05)
paste0("The number of genes passing Q<0.05 is ",nrow(combo))
combo = combo %>% filter(p_HEIDI>0.01)
paste0("The number of genes passing HEIDI > 0.05 is ",nrow(combo))
combo_tier1 = combo %>% filter(druggability_tier=="Tier 1")
paste0("The number of genes that are tier 1 druggable is ",nrow(combo_tier1))
write_csv(combo,"lcl_eqtl_sig_genes.csv")
table(combo$druggability_tier)
indiv_genes = combo %>% distinct(hgnc_names) 
paste0("The number of distinct genes is ",nrow(indiv_genes))

tbl = data.frame(table(combo$druggability_tier))
p=ggplot(tbl,aes(Var1,Freq,fill=Var1))+geom_col(colour="black")+
  theme_classic()+
  labs(x="Druggability Tier",y="Number of SMR-prioritised druggable genes")+
  scale_fill_manual(values = wes_palette("GrandBudapest1",4))+
  theme(legend.position = "none")

png("lcl_druggable_genes.png",height=8,width=8,units="in",res=300)
p
dev.off()

# plot direction of effect
plotting_genes = combo %>% select(hgnc_names, b_GWAS,b_eQTL,topSNP)

plotting_genes$b_eQTL = as.numeric(as.character(plotting_genes$b_eQTL))
plotting_genes$b_GWAS = as.numeric(as.character(plotting_genes$b_GWAS))
p=ggplot(plotting_genes,aes(b_eQTL,b_GWAS))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  geom_point()+geom_label_repel(aes(label=hgnc_names))+
  theme_classic()+
  labs(x="Effect on gene expression (Beta from eQTL data)",y="Effect on MS risk (Beta/log(OR) from MS GWAS)")

png("lcl_effect_of_genes_gwas_eqtl.png",height=8,width=8,units="in",res=300)
p
dev.off()

p2=ggplot(plotting_genes,aes(b_eQTL,b_GWAS))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  geom_point()+geom_label_repel(aes(label=topSNP))+
  theme_classic()+
  labs(x="Effect on gene expression (Beta from eQTL data)",y="Effect on MS risk (Beta/log(OR) from MS GWAS)")

png("lcl_effect_of_genes_gwas_eqtl_snplabels.png",height=8,width=8,units="in",res=300)
p2
dev.off()


expr_ms = eqtlgen %>% filter(ensembl_gene_id %in% druggable_genome$ensembl_gene_id) %>%
  filter(!(ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)) %>% 
  mutate(Q=p.adjust(p_SMR,method="fdr")) %>% 
  filter(Q<0.05) %>% 
  filter(p_HEIDI>0.01)

combo %>% filter(ensembl_gene_id %in% expr_ms$ensembl_gene_id) %>% left_join(expr_ms,by="ensembl_gene_id") %>% select(hgnc_names.x,b_SMR.x,b_SMR.y)



combo = geuvadis %>% left_join(eqtlgen,by="ensembl_gene_id") %>% filter(!is.na(b_SMR.x)) %>% filter(!is.na(b_SMR.y)) %>% filter(p_HEIDI.x>0.01) %>% filter(p_HEIDI.y>0.01) %>% select(hgnc_names,b_SMR.x,b_SMR.y,se_SMR.x,se_SMR.y,p_SMR.x,p_SMR.y)
cor.test(combo$b_SMR.x,combo$b_SMR.y)

plot_combo = combo  %>% mutate("diff"=b_SMR.x-b_SMR.y) %>% mutate("mean" = (b_SMR.x+b_SMR.y)/2) %>% select(diff,mean,hgnc_names)

ggplot(plot_combo,aes(mean,diff))+geom_point()
p=ggplot(combo,aes(b_SMR.x,b_SMR.y,col=-log10(p_SMR.y),size=-log10(p_SMR.x)))+geom_abline(slope=1,intercept=0)+
   geom_point()+theme_classic()+labs(x="Beta(SMR) - LCLs",y="Beta(SMR) - Blood")+geom_label_repel(aes(label=hgnc_names))
png(file="cage_eqtlgen_replicated.png",width=8,height=8,units="in",res=300)
p
dev.off()




############################
# b cells
############################
# merge with druggable genes
colnames(bcell)[1]="ensembl_gene_id"
combo = bcell %>% filter(ensembl_gene_id %in% druggable_genome$ensembl_gene_id) %>% 
  left_join(druggable_genome %>% select(hgnc_names,ensembl_gene_id,druggability_tier,start_b37,end_b37),by="ensembl_gene_id") %>% 
  filter(!is.na(b_SMR))

# get vars in numeric format

# exclude mhc
combo$topSNP_bp = as.numeric(as.character(combo$topSNP_bp))
combo$start_b37 = as.numeric(as.character(combo$start_b37))
combo$end_b37 = as.numeric(as.character(combo$end_b37))
combo$p_SMR = as.numeric(as.character(combo$p_SMR))

paste0("The number of tested genes before MHC filter is ",nrow(combo))

combo = combo %>% filter(!((ProbeChr == 6 & start_b37 >25000000 & start_b37 <35000000) |  (ProbeChr == 6 & end_b37 >25000000 & end_b37 <35000000)| (ProbeChr == 6 & topSNP_bp >25000000 & topSNP_bp <35000000)))

combo = combo %>% mutate(Q=p.adjust(p_SMR,method="fdr"))
paste0("The number of tested genes in total is ",nrow(combo))
combo = combo %>% filter(Q<0.05)
paste0("The number of genes passing Q<0.05 is ",nrow(combo))
combo = combo %>% filter(p_HEIDI>0.01)
paste0("The number of genes passing HEIDI > 0.05 is ",nrow(combo))
combo_tier1 = combo %>% filter(druggability_tier=="Tier 1")
paste0("The number of genes that are tier 1 druggable is ",nrow(combo_tier1))
write_csv(combo,"eqtl_sig_genes.csv")
table(combo$druggability_tier)
indiv_genes = combo %>% distinct(hgnc_names) 
paste0("The number of distinct genes is ",nrow(indiv_genes))

