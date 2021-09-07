library("phyloseq")
library("ggplot2")
#install.packages("ggplot")
library("vegan")
library("DESeq2")
#library("pathview")
#install.packages("RCurl")
#install.packages("igraph")
#library(igraph)
library("RCy3")
library("dplyr")
#install.packages("dplyr")
#install.packages("stringgaussnet")
library(microbiome)
library(ggnet)
#install.packages("remotes")
#remotes::install_github("steinbaugh/DESeqAnalysis")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel")

ainstalled.packages()[, c("Package", "LibPath")]
BiocManager::install("RCy3")
aaif(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
Sys.setenv(R_INSTALL_STAGED = FALSE)
library(RCy3)
library(igraph)

source('../utility/cytoscape_util.R')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

install.packages("viridis")
library(viridis)

devtools::install_github("briatte/ggnet")
library(ggnet)


source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
browseVignettes("DESeq2")

#setwd("/home/hl46161/new_investigation_living_mulch/")

#import otu table and convert to matrix. This is the otu table after mitochondira and chloroplast filtration, no blank
otu_table	=	read.csv("/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.txt",sep="\t",row.names =1,check.names = FALSE)

#import otu table and convert to matrix 
otu_table	=	as.matrix(otu_table)

#taxonomy come from qiime2 artifact and needs to seperate in domain, phylum,order,......
taxonomy = read.csv("/home/hl46161/publish_living_mulch/taxonomy.tsv",sep="\t",row.names =1)
#taxonomy = read.csv("taxonomy.tsv",sep="\t",row.names =1)
taxonomy	=	as.matrix(taxonomy)
#import phylogenetic tree table and convert to matrix
phy_tree	=	read_tree("/home/hl46161/publish_living_mulch/tree.nwk")
#import metadata table and convert to matrix 
metadata	=	read.table("/home/hl46161/publish_living_mulch/living_mulch_data_metafile_massalin.tsv",sep ="\t",header = 1,row.names=1)
#metadata <- metadata[-c(4),]

#convert the data into phyloseq format 
OTU	=	otu_table(otu_table,taxa_are_rows	=	TRUE)
TAX	=	tax_table(taxonomy)
META = sample_data(metadata)
#make sure taxonomy and otu table match. The differences come from filtered mitochondira and cloroplast 
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

sample_names(OTU)
sample_names(META)

##build phyloseq object 
ps	=	phyloseq(OTU,TAX,phy_tree,META)
## create rarefaction plot to see the read depth for each samples 
rarecurve(t(otu_table(ps)), step=50, cex=0.5)
## save the rarefraction graph if necessary 
ggsave("rarefraction")


## create rarefied phyloseq object  size = minimium sample size * 0.9 
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)
## seperate data set into different time points 
ps.mid = subset_samples(ps, Date_taken == "2018-06-28")
ps.early = subset_samples(ps, Date_taken == "2018-05-21")
ps.late = subset_samples(ps, Date_taken == "2018-08-31")

ps.lm = subset_samples(ps, treatment == "LivingMulch")

ps.late.rarefied = subset_samples(ps.rarefied, Date_taken == "2018-08-31")
sample_sums(ps)
sample_sums(ps.rarefied)
plot_bar(ps,fill = "phylum")

otu_table(ps)
#######################################
##alpha diversity
plot_richness(ps.rarefied, x="treatment", color="Date_taken", measures=c("Observed"))

plot_richness(ps.rarefied, x="Date_taken", color="treatment", measures=c("Observed","Shannon")) + geom_boxplot()

plot_richness(ps.rarefied, x="treatment", measures=c("Observed", "Shannon"),color = "Date_taken") + geom_boxplot()

#beta diversity test
rich = estimate_richness(ps.rarefied)
rich
pairwise.wilcox.test(rich$Observed, sample_data(ps.rarefied)$treatment)
pairwise.wilcox.test(rich$Shannon, sample_data(ps.rarefied)$treatment)
plot_richness(ps.rarefied, x="Date_taken",color="treatment", measures=c("Observed", "Shannon")) + geom_boxplot()

ps.late.rarefied <- subset_samples(ps.rarefied, Date_taken = "2018-08-31")
rich2 = estimate_richness(ps.late.rarefied)
pairwise.wilcox.test(rich2$Observed, sample_data(ps.late.rarefied)$treatment)
pairwise.wilcox.test(rich2$Shannon, sample_data(ps.late.rarefied)$treatment)

rich3 = estimate_richness(ps)
pairwise.wilcox.test(rich3$Observed, sample_data(ps)$treatment)
pairwise.wilcox.test(rich3$Shannon, sample_data(ps)$treatment)

##PCOA graph construction 
wunifrac_dist = phyloseq::distance(ps.late.rarefied, method="unifrac", weighted=F)
ordination = ordinate(ps.late.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied, ordination, color="treatment") + theme(aspect.ratio=1)
adonis(wunifrac_dist ~ sample_data(ps.rarefied)$treatment)


#######################################################################################

#differential abundance test 
#import phyloseq oject into deseq2
ds = phyloseq_to_deseq2(ps, ~ Treatment)
ds = DESeq(ds)
#set the detection threshold 
alpha = 0.05

##first calculate overall differential abundant OTU between No Cover and living mulch 
res = results(ds, contrast=c("Treatment","LivingMulch","NoCover"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
#merge the otu with taxonomy
lm_vs_nc_ps = cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))
#filter based on critieria
lm_vs_nc_filtered_1_ps <- subset(lm_vs_nc_ps, abs(log2FoldChange) > 0)
lm_vs_nc_filtered_1_ps
lm_vs_nc_filtered_1_ps <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
lm_vs_nc_filtered_1_ps

#output the results of lm vs nc 
write.csv(lm_vs_nc_filtered_1_ps,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/lm_VS_nocover_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and living mulch 
res2 = results(ds, contrast=c("Treatment","LivingMulch","CrimsonClover"), alpha=alpha)
res2 = res2[order(res2$padj, na.last=NA), ]
res_sig2 = res2[(res2$padj < alpha), ]
res_sig2


##first calculate overall differential abundant OTU between CerealRye and living mulch 
res4 = results(ds, contrast=c("Treatment","LivingMulch","CerealRye"), alpha=alpha)
res4 = res4[order(res4$padj, na.last=NA), ]
res_sig4 = res4[(res4$padj < alpha), ]
res_sig4
lm_vs_cr_4_ps = cbind(as(res_sig4, "data.frame"), as(tax_table(ps)[rownames(res_sig4), ], "matrix"))
lm_vs_cr_filtered_4_ps <- subset(lm_vs_cr_4_ps, baseMean > 0)
lm_vs_cr_filtered_4_ps
lm_vs_cr_filtered_4_ps <- subset(lm_vs_cr_filtered_4_ps,abs(log2FoldChange) > 0)
lm_vs_cr_filtered_4_ps

write.csv(lm_vs_cr_filtered_4_ps,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/lm_VS_cr_differential_OTU.csv")

##first calculate overall differential abundant OTU between CerealRye and Nocover 

res5 = results(ds, contrast=c("Treatment", "CerealRye","NoCover"), alpha=alpha)
res5 = res5[order(res5$padj, na.last=NA), ]
res_sig5 = res5[(res5$padj < alpha), ]
res_sig5
cr_vs_nc_5_ps = cbind(as(res_sig5, "data.frame"), as(tax_table(ps)[rownames(res_sig5), ], "matrix"))
cr_vs_nc__filtered_5_ps <- subset(cr_vs_nc_5_ps, baseMean > 0)
cr_vs_nc__filtered_5_ps
cr_vs_nc__filtered_5_ps <- subset(cr_vs_nc__filtered_5_ps,abs(log2FoldChange) > 0)
cr_vs_nc__filtered_5_ps

write.csv(cr_vs_nc__filtered_5_ps,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/cr_VS_nc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and Nocover 
res6 = results(ds, contrast=c("Treatment", "CrimsonClover","NoCover"), alpha=alpha)
res6 = res6[order(res6$padj, na.last=NA), ]
res_sig6 = res6[(res6$padj < alpha), ]
res_sig6
cc_vs_nc_6 = cbind(as(res_sig6, "data.frame"), as(tax_table(ps)[rownames(res_sig6), ], "matrix"))
cc_vs_nc__filtered_6_ps <- subset(cc_vs_nc_6, baseMean > 0)
cc_vs_nc__filtered_6_ps
cc_vs_nc__filtered_6_ps <- subset(cc_vs_nc__filtered_6_ps,abs(log2FoldChange) > 0)
cc_vs_nc__filtered_6_ps

write.csv(cc_vs_nc__filtered_6_ps,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/cc_VS_nc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and Nocover 
res7 = results(ds, contrast=c("Treatment", "CrimsonClover","CerealRye"), alpha=alpha)
res7 = res7[order(res7$padj, na.last=NA), ]
res_sig7 = res7[(res7$padj < alpha), ]
res_sig7
cc_vs_cr_7_ps = cbind(as(res_sig7, "data.frame"), as(tax_table(ps)[rownames(res_sig7), ], "matrix"))
cc_vs_cr__filtered_7_ps <- subset(cc_vs_cr_7_ps, baseMean > 0)
cc_vs_cr__filtered_7_ps
cc_vs_cr__filtered_7_ps <- subset(cc_vs_cr__filtered_7_ps,abs(log2FoldChange) > 0)
cc_vs_cr__filtered_7_ps 

write.csv(cc_vs_cr__filtered_7_ps,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/CC_VS_CR_differential_OTU.csv")


#generate a ven diagram
gg <- gplots::venn(list(LivingMulch=rownames(lm_vs_nc_ps),CerealRye = rownames(cc_vs_nc_6),CrimsonColver = rownames(cr_vs_nc_5_ps)))
ggsave(gg,"Venn_diagram", height=6, width=6, device="pdf")

####################################################################

ps.lm
ds2 = phyloseq_to_deseq2(ps, ~Date_taken)
ds2 = DESeq(ds2)
alpha = 0.05

res = results(ds2, contrast=c("Date_taken","2018-05-21","2018-08-31"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
lm_may_vs_Au = cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))
lm_may_vs_Au  <- subset(lm_may_vs_Au, abs(log2FoldChange) > 0)
lm_may_vs_Au 
lm_may_vs_Au  <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
lm_may_vs_Au 


res2 = results(ds2, contrast=c("Date_taken","2018-06-28","2018-08-31"), alpha=alpha)
res2 = res2[order(res2$padj, na.last=NA), ]
res_sig2 = res2[(res2$padj < alpha), ]
res_sig2
lm_june_vs_Au = cbind(as(res_sig2, "data.frame"), as(tax_table(ps)[rownames(res_sig2), ], "matrix"))
lm_june_vs_Au  <- subset(lm_june_vs_Au, abs(log2FoldChange) > 0)
lm_june_vs_Au 
lm_june_vs_Au  <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
lm_june_vs_Au 



#####################################################################


ds = phyloseq_to_deseq2(ps, ~ Date_taken)
ds = DESeq(ds)
otu_table(ps)
alpha = 0.05

res = results(ds2, contrast=c("Date_taken","2018-05-21","2018-08-31"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
may_vs_Au = cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))
may_vs_Au  <- subset(may_vs_Au, abs(log2FoldChange) > 0)
may_vs_Au 
may_vs_Au  <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
may_vs_Au 


res2 = results(ds2, contrast=c("Date_taken","2018-06-28","2018-08-31"), alpha=alpha)
res2 = res2[order(res2$padj, na.last=NA), ]
res_sig2 = res2[(res2$padj < alpha), ]
res_sig2
lm_june_vs_Au = cbind(as(res_sig2, "data.frame"), as(tax_table(ps)[rownames(res_sig2), ], "matrix"))
lm_june_vs_Au  <- subset(lm_june_vs_Au, abs(log2FoldChange) > 0)
lm_june_vs_Au 
lm_june_vs_Au  <- subset(lm_vs_nc_filtered_1_ps, baseMean > 0)
lm_june_vs_Au 


#####################################################calculate the differential abundant OTU based on taxnomic level 
?merge_taxa()

OrderFiltered <- subset_taxa(ps, Order != "NA")
ps.order = tax_glom(OrderFiltered , "Order")
otu_table(OrderFiltered)
otu_table(ps.order)
order_otu_table = as(otu_table(ps.order), "matrix")
order_otu_table  = as.data.frame(order_otu_table)
write.csv(order_otu_table,"/home/hl46161/publish_replicate_new_living_mulch/order_otu_table.csv")

FamilyFiltered <- subset_taxa(ps, Family != "NA")
ps.family = tax_glom(FamilyFiltered, "Family")
otu_table(FamilyFiltered)
otu_table(ps.family)
family_otu_table = as(otu_table(ps.family), "matrix")
family_otu_table  = as.data.frame(family_otu_table)
family_taxa_table <- as.data.frame(as( tax_table(ps.family), "matrix"))

#write.table(family_otu_table,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_otu_table.tsv",sep="\t")

#differential abundance test 

ds = phyloseq_to_deseq2(ps.family, ~ Treatment)
ds = DESeq(ds)
alpha = 0.05

##first calculate overall differential abundant OTU between No Cover and living mulch 
res = results(ds, contrast=c("Treatment","LivingMulch","NoCover"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
lm_vs_nc_ps.family = cbind(as(res_sig, "data.frame"), as(tax_table(ps.family)[rownames(res_sig), ], "matrix"))
lm_vs_nc_filtered_1_ps.family <- subset(lm_vs_nc_ps.family, abs(log2FoldChange) > 0)
lm_vs_nc_filtered_1_ps.family
lm_vs_nc_filtered_1_ps.family <- subset(lm_vs_nc_filtered_1_ps.family, baseMean > 0)


#write.csv(lm_vs_nc_filtered_1_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_VS_nc_differential_OTU.csv")


##calculate overall differential abundant OTU between CrimsonClover and living mulch 
res2 = results(ds, contrast=c("Treatment","LivingMulch","CrimsonClover"), alpha=alpha)
res2 = res2[order(res2$padj, na.last=NA), ]
res_sig2 = res2[(res2$padj < alpha), ]
res_sig2
lm_vs_cc_2_ps.family = cbind(as(res_sig2, "data.frame"), as(tax_table(ps.family)[rownames(res_sig2), ], "matrix"))
lm_vs_cc_2_ps.family <- subset(lm_vs_cc_2_ps.family, baseMean > 0)
lm_vs_cc_2_ps.family
lm_vs_cc_2_ps.family <- subset(lm_vs_cc_2_ps.family,abs(log2FoldChange) > 0)
lm_vs_cc_2_ps.family


#write.csv(lm_vs_cc_2_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_VS_cc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CerealRye and living mulch 
res4 = results(ds, contrast=c("Treatment","LivingMulch","CerealRye"), alpha=alpha)
res4 = res4[order(res4$padj, na.last=NA), ]
res_sig4 = res4[(res4$padj < alpha), ]
res_sig4
lm_vs_cr_4_ps.family = cbind(as(res_sig4, "data.frame"), as(tax_table(ps.family)[rownames(res_sig4), ], "matrix"))
lm_vs_cr_filtered_4_ps.family <- subset(lm_vs_cr_4_ps.family, baseMean > 0)
lm_vs_cr_filtered_4_ps.family
lm_vs_cr_filtered_4_ps.family <- subset(lm_vs_cr_filtered_4_ps.family,abs(log2FoldChange) > 0)
lm_vs_cr_filtered_4_ps.family

#write.csv(lm_vs_cr_filtered_4_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_VS_cr_differential_OTU.csv")


##first calculate overall differential abundant OTU between CerealRye and Nocover

res5 = results(ds, contrast=c("Treatment", "CerealRye","NoCover"), alpha=alpha)
res5 = res5[order(res5$padj, na.last=NA), ]
res_sig5 = res5[(res5$padj < alpha), ]
res_sig5
cr_vs_nc_5_ps.family = cbind(as(res_sig5, "data.frame"), as(tax_table(ps.family)[rownames(res_sig5), ], "matrix"))
cr_vs_nc__filtered_5_ps.family <- subset(cr_vs_nc_5_ps.family, baseMean > 0)
cr_vs_nc__filtered_5_ps.family
cr_vs_nc__filtered_5_ps.family <- subset(cr_vs_nc__filtered_5_ps.family,abs(log2FoldChange) > 0)
cr_vs_nc__filtered_5_ps.family

#write.csv(cr_vs_nc__filtered_5_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_cr_VS_nc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and Nocover 
res6 = results(ds, contrast=c("Treatment", "CrimsonClover","NoCover"), alpha=alpha)
res6 = res6[order(res6$padj, na.last=NA), ]
res_sig6 = res6[(res6$padj < alpha), ]
res_sig6
cc_vs_nc_6 = cbind(as(res_sig6, "data.frame"), as(tax_table(ps.family)[rownames(res_sig6), ], "matrix"))
cc_vs_nc__filtered_6_ps.family <- subset(cc_vs_nc_6, baseMean > 0)
cc_vs_nc__filtered_6_ps.family
cc_vs_nc__filtered_6_ps.family <- subset(cc_vs_nc__filtered_6_ps.family,abs(log2FoldChange) > 0)
cc_vs_nc__filtered_6_ps.family 

#write.csv(cc_vs_nc__filtered_6_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_cc_VS_nc_differential_OTU.csv")


##first calculate overall differential abundant OTU between CrimsonClover and Nocover 
res7 = results(ds, contrast=c("Treatment", "CrimsonClover","CerealRye"), alpha=alpha)
res7 = res7[order(res7$padj, na.last=NA), ]
res_sig7 = res7[(res7$padj < alpha), ]
res_sig7
cc_vs_cr_7_ps.family = cbind(as(res_sig7, "data.frame"), as(tax_table(ps.family)[rownames(res_sig7), ], "matrix"))
cc_vs_cr__filtered_7_ps.family <- subset(cc_vs_cr_7_ps.family, baseMean > 0)
cc_vs_cr__filtered_7_ps.family
cc_vs_cr__filtered_7_ps.family <- subset(cc_vs_cr__filtered_7_ps.family,abs(log2FoldChange) > 0)
cc_vs_cr__filtered_7_ps.family 

#write.csv(cc_vs_cr__filtered_7_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_cc_VS_cr_differential_OTU.csv")

#######################

#summarize differential abudant otu. First try to find out which dfOtu is shared between three treatment 
lm_vs_nc_filtered_1_ps.family$OTU_ID <- rownames(lm_vs_nc_filtered_1_ps.family)
cr_vs_nc__filtered_5_ps.family$OTU_ID <- rownames(cr_vs_nc__filtered_5_ps.family)
cc_vs_nc__filtered_6_ps.family$OTU_ID <- rownames(cc_vs_nc__filtered_6_ps.family)

lm_cr_cc_otu <- dplyr::inner_join(lm_vs_nc_filtered_1_ps.family,cc_vs_nc__filtered_6_ps.family,by = "OTU_ID")
lm_cr_cc_otu <- lm_cr_cc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,log2FoldChange.y,baseMean.x,padj.x)
lm_cr_cc_otu <- dplyr::inner_join(lm_cr_cc_otu, cr_vs_nc__filtered_5_ps.family,by = "OTU_ID")

#lm_cr_cc_otu <- lm_cr_cc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#lm_cr_cc_otu_100 <- subset(lm_cr_cc_otu,baseMean.x>=100)
#lm_cr_cc_otu_100 <- subset(lm_cr_cc_otu,log2FoldChange.x>0)

#lm_cr_cc_otu <- lm_cr_cc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,log2FoldChange.y,log2FoldChange,baseMean.x,padj.x)
colnames(lm_cr_cc_otu) <-c("Kingdom","Phylum","Class","Order","Family","OTU_ID","log2FoldChange.lm","log2FoldChange.cc","log2FoldChange.cr","baseMean","P_adjusted")


lm_cr_cc_otu$log2FoldChange.lm <- round(lm_cr_cc_otu$log2FoldChange.lm, digits = 2)
lm_cr_cc_otu$log2FoldChange.cc <- round(lm_cr_cc_otu$log2FoldChange.cc, digits = 2)
lm_cr_cc_otu$log2FoldChange.cr <- round(lm_cr_cc_otu$log2FoldChange.cr, digits = 2)
#lm_cr_cc_otu$P_adjusted <- round(lm_cr_cc_otu$P_adjusted, digits = 2)

#lm_cr_cc_otu %>%
#mutate(Taxa = if_else( grepl("[[:digit:]]", lm_cr_cc_otu$Order),Order,Family))
#grepl("[[:digit:]]", lm_cr_cc_otu$Family)
#gsub("[[:digit:]]",lm_cr_cc_otu$Order,lm_cr_cc_otu$Family, fixed = TRUE)
?gsub()

lm_cr_cc_otu$Taxa <- as.character(lm_cr_cc_otu$Family)
lm_cr_cc_otu$Taxa[c(2,3,18,19,22,25,26,27,28,29,30,31,32,34,36,37,38,43)] <- c("c__Ktedonobacteria","o__Ktedonobacterales","p__Chloroflexi","o__Tepidisphaerales",
                                                                                " o__Rokubacteriales"," c__Planctomycetes"," c__Ktedonobacteria"," p__WPS-2","c__Phycisphaerae	",
                                                                                "c__Ktedonobacteria","p__Acidobacteriota"," c__Gammaproteobacteria","o__Acidobacteriales",
                                                                                "o__Gaiellales","c__Vicinamibacteria","p__Acidobacteriota","o__Elsterales",
                                                                               "c__Anaerolineae")



lm_cr_cc_otu_table <- lm_cr_cc_otu %>%
  mutate(Abudance =
           case_when(baseMean <= 100 ~ "*",
                     baseMean > 1000 ~ "***",
                     baseMean >100 ~ "**",
                     ))

lm_cr_cc_otu_table  <- lm_cr_cc_otu_table  %>% select(log2FoldChange.lm:log2FoldChange.cr,Taxa,Abudance)
 


#output the shared OTU between three cover crop treatment
write.csv(lm_cr_cc_otu,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_cc_cr_VS_nc_differential_OTU.csv",row.names = FALSE)

write.csv(lm_cr_cc_otu_table,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_table_lm_cr_cc_vs_nc.csv",row.names = FALSE)
  
  
lm_cr_cc_otu$CerealRye <- 1
lm_cr_cc_otu$LivingMulch <- 1
lm_cr_cc_otu$CrimsonClover <- 1

lm_cc_otu <- dplyr::inner_join(lm_vs_nc_filtered_1_ps.family, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
lm_cc_otu <- lm_cc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cc_otu <- dplyr::left_join(lm_cc_otu,cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
lm_cc_otu <- lm_cc_otu[rowSums(is.na(lm_cc_otu)) > 4,]
lm_cc_otu <- lm_cc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cc_otu$CerealRye <- 0
lm_cc_otu$LivingMulch <- 1
lm_cc_otu$CrimsonClover <- 1



lm_cr_otu <- dplyr::inner_join(lm_vs_nc_filtered_1_ps.family, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
lm_cr_otu <- lm_cr_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cr_otu <- dplyr::left_join(lm_cr_otu, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
lm_cr_otu <- lm_cr_otu[rowSums(is.na(lm_cr_otu)) > 4,]
lm_cr_otu <- lm_cr_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cr_otu$CerealRye <- 1
lm_cr_otu$LivingMulch <- 1
lm_cr_otu$CrimsonClover <- 0



cc_cr_otu <- dplyr::inner_join(cc_vs_nc__filtered_6_ps.family, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
cc_cr_otu <- cc_cr_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cc_cr_otu <- dplyr::left_join(cc_cr_otu,lm_vs_nc_filtered_1_ps.family, by = "OTU_ID")
cc_cr_otu <- cc_cr_otu[rowSums(is.na(cc_cr_otu)) > 4,]
#cc_cr_otu <- cc_cr_otu %>% select(Domain.x:Damily.x,OTU_ID)
cc_cr_otu <- cc_cr_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cc_cr_otu$CerealRye <- 1
cc_cr_otu$LivingMulch <- 0
cc_cr_otu$CrimsonClover <- 1




lm_nc_otu <- dplyr::left_join(lm_vs_nc_filtered_1_ps.family, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
lm_nc_otu <- lm_nc_otu[rowSums(is.na(lm_nc_otu)) > 7,]
lm_nc_otu <- lm_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_nc_otu <- dplyr::left_join(lm_nc_otu, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
lm_nc_otu <- lm_nc_otu[rowSums(is.na(lm_nc_otu)) > 7,]
lm_nc_otu <- lm_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#lm_nc_otu <- lm_nc_otu [,c(7:ncol(lm_nc_otu))]
lm_nc_otu$CerealRye <- 0
lm_nc_otu$LivingMulch <- 1
lm_nc_otu$CrimsonClover <- 0

cc_nc_otu <- dplyr::left_join(cc_vs_nc__filtered_6_ps.family, lm_vs_nc_filtered_1_ps.family, by = "OTU_ID")
cc_nc_otu <- cc_nc_otu[rowSums(is.na(cc_nc_otu)) > 7,]
cc_nc_otu <- cc_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cc_nc_otu <- dplyr::left_join(cc_nc_otu, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
cc_nc_otu <- cc_nc_otu[rowSums(is.na(cc_nc_otu)) > 7,]
cc_nc_otu <- cc_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#cc_nc_otu <- cc_nc_otu [,c(7:ncol(cc_nc_otu))]
cc_nc_otu$CerealRye <- 0
cc_nc_otu$LivingMulch <- 0
cc_nc_otu$CrimsonClover <- 1

cr_nc_otu <- dplyr::left_join(cr_vs_nc__filtered_5_ps.family, lm_vs_nc_filtered_1_ps.family, by = "OTU_ID")
cr_nc_otu <- cr_nc_otu[rowSums(is.na(cr_nc_otu)) > 7,]
cr_nc_otu <- cr_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cr_nc_otu <- dplyr::left_join(cr_nc_otu, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
cr_nc_otu <- cr_nc_otu[rowSums(is.na(cr_nc_otu)) > 7,]
cr_nc_otu <- cr_nc_otu %>% select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#cr_nc_otu <- cr_nc_otu [,c(7:ncol(cr_nc_otu))]
cr_nc_otu$CerealRye <- 1
cr_nc_otu$LivingMulch <- 0
cr_nc_otu$CrimsonClover <- 0


otu_summary_matrix <- dplyr::bind_rows(lm_nc_otu,cc_nc_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,cr_nc_otu)

otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,lm_cr_cc_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,lm_cr_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,lm_cc_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,cc_cr_otu)
otu_summary_matrix$Confidence.x <- NULL

upset_data <- otu_summary_matrix %>% select(CerealRye:CrimsonClover)

upset(upset_data,
      point.size = 3.2, line.size = 2,
      mainbar.y.label = "Distrubution of differential family", sets.x.label = "Total # of differential family",
      text.scale =c(2, 2, 2, 2, 2, 2)
)

ggsave("differential_otu_upset_plot.pdf",width = 8,height = 8,device = "pdf")

grid.text(
  "@littlemissdata",
  x = 0.90,
  y = 0.05,
  gp = gpar(
    fontsize = 10,
    fontface = 10
  )
)

