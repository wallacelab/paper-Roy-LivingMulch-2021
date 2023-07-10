if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") 
#install.packages ("rlang")
#remove.packages("data.table")
#remove.packages("ggplot2")
#install.packages ("data.table")
#install.packages ("ggplot2")
install_github("vqv/ggbiplot")
install.packages ("svglite")
library(qiime2R)
library(ggplot2)
library(tidyverse)
library(devtools)
library(gridExtra)
library(ggpubr)
library(rstatix)
library(emmeans)
library(phyloseq)
library(svglite)
library(vegan)
vegan::version
version

#set the working directory 
setwd("/home/hl46161/publish_living_mulch/")

#read in metadata and shannon index 
metadata <- read.table("living_mulch_data_metafile_massalin.tsv", sep="\t",header = TRUE)
shannon<-read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") 
shannon

##read in evenness artifact
eveness <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000/evenness_vector.qza")
eveness <-eveness$data %>% rownames_to_column("SampleID") 

##read in observed otu artifact
observed_otu <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000/observed_features_vector.qza")
observed_otu <- observed_otu$data %>% rownames_to_column("SampleID")
#merge the data 
metadata <- merge(shannon,metadata, by = "SampleID")
metadata <- merge(eveness,metadata, by = "SampleID")
metadata <- merge(observed_otu,metadata, by = "SampleID")
gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))
shannon$data 

#################

# install.packages("rJava")
#install.packages("UpSetR")
#install.packages("tidyverse")
#install.packages("venneuler")
#install.packages("grid")

library(UpSetR)
library(tidyverse)
library(grid)


############plot the pielou evenness of each treatment seperate by sampling date 

#tiff("Pielou_Evenness.tiff", width = 8, height = 8, units = 'in', res = 300)

gg_pielou <-ggplot(metadata,aes(x=Treatment, y=pielou_evenness)) + 
  facet_grid(~`Date_taken`) +
  xlab("Treatment") + 
  facet_grid(~`Date_taken`) +
  ylab("Pielou_Evenness") + 
  geom_point(aes(color=Treatment),size=4) +
  theme(
    legend.text = element_text(color = "black", size = 20),   #set the text size to be 20 and bold the text
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20, face="bold")
  ) + 
  theme(
    panel.grid.major = element_blank(),   ##remove grid line and background 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 16)) +
  geom_line() 

gg_pielou

#dev.off()

#save the plot in /home/hl46161/publish_living_mulch/ directroy 
ggsave("Pielou_Evenness.pdf", height=8, width=8, device="pdf",dpi = 600)
ggsave("Pielou_Evenness.png", height=8, width=8, device="png",dpi = 600)
ggsave("Pielou_Evenness.svg", height=8, width=8, device="svg",dpi = 600)


############plot the shannon diversity of each treatment seperate by sampling date 

#tiff("shannon_entropy.tiff", width = 8, height = 8, units = 'in', res = 300)

gg_shannon <-ggplot(metadata,aes(x=Treatment, y=shannon_entropy)) + 
  facet_grid(~`Date_taken`) +
  xlab("Treatment") + 
  facet_grid(~`Date_taken`) +
  ylab("Shannon_Entropy") + 
  geom_point(aes(color=Treatment),size=4) +
  theme(
    legend.text = element_text(color = "black", size = 20),       #set the text size to be 20 and bold the text
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20, face="bold")
  ) + 
  theme(
    panel.grid.major = element_blank(),   ##remove grid line and background
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 20)) +
  geom_line() 

gg_shannon

#dev.off()

#save the plot in /home/hl46161/publish_living_mulch/ directroy 
ggsave("shannon_entropy.pdf", height=8, width=8, device="pdf",dpi = 600)
ggsave("shannon_entropy.png", height=8, width=8, device="png",dpi = 600)
ggsave("shannon_entropy.svg", height=8, width=8, device="svg",dpi = 600)

###############################################
#merge shannon plot and pielou plot into one 
alpha_diversity = ggarrange(gg_shannon,gg_pielou,
                            legend = "top",
                            nrow = 1,
                            common.legend = TRUE)

alpha_diversity

ggsave("/home/hl46161/publish_living_mulch/alpha_diversity.png", height=8, width=16, device="png",dpi = 600)
ggsave("/home/hl46161/publish_living_mulch/alpha_diversity.pdf", height=8, width=16, device="pdf",dpi = 600)
ggsave("/home/hl46161/publish_living_mulch/alpha_diversity.svg", height=8, width=16, device="svg",dpi = 600)
ggsave("/home/hl46161/publish_living_mulch/alpha_diversity.tiff", height=8, width=16, device="tiff",dpi = 400)

?ggarrange()
################################ Try to use linear regression to determine the effect of the treatment and sampling date on alpha diversity 

treatment_order <- metadata$Treatment
#change the reference level of treatment
treatment_order <- relevel(metadata$Treatment, ref="NoCover")

#run linear regression of shannon index on treatment and date taken and their interaction
lm_treatment_date_relation_shannon <- lm(formula = shannon_entropy ~ treatment_order + Date_taken, data = metadata)

#summary the model 
summary(lm_treatment_date_relation_shannon)
### assumption check 
gvlma::gvlma(lm_treatment_date_relation_shannon)
plot(lm_treatment_date_relation_shannon)


#######################################################################
###assumption of model above is not satisfied, drop 11th and 12th samples according to plots of previous function
### 11th and 12th sample are GA18-21S-LM1400 and GA18-22S-600
treatment_order
treatment_order <- treatment_order[-c(11,12)]

#### rerun the model after filteration 
lm_treatment_date_relation_shannon <- lm(formula = shannon_entropy ~ treatment_order*Date_taken, data = metadata[-c(11,12),])

######summary the model
summary(lm_treatment_date_relation_shannon)
#anova(lm_treatment_date_relation_shannon)
Anova(lm_treatment_date_relation_shannon,type=3)
### assumption check 
gvlma::gvlma(lm_treatment_date_relation_shannon)
plot(lm_treatment_date_relation_shannon)

##since the iteraction is not significant, we can test treatment and Date_taken sperately using type 2 anova from car R Package
lm_treatment_date_relation_shannon <- lm(formula = shannon_entropy ~ treatment_order  + Date_taken, data = metadata[-c(11,12),])

######summary the model 
summary(lm_treatment_date_relation_shannon)
Anova(lm_treatment_date_relation_shannon,type=2)

### assumption check 
gvlma::gvlma(lm_treatment_date_relation_shannon)
#plot(lm_treatment_date_relation_shannon)

emm1 = emmeans(lm_treatment_date_relation_shannon, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts


################################### run the model of evenness on treatment and date_taken, and their interaction
lm_treatment_date_relation_evenness <- lm(formula = pielou_evenness ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_treatment_date_relation_evenness)
##########significant differences exist in iteraction of living mulch and June 28th and August 31 
### change the model to run type three Anova
lm_treatment_date_relation_evenness <- lm(formula = pielou_evenness ~ treatment_order*Date_taken, data = metadata[-c(11,12),],contrasts=list(treatment_order=contr.sum, Date_taken=contr.sum))
summary(lm_treatment_date_relation_evenness)
Anova(lm_treatment_date_relation_evenness,type=3)
gvlma::gvlma(lm_treatment_date_relation_evenness)
#########
lm_treatment_date_relation_evenness <- lm(formula = pielou_evenness ~ treatment_order + Date_taken, data = metadata[-c(11,12),])
summary(lm_treatment_date_relation_evenness)
Anova(lm_treatment_date_relation_evenness,type=2)
gvlma::gvlma(lm_treatment_date_relation_evenness)
#plot(lm_treatment_date_relation_evenness)



#######################################################
###no significant effect of either treatment and sampling date on # of ASV observed 
lm_treatment_date_relation_observed_otu <- lm(formula = observed_features ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_treatment_date_relation_observed_otu)

lm_treatment_date_relation_observed_otu <- lm(formula = observed_features ~ treatment_order + Date_taken, data = metadata[-c(11,12),])
summary(lm_treatment_date_relation_observed_otu)


################################################################
library(emmeans)
library(multcomp)
########## run a linear regression of soil physical characteristic on treatment and date_taken 
### The goal is to see whether treatment affect the soil physical characteritics 

treatment_order <- metadata[-c(11,12),]$Treatment
lm_pH <- lm(formula = pH ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_pH) # no significant iteraction, so test treatment and date seperately 

lm_pH <- lm(formula = pH ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_pH)

gvlma::gvlma(lm_pH)

Anova(lm_pH,type=3)

#Date, treatment, interaction
ph_vector <- c(1,0,0)

##### no significant differnce between treatment
emm1 = emmeans(lm_pH, specs = pairwise ~ treatment_order,adjust="FDR")
emm1

cld <- cld(emmeans(lm_pH, ~ treatment_order))
cld

cld <- cld(glht(lm_pH, linfct=mcp(treatment_order = "Tukey")))
cld

#find out the max value of each treatment of pH
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(pH))
value_max
### add the compact letter display to max value and organize into data frame 
value_max$letters <- c("A","A","A","A")
value_max <- as.data.frame(value_max)
###### create a treatment column that match with metadata so that the value max dataframe can match and be plotted on the graph
rownames(value_max) <- value_max$Treatment
value_max

#plot the pH value based on treatment and colored by date 
gg_pH <- ggplot(metadata,aes(x=Treatment,y=pH)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + 
  ylab("pH") + scale_x_discrete(
                                labels=c("CR", "CC", "LM","NC")) +
  theme(
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text = element_text(size=18, face="bold")
  ) + geom_text(data = value_max, aes(x=Treatment, y = 0.1+max_value,vjust=0),label=value_max$letters, size = 8) # input the compact letter display on the top of max value

gg_pH 

ggsave("pH.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("pH.png", height=8, width=10, device="png",dpi = 600)
ggsave("pH.svg", height=8, width=10, device="svg",dpi = 600)



#######################################################

lm_LBCEQ <- lm(formula = LBCEQ ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_LBCEQ) # significant iteraction exist between CrimsonClove and 2018-06-28

lm_LBCEQ <- lm(formula = LBCEQ ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_LBCEQ)

data = metadata

interaction.plot(x.factor = data$Date_taken, #x-axis variable
                 trace.factor = data$Treatment, #variable for lines
                 response = data$LBCEQ, #y-axis variable
                 fun = median, #metric to plot
                 ylab = "LBCEQ",
                 xlab = "Treatment",
                 col = c("pink", "blue","green","red"),
                 #lty = 1, #line type
                 lwd = 2, #line width
                 trace.label = "Date_taken")

######## All three cover crop treatment significant different from nocover 
emm1 = emmeans(lm_LBCEQ, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts

LBCEQ_vector <- c(1,1,1)

## generate compact letter display 
cld <- cld(glht(lm_LBCEQ, linfct=mcp(treatment_order = "Tukey")))
cld

cld <- cld(emmeans(lm_LBCEQ, ~ treatment_order))
cld

#find out the max value of each treatment of LBCEQ
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(LBCEQ))
value_max
### add the compact letter display to max value and organize into data frame 
value_max$letters <- c("A","A","A","B")
value_max <- as.data.frame(value_max)
###### create a treatment column that match with metadata so that the value max dataframe can match and be plotted on the graph
rownames(value_max) <- value_max$Treatment
value_max

#plot the LBCEQ value based on treatment and colored by date 
gg_LBCEQ <- ggplot(metadata,aes(x=Treatment,y=LBCEQ)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + 
  ylab("LBCEQ") + scale_x_discrete(
    labels=c("CR", "CC", "LM","NC")) +
  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 10+max_value,vjust=0),label=value_max$letters,size = 8) # input the compact letter display on the top of max value

gg_LBCEQ

ggsave("LBCEQ.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("LBCEQ.png", height=8, width=10, device="png",dpi = 600)
ggsave("LBCEQ.svg", height=8, width=10, device="svg",dpi = 600)

##################################################################################
### no significant difference of treatment in CEC  
lm_CEC <- lm(formula = CEC ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_CEC)
lm_CEC <- lm(formula = CEC ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_CEC)

CEC_vector <- c(0,0,0)
############################ 
emm1 = emmeans(lm_CEC, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts

cld <- cld(glht(lm_CEC, linfct=mcp(treatment_order = "Tukey")))
cld

#find out the max value of each treatment of LBCEQ
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(CEC))
value_max
### add the compact letter display to max value and organize into data frame 
value_max$letters <- c("A","AB","B","AB")
value_max <- as.data.frame(value_max)
###### create a treatment column that match with metadata so that the value max dataframe can match and be plotted on the graph
rownames(value_max) <- value_max$Treatment
value_max

#plot the LBCEQ value based on treatment and colored by date 
gg_CEC <- ggplot(metadata,aes(x=Treatment,y=CEC)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + 
  ylab("CEC") + scale_x_discrete(
    labels=c("CR", "CC", "LM","NC")) +
  theme(
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text = element_text(size=18, face="bold")
  ) + geom_text(data = value_max, aes(x=Treatment, y = 0.1+max_value,vjust=0),label=value_max$letters,size = 8) # input the compact letter display on the top of max value

gg_CEC

ggsave("CEC.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("CEC.png", height=8, width=10, device="png",dpi = 600)
ggsave("CEC.svg", height=8, width=10, device="tiff",dpi = 400)

#######

#################
####the model is not significant 
###no need to further explore  
lm_Ca <- lm(formula = Ca ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Ca)
lm_Ca <- lm(formula = Ca ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Ca)

Ca_vector <- c(0,1,0)

emm1 = emmeans(lm_Ca, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts
emm1$emmeans

cld(glht(lm_Ca, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of K
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(Ca))
value_max
value_max$letters <- c("A","A","A","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the K value based on treatment and colored by date 
gg_Ca <- ggplot(metadata,aes(x=Treatment,y=Ca)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + 
  ylab("Ca (ppm) ") + scale_x_discrete(
    labels=c("CR", "CC", "LM","NC")) +
  theme(
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text = element_text(size=18, face="bold")
  ) + geom_text(data = value_max, aes(x=Treatment, y = 1.05*max_value,vjust=0),label=value_max$letters,size = 8)

gg_Ca

ggsave("Ca.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("Ca.png", height=8, width=10, device="png",dpi = 600)
ggsave("Ca.svg", height=8, width=10, device="svg",dpi = 600)

####################
lm_K <- lm(formula = K ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_K)
####### no significant iteraction between treatment and date taken
lm_K <- lm(formula = K ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_K)

K_vector <- c(1,1,0)
######## All three cover crop treatment significant different from nocover 
emm1 = emmeans(lm_K, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts

cld(glht(lm_K, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of K
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(K))
value_max
value_max$letters <- c("AB","A","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the K value based on treatment and colored by date 
gg_K <- ggplot(metadata,aes(x=Treatment,y=K)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + 
  ylab("K (ppm) ") + scale_x_discrete(
    labels=c("CR", "CC", "LM","NC")) +
  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 1.05*max_value,vjust=0),label=value_max$letters, size = 8)

gg_K

ggsave("K.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("K.png", height=8, width=10, device="png",dpi = 600)
ggsave("K.svg", height=8, width=10, device="svg",dpi = 600)

################

lm_Mg <- lm(formula = Mg ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Mg)  ### significant iteraction exist between LivingMulch and 2018-06-28

lm_Mg <- lm(formula = Mg ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Mg) ### Living mulch and Crimsonclover significantly different from Nocover 

Mg_vector <- c(0,1,1)

data = metadata[-c(11,12),]

interaction.plot(x.factor = data$Date_taken, #x-axis variable
                 trace.factor = data$Treatment, #variable for lines
                 response = data$Mg, #y-axis variable
                 fun = median, #metric to plot
                 ylab = "Magesium",
                 xlab = "System",
                 col = c("pink", "blue","green","red"),
                 lty = 1, #line type
                 lwd = 2, #line width
                 trace.label = "Date_taken")


emm1 = emmeans(lm_Mg, specs = pairwise ~ treatment_order, adjust = "fdr")
emm1$contrasts

## generate compact letter display 
cld(glht(lm_Mg, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of K
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(Mg))
value_max$letters <- c("A","AB","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the Mg value based on treatment and colored by date 
gg_Mg <- ggplot(metadata,aes(x=Treatment,y=Mg)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + 
  ylab("Mg (ppm)") + scale_x_discrete(
    labels=c("CR", "CC", "LM","NC")) +
  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 5 + max_value,vjust=0),label=value_max$letters, size = 8)

gg_Mg
ggsave("Mg.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("Mg.png", height=8, width=10, device="png",dpi = 600)
ggsave("Mg.svg", height=8, width=10, device="svg",dpi = 600)

###################################################

lm_P <- lm(formula = P ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_P) 

lm_P <- lm(formula = P ~ treatment_order + Date_taken, data = metadata[-c(11,12),])
summary(lm_P) 

p_vector <- c(0,0,0)
  
emm1 = emmeans(lm_P, specs = pairwise ~ treatment_order,adjust = "fdr")
emm1$contrasts

## generate compact letter display 
cld(glht(lm_P, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of K
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(P))
value_max$letters <- c("A","A","A","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

gg_P <- ggplot(metadata,aes(x=Treatment,y=P)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + 
  ylab("P (ppm)") + scale_x_discrete(
    labels=c("CR", "CC", "LM","NC")) +
  theme(
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text = element_text(size=18, face="bold")
  ) + geom_text(data = value_max, aes(x=Treatment, y = 1.05*max_value,vjust=0),label=value_max$letters,size = 8)

gg_P
ggsave("P.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("P.png", height=8, width=10, device="png",dpi = 600)
ggsave("P.svg", height=8, width=10, device="svg",dpi = 600)
###################################

lm_Zn <- lm(formula = Zn ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Zn)  ### no significant iteraction exist 

lm_Zn <- lm(formula = Zn ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Zn)  ######## all three cover crop significantly different 

emm1 = emmeans(lm_Zn, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts

Zn_vector <- c(0,1,0)

## generate compact letter display 
cld(glht(lm_Zn, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of K
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(Zn))
value_max$letters <- c("B","B","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the Zn value based on treatment and colored by date 
gg_Zn <- ggplot(metadata,aes(x=Treatment,y=Zn)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + 
  ylab("Zn (ppm)") + scale_x_discrete(
    labels=c("CR", "CC", "LM","NC")) +
  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 0.2 + max_value,vjust=0),label=value_max$letters,size = 8)

gg_Zn
ggsave("Zn.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("Zn.png", height=8, width=10, device="png",dpi = 600)
ggsave("Zn.svg", height=8, width=10, device="svg",dpi = 600)

####################################

lm_Amon <- lm(formula = Amon ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Amon)  ### no significant iteraction exist 

lm_Amon <- lm(formula = Amon ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Amon)  #CereaRye amd Living mulch significant different 

emm1 = emmeans(lm_Amon, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts
## generate compact letter display 
cld(glht(lm_Amon, linfct=mcp(treatment_order = "Tukey")))

Amon_vector <- c(0,1,0)

#find out the max value of each treatment of Amon
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(Amon))
value_max$letters <- c("AB","AB","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the Amon value based on treatment and colored by date 
gg_Amon <- ggplot(metadata,aes(x=Treatment,y=Amon)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + ylab("Amon (ppm)") + scale_x_discrete(
  labels=c("CR", "CC", "LM","NC")) +
  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y = 1 + max_value,vjust=0),label=value_max$letters,size = 8)

gg_Amon
ggsave("Amon.pdf", height=8, width=10, device="pdf")
ggsave("Amon.png", height=8, width=10, device="png")
ggsave("Amon.tiff", height=8, width=10, device="tiff",dpi = 300)

###########################
lm_Nit <- lm(formula = Nit ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_Nit) ## no significant correlatiob between treatment and sampling date 
lm_Nit <- lm(formula = Nit ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_Nit) ### no tratment is significantly different, but sampling date are significantly different, which is understandable  

Nit_vector <- c(1,0,0)
  
emm1 = emmeans(lm_Nit, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts

cld(glht(lm_Nit, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of Nitrate
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(Nit))
value_max$letters <- c("A","A","A","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the Amon value based on treatment and colored by date 
gg_Nit <- ggplot(metadata,aes(x=Treatment,y=Nit)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + ylab("Nitrate (ppm)") + scale_x_discrete(
  labels=c("CR", "CC", "LM","NC")) +
  theme(
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text = element_text(size=18, face="bold")
  ) + geom_text(data = value_max, aes(x=Treatment, y = 1 + max_value,vjust=0,size = 8),label=value_max$letters,size = 8)

gg_Nit
ggsave("Nit.pdf", height=8, width=10, device="pdf")
ggsave("Nit.png", height=8, width=10, device="png")
ggsave("Nit.tiff", height=8, width=10, device="tiff",dpi = 300)


##############################

lm_N <- lm(formula = N ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_N)  #no significant iteraction exist 

lm_N <- lm(formula = N ~ treatment_order+Date_taken, data = metadata[-c(11,12),])
summary(lm_N) # all three treatment are significantly different from no cover 

N_vector <- c(0,1,0)

emm1 = emmeans(lm_N, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts
## generate compact letter display 
cld(glht(lm_N, linfct=mcp(treatment_order = "Tukey")))

#find out the max value of each treatment of N
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(N))
value_max$letters <- c("B","BC","C","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the N value based on treatment and colored by date 
gg_N <- ggplot(metadata,aes(x=Treatment,y=N)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + ylab("N (ppm)") + scale_x_discrete(
  labels=c("CR", "CC", "LM","NC")) +
  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y =  1.05*max_value,vjust=0),label=value_max$letters, size = 8)

gg_N
ggsave("N.pdf", height=8, width=10, device="pdf")
ggsave("N.png", height=8, width=10, device="png")
ggsave("N.tiff", height=8, width=10, device="tiff",dpi = 300)


#################################

lm_TOC <- lm(formula = TOC ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_TOC) # no significant iteraction

lm_TOC <- lm(formula = TOC ~ treatment_order + Date_taken, data = metadata[-c(11,12),])
summary(lm_TOC) #Living mulch is significant on TOC

TOC_vector <- c(0,1,0)
  
emm1 = emmeans(lm_TOC, specs = pairwise ~ treatment_order,adjust="FDR")
emm1$contrasts

## generate compact letter display 
cld(glht(lm_TOC, linfct=mcp(treatment_order = "Tukey")))
#find out the max value of each treatment of TOC
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(TOC))
value_max$letters <- c("A","AB","B","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the TOC value based on treatment and colored by date 
gg_TOC <- ggplot(metadata,aes(x=Treatment,y=TOC)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + ylab("TOC (ppm)") + scale_x_discrete(
  labels=c("CR", "CC", "LM","NC")) +
  theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text = element_text(size=18, face="bold")
) + geom_text(data = value_max, aes(x=Treatment, y =  1.05*max_value,vjust=0),label=value_max$letters, size = 8)

gg_TOC

ggsave("TOC.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("TOC.png", height=8, width=10, device="png",dpi = 600)
ggsave("TOC.svg", height=8, width=10, device="svg",dpi = 600)

#########################################################3

lm_BS <- lm(formula = BS ~ treatment_order*Date_taken, data = metadata[-c(11,12),])
summary(lm_BS) # no significant iteraction

lm_BS <- lm(formula = BS ~ treatment_order + Date_taken, data = metadata[-c(11,12),])
summary(lm_BS) #Living mulch is significant on BS

emm1 = emmeans(lm_BS, specs = pairwise ~ treatment_order, adjust = "FDR")
emm1$contrasts

BS_vector <- c(0,1,0)

?emmeans()
## generate compact letter display 
cld(glht(lm_BS, linfct=mcp(treatment_order = "Tukey")))
#find out the max value of each treatment of BS
value_max = metadata %>% group_by(Treatment) %>% summarize(max_value = max(BS))
value_max$letters <- c("A","A","A","A")
value_max <- as.data.frame(value_max)
rownames(value_max) <- value_max$Treatment
value_max

#plot the BS value based on treatment and colored by date 
gg_BS <- ggplot(metadata,aes(x=Treatment,y=BS)) +  geom_line() +  geom_point(aes(col=soil_test_date),size=6) + ylab("BS (%)") + scale_x_discrete(
  labels=c("CR", "CC", "LM","NC")) +
  theme(
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text = element_text(size=18, face="bold")
  ) + geom_text(data = value_max, aes(x=Treatment, y =  1.05*max_value,vjust=0),label=value_max$letters, size = 8)

gg_BS

ggsave("BS.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("BS.png", height=8, width=10, device="png",dpi = 600)
ggsave("BS.svg", height=8, width=10, device="svg",dpi = 600)

#################################################################

soil_physical_data <- ggarrange( gg_N,gg_LBCEQ,gg_TOC,gg_Amon,gg_K,gg_Mg,gg_Zn,gg_CEC,
                                 font.label = list(size = 20, color = "black", face = "bold", family = NULL),
                            ncol = 4, 
                            nrow = 2,
                            common.legend = TRUE)
soil_physical_data 

ggsave("soil_physical_data.pdf",height=16, width=20, device="pdf")
ggsave("soil_physical_data.png",height=16, width=20, device="png")
ggsave("soil_physical_data.svg",height=16, width=20, device="svg")
ggsave("soil_physical_data.tiff",height=16, width=20, device="tiff")

chemical_variable_significant <- rbind(ph_vector,p_vector,TOC_vector,Nit_vector,N_vector,Mg_vector,LBCEQ_vector,K_vector,CEC_vector,
                                       Ca_vector,BS_vector,Amon_vector)
chemical_variable_significant

colnames(chemical_variable_significant) <- c("Date_taken","Treatment","interaction")

chemical_variable_significant

###############################################################


Pvaluecal <- function(modelsummary){
  modelsummary <- unlist(modelsummary)
  f_value <- modelsummary$fstatistic.value
  numdf <- modelsummary$fstatistic.numdf
  dendf <- modelsummary$fstatistic.dendf
  p <- pf(f_value,numdf,dendf,lower.tail = F)
  return(p)
}

##################################################

p_value_adjust <- function(modelsummary,variable_name){
  
  p_value_list <- c()
  variable_name_list <- c()
  
  for(individual_result in modelsummary){
    
    P_value_result <- Pvaluecal(individual_result)
    P_value_result
    p_value_list <- append(p_value_list,P_value_result)
    
  }
  
  adjuested_p_value_list <- p.adjust(p_value_list, method = "BH", n = length(p_value_list))
  adjuested_p_value_list
  
  df <- data.frame(variable_name,adjuested_p_value_list,p_value_list)
  colnames(df) <- c("varaible_names","adjusted_p_value","raw_p_value")
  df
  
  return(df)
  
}

############################################################


p_value_adjust_cor <- function(modelsummary,variable_name){
  
  p_value_list <- c()
  
  for(individual_result in modelsummary){
    
    
    print(individual_result)
    p_value_list <- append(p_value_list,individual_result$p.value)
    
  }
  
  adjuested_p_value_list <- p.adjust(p_value_list, method = "BH", n = length(p_value_list))
  adjuested_p_value_list
  
  df <- data.frame(variable_name,adjuested_p_value_list,p_value_list)
  colnames(df) <- c("varaible_names","adjusted_p_value","raw_p_value")
  df
  
  return(df)
  
}

###########################################################


metadata_late = subset(metadata[-c(11,12),],Date_taken == "2018-08-31")

n <- 20

late_environment_group <- metadata_late[,c(8:20)]
#late_environment_group <- late_environment_group[,c("LBCEQ","pH","Ca","BS","CEC","K","Mg","Zn","Amon","TOC","N")]
late_environment_group
####### for loop the soil physical data column with evenness entropy 
lm_evenness_soil_august <-G2F_2019_environment_group
######### summarize the model
evenness_soil_august_summaries <- lapply(lm_evenness_soil_august, summary)
evenness_soil_august_summaries

########## The 2 3 5 7 13 column show significant correlation 
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"   
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
##############Therefore pH BS Ca Mg N significantly correlates with evenness 
########## check assumption 
assumption_check <- lapply(lm_evenness_soil_august, gvlma::gvlma)
assumption_check
########### 

lm_August_eveness_adjusted_p_value <- p_value_adjust(evenness_soil_august_summaries,colnames(late_environment_group))
lm_August_eveness_adjusted_p_value

cor_evenness_soil_august <- lapply(8:n, function(x) cor.test(metadata_late$pielou_evenness,metadata_late[,x],method="spearman"))
cor_evenness_soil_august

cor_August_eveness_adjusted_p_value <- p_value_adjust_cor(cor_evenness_soil_august,colnames(late_environment_group))
cor_August_eveness_adjusted_p_value
lm_August_eveness_adjusted_p_value


####### for loop the soil physical data column with shannon entropy 

lm_shannon_soil_august <- lapply(late_environment_group, function(x) lm(shannon_entropy ~ x, data = metadata_late))
######### summarize the model
shannon_soil_august_summaries <- lapply(lm_shannon_soil_august, summary)
shannon_soil_august_summaries
## pH is significant 
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"   
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
########## check assumption 
assumption_check <- lapply(lm_shannon_soil_august, gvlma::gvlma)
assumption_check
########## For august no soil physical data is associated with shannon entropy 

August_shannon_adjusted_p_value <- p_value_adjust(shannon_soil_august_summaries,colnames(late_environment_group))
August_shannon_adjusted_p_value



##############################################################

metadata_mid = subset(metadata[-c(11,12),],Date_taken == "2018-06-28")
n <- 20

mid_environment_group <- metadata_mid[,c(8:20)]
#mid_environment_group <- mid_environment_group[,c("LBCEQ","pH","Ca","BS","CEC","K","Mg","Zn","Amon","TOC","N")]

####### for loop the soil physical data column with shannon entropy 
lm_shannon_soil_mid <- lapply(mid_environment_group, function(x) lm(shannon_entropy ~ x, data = metadata_mid))
######### summarize the model
shannon_soil_mid_summaries <- lapply(lm_shannon_soil_mid, summary)
shannon_soil_mid_summaries
######### 10 13 models are significant
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"
### therefore nitrogen and ammonia is signifiant 
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
############## N and AMon is strongly associated with shannon entropy in June 
########## check assumption 
assumption_check <- lapply(lm_shannon_soil_mid, gvlma::gvlma)
assumption_check
########## 
mid_shannon_adjusted_p_value <- p_value_adjust(shannon_soil_mid_summaries,colnames(mid_environment_group))
mid_shannon_adjusted_p_value


####### for loop the soil physical data column with pielou evenness 
lm_evenness_soil_mid <- lapply(mid_environment_group, function(x) lm(pielou_evenness ~ x, data = metadata_mid))
######### summarize the model 
evenness_soil_mid_summaries <- lapply(lm_evenness_soil_mid, summary)
evenness_soil_mid_summaries
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
#"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"   
########## N and AMon is strongly associated with pielou evenness in June 
########## check assumption 
assumption_check <- lapply(lm_evenness_soil_mid, gvlma::gvlma)
assumption_check

mid_eveness_adjusted_p_value <- p_value_adjust(evenness_soil_mid_summaries,colnames(mid_environment_group))
mid_eveness_adjusted_p_value

cor_evenness_soil_mid <- lapply(8:n, function(x) cor.test(metadata_mid$pielou_evenness,metadata_mid[,x],method="spearman"))
cor_evenness_soil_mid

cor_mid_eveness_adjusted_p_value <- p_value_adjust_cor(cor_evenness_soil_mid,colnames(mid_environment_group))
cor_mid_eveness_adjusted_p_value

mid_eveness_adjusted_p_value

###############################################################

mid_environment_group <- metadata_mid[,c(8:20)]
mid_environment_group <- mid_environment_group[,c("LBCEQ","pH","CEC","K","Mg","Zn","Amon","TOC","N")]

####### for loop the soil physical data column with shannon entropy 
lm_shannon_soil_mid <- lapply(mid_environment_group, function(x) lm(shannon_entropy ~ x, data = metadata_mid))
######### summarize the model
shannon_soil_mid_summaries <- lapply(lm_shannon_soil_mid, summary)
shannon_soil_mid_summaries

########################################################################
metadata_early = subset(metadata[-c(11,12),],Date_taken == "2018-05-21")
n <- 20

early_environment_group <- metadata_early[,c(8:20)]


####### for loop the soil physical data column with shannon entropy 
lm_shannon_soil_early <- lapply(early_environment_group, function(x) lm(shannon_entropy ~ x, data = metadata_early))

######### summarize the model
shannon_soil_early_summaries <- lapply(lm_shannon_soil_early, summary)
shannon_soil_early_summaries
######### 9th models are significant
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N" 
#### therefore Zn is significant in shannon
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
########## check assumption 
assumption_check <- lapply(lm_shannon_soil_early, gvlma::gvlma)
assumption_check
########## 

lm_evenness_soil_early <- lapply(early_environment_group, function(x) lm(pielou_evenness ~ x, data = metadata_early))

######### summarize the model
evenness_soil_early_summaries <- lapply(lm_evenness_soil_early, summary)
evenness_soil_early_summaries
######### 9th models are significant
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N" 
#### therefore Zn is significant in shannon
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
########## check assumption 
assumption_check <- lapply(lm_evenness_soil_early, gvlma::gvlma)
assumption_check
########## 


###############################################################
#Test the association between soil physical data and alpha diversity, first start with evenness

metadata_late = subset(metadata[-c(11,12),],Date_taken == "2018-08-31")

n <- 20

late_environment_group <- metadata_late[,c(8:20)]
late_environment_group

####### for loop the soil physical data column with evenness entropy 
lm_evenness_soil_august <- lapply(8:n, function(x) lm(pielou_evenness ~ metadata_late[,x], data = metadata_late))
######### summarize the model
evenness_soil_august_summaries <- lapply(lm_evenness_soil_august, summary)
evenness_soil_august_summaries


########## The 2 3 5 7 13 column show significant correlation 
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"   
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
##############Therefore pH BS Ca Mg N significantly correlates with evenness 
########## check assumption 
assumption_check <- lapply(lm_evenness_soil_august, gvlma::gvlma)
assumption_check
########### 

####### for loop the soil physical data column with shannon entropy 
lm_shannon_soil_august <- lapply(8:n, function(x) lm(shannon_entropy ~ metadata_late[,x], data = metadata_late))
######### summarize the model
shannon_soil_august_summaries <- lapply(lm_shannon_soil_august, summary)
shannon_soil_august_summaries
## pH is significant 
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"   
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
########## check assumption 
assumption_check <- lapply(lm_shannon_soil_august, gvlma::gvlma)
assumption_check
########## For august no soil physical data is associated with shannon entropy 

#######################################################################

metadata_mid = subset(metadata[-c(11,12),],Date_taken == "2018-06-28")
n <- 20
####### for loop the soil physical data column with shannon entropy 
lm_shannon_soil_mid <- lapply(8:n, function(x) lm(shannon_entropy ~ metadata_mid[,x], data = metadata_mid))
######### summarize the model
shannon_soil_mid_summaries <- lapply(lm_shannon_soil_mid, summary)
shannon_soil_mid_summaries
######### 10 13 models are significant
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"
### therefore nitrogen and ammonia is signifiant 
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
############## N and AMon is strongly associated with shannon entropy in June 
########## check assumption 
assumption_check <- lapply(lm_shannon_soil_mid, gvlma::gvlma)
assumption_check
########## 

####### for loop the soil physical data column with pielou evenness 
lm_evenness_soil_mid <- lapply(8:n, function(x) lm(pielou_evenness ~ metadata_mid[,x], data = metadata_mid))
######### summarize the model 
evenness_soil_mid_summaries <- lapply(lm_evenness_soil_mid, summary)
evenness_soil_mid_summaries
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
#"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N"   
########## N and AMon is strongly associated with pielou evenness in June 
########## check assumption 
assumption_check <- lapply(lm_evenness_soil_mid, gvlma::gvlma)
assumption_check

###############################################

metadata_early = subset(metadata[-c(11,12),],Date_taken == "2018-05-21")
n <- 20
####### for loop the soil physical data column with shannon entropy 
lm_shannon_soil_early <- lapply(8:n, function(x) lm(shannon_entropy ~ metadata_early[,x], data = metadata_early))
lm_shannon_soil_early <- lapply(8:n, function(x) lm(pielou_evenness ~ metadata_early[,x], data = metadata_early))
######### summarize the model
shannon_soil_early_summaries <- lapply(lm_shannon_soil_early, summary)
shannon_soil_early_summaries
######### 9th models are significant
##"LBCEQ" "pH"    "BS"    "CEC"   "Ca"    "K"     "Mg"    "P"     "Zn"    "Amon"  "Nit"   "TOC"   "N" 
#### therefore Zn is significant in shannon
####### create soil column name vector so that it is easier to interpret the previous summarize list 
soil_column_name <- colnames(metadata[8:20])
soil_column_name
########## check assumption 
assumption_check <- lapply(lm_shannon_soil_early, gvlma::gvlma)
assumption_check
########## 


###############################################################
##plot mantel test matrix 

calc_distance_vector <- function(DM){
  
  rownames(DM) <- DM$X
  ##delete previous sample names column
  DM$X <-NULL
  ## reorder the dataframe by rownames and column names 
  DM <- DM[order(row.names(DM)),  order(colnames(DM))]
  
  #set up loop numbers 
  row_loop=nrow(DM)
  print(row_loop)
  column_loop=ncol(DM)-1
  print(column_loop)
  ##set up an emptry vector
  unifrac_distance_vector <- c(0)
  ###

  
  ###### loop trhough the lower triangle of the matrix and add elements to the vector 
  while(row_loop <=nrow(DM) & row_loop >1){
    for(i in column_loop:1){
      #print(i)
      #print(DM[row_loop,i])
      unifrac_distance_vector <- append(unifrac_distance_vector,DM[row_loop,i])
    }
    column_loop=column_loop-1
    row_loop = row_loop-1
  }
  ### remove first zero 
  unifrac_distance_vector <- unifrac_distance_vector[2:length(unifrac_distance_vector)]
  unifrac_distance_vector
  return(unifrac_distance_vector)
}
  
######################################################################3
## read in august weighted unifrac distance 
August_weighted_unifrac_distance <- read.csv("/home/hl46161/publish_living_mulch/core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-weighted_unifrac_distance-matrix.tsv",sep="\t",header=TRUE)
August_weighted_unifrac_distance_vector <- calc_distance_vector(August_weighted_unifrac_distance)
August_weighted_unifrac_distance_vector

August_unweighted_unifrac_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-unweighted_unifrac_distance-matrix.tsv",sep="\t",header=TRUE)
August_unweighted_unifrac_distance_vector <- calc_distance_vector(August_unweighted_unifrac_distance)
August_unweighted_unifrac_distance_vector

August_bray_curtis_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-bray-curtis_distance.tsv",sep="\t",header=TRUE)
August_bray_curtis_distance_vector <- calc_distance_vector(August_bray_curtis_distance)
August_bray_curtis_distance_vector

August_aitchison_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-aitchison-distance-matrix.tsv",sep="\t",header=TRUE)
August_aitchison_distance_vector <- calc_distance_vector(August_aitchison_distance)
August_aitchison_distance_vector

#####
## read in august N unifrac distance 
August_N_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-N-distance-matrix.tsv",sep="\t",header=TRUE)
August_N_distance_vector <- calc_distance_vector(August_N_distance)
August_N_distance_vector

##############################
## read in august N unifrac distance 
August_LBCEQ_distance <- read.csv("./core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August-LBCEQ-distance-matrix.tsv",sep="\t",header=TRUE)
August_LBCEQ_distance_vector <- calc_distance_vector(August_LBCEQ_distance)
August_LBCEQ_distance_vector

########################################

August_data <- data.frame(August_N_distance_vector,August_LBCEQ_distance_vector,August_weighted_unifrac_distance_vector,August_unweighted_unifrac_distance_vector,August_bray_curtis_distance_vector,August_bray_curtis_distance_vector)

summary(lm(August_weighted_unifrac_distance_vector ~ August_N_distance_vector,data=August_data))
summary(lm(August_unweighted_unifrac_distance_vector ~ August_N_distance_vector,data=August_data))
summary(lm(August_bray_curtis_distance_vector ~ August_N_distance_vector,August_data))
summary(lm(August_aitchison_distance_vector ~ August_N_distance_vector,August_data))

summary(lm(August_weighted_unifrac_distance_vector ~ August_LBCEQ_distance_vector,August_data))
summary(lm(August_unweighted_unifrac_distance_vector ~ August_LBCEQ_distance_vector,August_data))
summary(lm(August_bray_curtis_distance_vector ~ August_LBCEQ_distance_vector,August_data))
summary(lm(August_aitchison_distance_vector ~ August_LBCEQ_distance_vector,August_data))

#########################################
library(ggpmisc)
#plot two vector 
###
August_N <- ggplot(data=August_data, aes(x=August_N_distance_vector, y=August_weighted_unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  scale_color_manual(values=c("cyan")) +
  xlab("distance in soil N concentration") + ggtitle("August") +
  ylab("weighted unifrac distance") +
  theme(
    legend.text = element_text(color = "black", size = 30),
    legend.title = element_text(color = "black", size = 30),
    axis.title.x = element_text(size=30, face="bold"),
    axis.title.y = element_text(size=30, face="bold"),
    plot.title = element_text(size=30, face="bold"),
    axis.text.y = element_text(size=30, face="bold"),
    axis.text.x = element_text(size=30, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(0.042,0.04),
           y = c(0.14,0.12), 
           label = c("p-value = 6.34e-08","R ^ 2 = 0.3592"),size=10) 
August_N

ggsave("August_N_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("August_N_weighted_unifrac_distance.png", height=8, width=10, device="png")

August_LBCEQ <- ggplot(data=August_data, aes(x=August_LBCEQ_distance_vector, y=August_weighted_unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  scale_color_manual(values=c("cyan")) +
  xlab("distance in LBCEQ value") + 
  ylab("weighted unifrac distance") + ggtitle("August") +
  theme(
    legend.text = element_text(color = "black", size = 30),
    legend.title = element_text(color = "black", size = 30),
    axis.title.x = element_text(size=30, face="bold"),
    axis.title.y = element_text(size=30, face="bold"),
    plot.title = element_text(size=30, face="bold"),
    axis.text.y = element_text(size=30, face="bold"),
    axis.text.x = element_text(size=30, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(160,160),
           y = c(0.14,0.12), 
           label = c("p-value = 1.817e-05","R ^ 2 = 0.2395"),size=10) 

August_LBCEQ 

ggsave("August_LBCEQ_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("August_LBCEQ_weighted_unifrac_distance.png", height=8, width=10, device="png")

###########################################################

## read in august weighted unifrac distance 
may_weighted_unifrac_distance <- read.csv("./core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/may-weighted-unifrac-distance-matrix.tsv",sep="\t",header=TRUE)
may_weighted_unifrac_distance_vector <- calc_distance_vector(may_weighted_unifrac_distance)
may_weighted_unifrac_distance_vector

##
may_unweighted_unifrac_distance <- read.csv("./core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/may-unweighted-unifrac-distance-matrix.tsv",sep="\t",header=TRUE)
may_unweighted_unifrac_distance_vector <- calc_distance_vector(may_unweighted_unifrac_distance)
may_unweighted_unifrac_distance_vector

may_bray_curtis_distance <- read.csv("./core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/may-bray-curtis-distance-matrix.tsv",sep="\t",header=TRUE)
may_bray_curtis_distance_vector <- calc_distance_vector(may_bray_curtis_distance)
may_bray_curtis_distance_vector

may_aitchison_distance <- read.csv("/home/hl46161/publish_living_mulch/core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/May-aitchison-distance-matrix.tsv",sep="\t",header=TRUE)
may_aitchison_distance_vector <- calc_distance_vector(may_aitchison_distance)
may_aitchison_distance_vector

length(may_bray_curtis_distance_vector)

## read in august N unifrac distance 
may_Zn_distance <- read.csv("./core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/May-Zn-distance-matrix.tsv",sep="\t",header=TRUE)
##delete previous sample names column
rownames(may_Zn_distance) <- may_Zn_distance$X
may_Zn_distance$X <-NULL
may_Zn_distance <-may_Zn_distance[order(row.names(may_Zn_distance)),  order(colnames(may_Zn_distance))]
may_Zn_distance <- may_Zn_distance[-c(3),-c(3)]

ncol(may_Zn_distance)

#set up loop numbers 
row_loop=11
column_loop=10
##set up an emptry vector
may_Zn_distance_vector <- c(0)

######

###### loop trhough the lower triangle of the matrix and add elements to the vector 
while(row_loop <=11 & row_loop >1){
  for(i in column_loop:1){
    print(i)
    print(may_Zn_distance[row_loop,i])
    may_Zn_distance_vector <- append(may_Zn_distance_vector,may_Zn_distance[row_loop,i])
  }
  column_loop=column_loop-1
  row_loop = row_loop-1
}
### remove first zero 
may_Zn_distance_vector <- may_Zn_distance_vector[2:length(may_Zn_distance_vector)]
may_Zn_distance_vector
length(may_Zn_distance_vector)


############################

May_data <- data.frame(may_weighted_unifrac_distance_vector,may_unweighted_unifrac_distance_vector,may_bray_curtis_distance_vector,may_Zn_distance_vector,may_aitchison_distance_vector)
summary(lm(may_weighted_unifrac_distance_vector ~ may_Zn_distance_vector, data=May_data))
summary(lm(may_unweighted_unifrac_distance_vector ~ may_Zn_distance_vector, data=May_data))
summary(lm(may_bray_curtis_distance_vector ~ may_Zn_distance_vector, data=May_data))
summary(lm(may_aitchison_distance_vector ~ may_Zn_distance_vector, data=May_data))


###################################################################################################

ggplot(data=May_data, aes(x=may_Zn_distance_vector, y=may_weighted_unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  scale_color_manual(values=c("cyan")) +
  xlab("distance in Zn level") + 
  ylab("weighted unifrac distance") +
  theme(
    legend.text = element_text(color = "black", size = 30),
    legend.title = element_text(color = "black", size = 30),
    axis.title.x = element_text(size=30, face="bold"),
    axis.title.y = element_text(size=30, face="bold"),
    axis.text.y = element_text(size=30, face="bold"),
    axis.text.x = element_text(size=30, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(1.75,1.75),
           y = c(0.13,0.12), 
           label = c("p-value = 5.293e-06","R ^ 2 = 0.3134"),size=10) 

ggsave("Zn_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("Zn_weighted_unifrac_distance.png", height=8, width=10, device="png")


################

## read in august weighted unifrac distance 
June_weighted_unifrac_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June-weighted-distance-matrix.tsv",sep="\t",header=TRUE)
June_weighted_unifrac_distance_vector <- calc_distance_vector(June_weighted_unifrac_distance) 

June_unweighted_unifrac_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June-unweighted-distance-matrix.tsv",sep="\t",header=TRUE)
June_unweighted_unifrac_distance_vector <- calc_distance_vector(June_unweighted_unifrac_distance) 

June_bray_curtis_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June-bray-curtis-distance-matrix.tsv",sep="\t",header=TRUE)
June_bray_curtis_distance_vector <- calc_distance_vector(June_bray_curtis_distance) 

June_aitchison_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June-distance-matrix.tsv",sep="\t",header=TRUE)
June_aitchison_distance_vector <- calc_distance_vector(June_aitchison_distance)
June_aitchison_distance_vector


## read in august N unifrac distance 
June_N_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/N-June-distance-matrix.tsv",sep="\t",header=TRUE)
June_N_distance_vector <- calc_distance_vector(June_N_distance) 

## read in august N unifrac distance 
June_LBCEQ_distance <- read.csv("./core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June-LBCEQ-distance-matrix.tsv",sep="\t",header=TRUE)
June_LBCEQ_distance_vector <- calc_distance_vector(June_LBCEQ_distance) 

#########################################################

June_data <- data.frame(June_N_distance_vector,June_LBCEQ_distance_vector,June_weighted_unifrac_distance_vector,June_unweighted_unifrac_distance_vector,June_bray_curtis_distance_vector)

summary(lm(June_weighted_unifrac_distance_vector ~ June_N_distance_vector,data=June_data))
summary(lm(June_unweighted_unifrac_distance_vector ~ June_N_distance_vector,data=June_data))
summary(lm(June_bray_curtis_distance_vector ~ June_N_distance_vector,June_data))
summary(lm(June_aitchison_distance_vector ~ June_N_distance_vector,June_data))

summary(lm(June_weighted_unifrac_distance_vector ~ June_LBCEQ_distance_vector,June_data))
summary(lm(June_unweighted_unifrac_distance_vector ~ June_LBCEQ_distance_vector,June_data))
summary(lm(June_bray_curtis_distance_vector ~ June_LBCEQ_distance_vector,June_data))
summary(lm(June_aitchison_distance_vector ~ June_N_distance_vector,June_data))


####################################################

June_N <- ggplot(data=June_data, aes(x=June_N_distance_vector, y=June_weighted_unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  scale_color_manual(values=c("cyan")) +
  xlab("distance in soil N concentration") + 
  ylab("weighted unifrac distance") + ggtitle("June") +
  theme(
    legend.text = element_text(color = "black", size = 30),
    legend.title = element_text(color = "black", size = 30),
    axis.title.x = element_text(size=30, face="bold"),
    axis.title.y = element_text(size=30, face="bold"),
    plot.title = element_text(size=30, face="bold"),
    axis.text.y = element_text(size=30, face="bold"),
    axis.text.x = element_text(size=30, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(0.05,0.05),
           y = c(0.14,0.13), 
           label = c("p-value = 0.012","R ^ 2 = 0.1283"),size=10) 

June_N
ggsave("June_N_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("June_N_unifrac_distance.png", height=8, width=10, device="png")


June_LBCEQ <- ggplot(data=June_data, aes(x=June_LBCEQ_distance_vector, y=June_weighted_unifrac_distance_vector)) + geom_point(size=3,color="steelblue") + 
  xlab("distance in LBCEQ value") + 
  ylab("weighted unifrac distance") + ggtitle("June") +
  theme(
    legend.text = element_text(color = "black", size = 30),
    legend.title = element_text(color = "black", size = 30),
    axis.title.x = element_text(size=30, face="bold"),
    axis.title.y = element_text(size=30, face="bold"),
    plot.title = element_text(size=30, face="bold"),
    axis.text.y = element_text(size=30, face="bold"),
    axis.text.x = element_text(size=30, face="bold")
  ) + geom_smooth(method="lm",color="blue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  annotate("text",
           x = c(200,200),
           y = c(0.14,0.13), 
           label = c("p-value = 3.741e-05","R ^ 2 = 0.2229"),size=10) 

June_LBCEQ

ggsave("June_LBCEQ_weighted_unifrac_distance.pdf", height=8, width=10, device="pdf")
ggsave("June_LBCEQ_unifrac_distance.png", height=8, width=10, device="png")


#############################################
meta_soil_beta <- ggarrange(June_N,August_N,June_LBCEQ,August_LBCEQ,
          ncol = 2, 
          nrow = 2)
meta_soil_beta 

ggsave("soil_physical_data_correlation_weighted_unifrac_distance.pdf",meta_soil_beta,height=16, width=20, device="pdf",dpi = 600)
ggsave("soil_physical_data_correlation_weighted_unifrac_distance.png",meta_soil_beta,height=16, width=20, device="png",dpi = 600)
ggsave("soil_physical_data_correlation_weighted_unifrac_distance.svg",meta_soil_beta,height=16, width=20, device="svg",dpi = 600)
ggsave("soil_physical_data_correlation_weighted_unifrac_distance.tiff",meta_soil_beta,height=8, width=10, device="tiff",dpi = 300)
ggsave("soil_physical_data_correlation_weighted_unifrac_distance.tif",meta_soil_beta,height=16, width=20, device="tif",dpi = 300)


ggsave()
aes(color="blue")
?annotate()
?label()
?grid.arrange()

gg <- grid.arrange(June_N,August_N,June_LBCEQ,August_LBCEQ, name= c("June","August","June","August"))
gg
ggsave(gg,"soil_physical_data_correlation_weighted_unifrac_distance.pdf",height=16, width=20, device="pdf")
ggsave(gg,"soil_physical_data_correlation_weighted_unifrac_distance.png",height=16, width=20, device="png")
ggsave(gg,"soil_physical_data_correlation_weighted_unifrac_distance.svg",height=16, width=20, device="svg")

?ggsave()

##########################################
#make supplement graph for 
sample_frequency <- read.table("/home/hl46161/Downloads/sample-frequency-detail(2).csv",sep=",")
colnames(sample_frequency) <- c("SampleID","Frequency")
?geom_histogram()
ggplot(data=sample_frequency,aes(Frequency)) + geom_histogram(fill='blue',col="black", alpha=0.2) + 
  ggtitle("                    Histogram for sample reads size distribution") +
  scale_x_continuous(breaks=c(0,50000,100000,150000,200000,250000),limits = c(0,280000)) +
  scale_y_continuous(breaks=c(2,4,6,8,10,12,14),limits = c(0,14)) +
  xlab('Number of Reads') +
  ylab("Number of Sample") +
  theme(
    legend.text = element_text(color = "black", size = 18),
    legend.title = element_text(color = "black", size = 18),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold"),
    axis.text.y = element_text(size=18, face="bold"),
    axis.text.x = element_text(size=18, face="bold"),
    plot.title = element_text(size=18, face="bold")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("reads_distribution.pdf", height=8, width=10, device="pdf",dpi = 600)
ggsave("reads_distribution.png", height=8, width=10, device="png",dpi = 600)
ggsave("reads_distribution.svg", height=8, width=10, device="svg",dpi = 600)

? ggtitle()
########################################### make rarefraction graph
install.packages ("rjson")
library("reshape2")
library("rjson")
library(data.table)
library(tidyr)
library(stringr)
library(dplyr)

### read in the ASV count file 
observed_features_rarefration_data <- read.csv("/home/hl46161/publish_living_mulch/observed_features.csv")
### summarize the file based on sample ID
observed_features_rarefration_data_melt <- melt(observed_features_rarefration_data,id.vars = "sample.id")
### order the file based on sample ID 
observed_features_rarefration_data_melt <- observed_features_rarefration_data_melt[order(observed_features_rarefration_data_melt$sample.id),]
n <- 10
##### summarize mean value of AVS count of every 10 iteration at one sequence septh
observed_features_rarefration_data_melt <- observed_features_rarefration_data_melt %>% group_by(mean = (row_number() -1) %/% n) %>%
  mutate(mean = mean(value))

#### sustitute depth.1_iter.1 to contain only sequence depth 
observed_features_rarefration_data_melt$variable <- gsub("depth.","",observed_features_rarefration_data_melt$variable)
observed_features_rarefration_data_melt$variable <- gsub("_iter.10","",observed_features_rarefration_data_melt$variable)
observed_features_rarefration_data_melt$variable <- gsub("_iter.[0-9]","",observed_features_rarefration_data_melt$variable)
##### reselect the column 
observed_features_rarefration_data_melt <- observed_features_rarefration_data_melt[,c("sample.id","variable","mean")]
##### remove duplicate rows 
observed_features_rarefration_data_melt <- unique(observed_features_rarefration_data_melt)
##### rename the colulmn 
colnames(observed_features_rarefration_data_melt) <- c("SampleID","depth","ASV_count")
########
observed_features_rarefration_data_melt <- observed_features_rarefration_data_melt[order(observed_features_rarefration_data_melt$depth),]
###### set the depth to bu numeric. this step is important. you can try not to do it to see what happen to the graph
observed_features_rarefration_data_melt$depth <- as.numeric(observed_features_rarefration_data_melt$depth)
ggplot(data=observed_features_rarefration_data_melt,aes(x =depth, y = ASV_count,group=SampleID,colour=SampleID)) + geom_line() + geom_point() +
  ylim(0, 5000) +
  scale_x_continuous("Read Repth",breaks=c(0,5000,10000,15000,20000,25000,30000),limits = c(0,32000)) +
  ylab("ASV Count") +
  geom_text(x=25000, y=600, label="GA18-20S-600",color="brown") +
  theme(legend.position='none',axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text = element_text(size=15, face="bold")
  )

ggsave("rarefraction_curve_observed_ASV.pdf",height=10, width=8, device="pdf")
ggsave("rarefraction_curve_observed_ASV.png",height=10, width=8, device="png")
ggsave("rarefraction_curve_observed_ASV.svg",height=10, width=8, device="svg")

#############################################################################

### read in the shannon count file 
shannon_rarefration_data <- read.csv("/home/hl46161/publish_living_mulch/shannon.csv")
### summarize the file based on sample ID
shannon_rarefration_data_melt <- melt(shannon_rarefration_data,id.vars = "sample.id")
### order the file based on sample ID 
shannon_rarefration_data_melt <- shannon_rarefration_data_melt[order(shannon_rarefration_data_melt$sample.id),]
n <- 10
##### summarize mean value of AVS count of every 10 iteration at one sequence septh
shannon_rarefration_data_melt <- shannon_rarefration_data_melt %>% group_by(mean = (row_number() -1) %/% n) %>%
  mutate(mean = mean(value))

#### sustitute depth.1_iter.1 to contain only sequence depth 
shannon_rarefration_data_melt$variable <- gsub("depth.","",shannon_rarefration_data_melt$variable)
shannon_rarefration_data_melt$variable <- gsub("_iter.10","",shannon_rarefration_data_melt$variable)
shannon_rarefration_data_melt$variable <- gsub("_iter.[0-9]","",shannon_rarefration_data_melt$variable)
shannon_rarefration_data_melt <- shannon_rarefration_data_melt[,c("sample.id","variable","mean")]
shannon_rarefration_data_melt <- unique(shannon_rarefration_data_melt)
##### rename the colulmn 
colnames(shannon_rarefration_data_melt) <- c("SampleID","depth","Shannon_index")
shannon_rarefration_data_melt <- shannon_rarefration_data_melt[order(shannon_rarefration_data_melt$depth),]
###### set the depth to bu numeric. this step is important. you can try not to do it to see what happen to the graph
shannon_rarefration_data_melt$depth <- as.numeric(shannon_rarefration_data_melt$depth)
ggplot(data=shannon_rarefration_data_melt,aes(x =depth, y = Shannon_index,group=SampleID,color=SampleID)) + geom_line() + geom_point() + 
  ylim(0, 12) + 
  ylab("Shannon Index") +
  geom_text(x=25000, y=7.8, label="GA18-20S-600",color="brown") +
  scale_x_continuous("Read Depth",breaks=c(0,5000,10000,15000,20000,25000,30000),limits = c(0,32000)) +
  theme(legend.position='none',axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text = element_text(size=15, face="bold")
  )
ggsave("rarefraction_curve_shannon.pdf",height=10, width=8,dpi=600, device="pdf")
ggsave("rarefraction_curve_shannon.png",height=10, width=8,dpi=600, device="png")
ggsave("rarefraction_curve_shannon.svg",height=10, width=8,dpi=600, device="svg")


###########################################
library(vegan)
#### the txt filed need to remove the # created from biom header before imported into R
otu_table <- read.table("/home/hl46161/publish_living_mulch/exported-feature-table-for-NMDS/living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000-rarefied.txt",sep="\t",header=TRUE,row.names = 1,check.names = FALSE)
otu_table <- t(otu_table)
otu_table <- as.data.frame(otu_table)
otu_table$SampleID <- rownames(otu_table)
otu_table$SampleID

#import metadata file 
metadata <- read.table("/home/hl46161/publish_living_mulch/living_mulch_data_metafile_massalin.tsv",sep="\t",header=TRUE)
merge_file <- merge(metadata,otu_table,by="SampleID")

#subset samples based on sampling date 
August_samples <- merge_file[merge_file$Date_taken == "2018-08-31",]
June_samples <- merge_file[merge_file$Date_taken == "2018-06-28",]
May_samples <- merge_file[merge_file$Date_taken == "2018-05-21",]

#collect otu table 
August_abund <- August_samples[,c(21:ncol(August_samples))]
June_abund <- June_samples[,c(21:ncol(June_samples))]
May_abund <- May_samples[,c(21:ncol(May_samples))]

#collect soil environmental table 
August_soil_physical_value <- August_samples[,c(5:18)]
June_soil_physical_value <- June_samples[,c(5:18)]
May_soil_physical_value <- May_samples[,c(5:18)]

#############
## run NMDS based on bray curtis distance 
May_nmds = metaMDS(May_abund, distance = "bray")
May_nmds
## fit  enviromental varibales into NMDS
en = envfit(May_nmds,May_soil_physical_value, permutations = 999, na.rm = TRUE)
en
plot(May_nmds)
plot(en)

May_nmds$points

data.scores = as.data.frame(scores(May_nmds)) ## extract NMDS information 
data.scores$treatment = May_samples$Treatment ## add treatment inforamtion the sample

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #
en_coord_cont

env.scores <- as.data.frame(scores(en, display = "vectors")) #extracts relevant scores from envifit
env.scores
en$vectors$pvals

env.scores <- cbind(env.scores, pval = en$vectors$pvals) # add pvalues to dataframe
sig.env.scores <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

env.scores 
head(sig.env.scores)

may_gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = treatment, shape = treatment), size = 6, alpha = 0.5) + 
  scale_colour_manual(values = c("darkgreen", "steelblue","red","black"))  + ggtitle("May") +
  scale_shape_manual(values=c(15, 16, 17,18)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scores, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = sig.env.scores, aes(x = NMDS1, y = NMDS2), colour = "grey30",size=10, 
            fontface = "bold", label = row.names(sig.env.scores)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 30, colour = "grey30",face = "bold"),
        axis.title.x = element_text(size=30, face="bold"),
        axis.title.y = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"),
        axis.text.x = element_text(size=30, face="bold"),
        plot.title = element_text(size=30, face="bold")) + 
  labs(colour = "Treatment",shape ="Treatment")



may_gg

ggsave("/home/hl46161/publish_living_mulch/May_NMDS.pdf", height=8, width=8, device="pdf")

###############################################################################################
## run NMDS based on bray curtis distance 
June_nmds = metaMDS(June_abund, distance = "bray")
June_nmds
## fit  enviromental varibales into NMDS
en = envfit(June_nmds,June_soil_physical_value, permutations = 999, na.rm = TRUE)
en
plot(June_nmds)
plot(en)

June_nmds$points

data.scores = as.data.frame(scores(June_nmds)) ## extract NMDS information 
data.scores$treatment = June_samples$Treatment ## add treatment inforamtion the sample

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #
en_coord_cont

env.scores <- as.data.frame(scores(en, display = "vectors")) #extracts relevant scores from envifit
env.scores
en$vectors$pvals

env.scores <- cbind(env.scores, pval = en$vectors$pvals) # add pvalues to dataframe
sig.env.scores <- subset(env.scores, pval<=0.045) #subset data to show variables significant at 0.05

env.scores 
head(sig.env.scores)

June_gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = treatment, shape = treatment), size = 6, alpha = 0.5) + 
  
  scale_colour_manual(values = c("darkgreen", "steelblue","red","black"))  + 
  scale_shape_manual(values=c(15, 16, 17,18)) + ggtitle("June") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scores, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = sig.env.scores, aes(x = NMDS1, y = NMDS2), colour = "grey30",size=10, 
            fontface = "bold", label = row.names(sig.env.scores)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 30, colour = "grey30",face = "bold"),
        axis.title.x = element_text(size=30, face="bold"),
        axis.title.y = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"),
        axis.text.x = element_text(size=30, face="bold"),
        plot.title = element_text(size=30, face="bold"),
        ) + 
  labs(colour = "Treatment",shape ="Treatment")

June_gg

ggsave("/home/hl46161/publish_living_mulch/June_NMDS.pdf", height=8, width=8, device="pdf")

#,legend.position='none'

###########################################################################################
## run NMDS based on bray curtis distance 
August_nmds = metaMDS(August_abund, distance = "bray")
August_nmds
## fit  enviromental varibales into NMDS
en = envfit(August_nmds,August_soil_physical_value, permutations = 999, na.rm = TRUE)
en
plot(August_nmds)
plot(en)

August_nmds$points

data.scores = as.data.frame(scores(August_nmds)) ## extract NMDS information 
data.scores$treatment = August_samples$Treatment ## add treatment inforamtion the sample

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #
en_coord_cont

env.scores <- as.data.frame(scores(en, display = "vectors")) #extracts relevant scores from envifit
env.scores
en$vectors$pvals

env.scores <- cbind(env.scores, pval = en$vectors$pvals) # add pvalues to dataframe
sig.env.scores <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

env.scores 
head(sig.env.scores)

August_gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = treatment, shape = treatment), size = 6, alpha = 0.5) + 
  scale_colour_manual(values = c("darkgreen", "steelblue","red","black"))  + 
  scale_shape_manual(values=c(15, 16, 17,18)) + ggtitle("August") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scores, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = sig.env.scores, aes(x = NMDS1, y = NMDS2), colour = "grey30",size=10, 
            fontface = "bold", label = row.names(sig.env.scores)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 30, colour = "grey30",face = "bold"),
        axis.title.x = element_text(size=30, face="bold"),
        axis.title.y = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"),
        axis.text.x = element_text(size=30, face="bold"),
        plot.title = element_text(size=30, face="bold")) + 
  labs(colour = "Treatment",shape ="Treatment")

August_gg

ggsave("/home/hl46161/publish_living_mulch/August_NMDS.pdf", height=8, width=8, device="pdf")

########################################################
library(gridExtra)
library(ggpubr)

mega_NMDS = ggarrange(may_gg,June_gg,August_gg,
                      ncol = 3,
                      nrow = 1,
                      common.legend = TRUE)

#?ggarrange()
mega_NMDS

ggsave("/home/hl46161/publish_living_mulch/mega_NMDS.pdf", height=12, width=25, device="pdf")
ggsave("/home/hl46161/publish_living_mulch/mega_NMDS.png", height=12, width=25, device="png")
ggsave("/home/hl46161/publish_living_mulch/mega_NMDS.svg", height=12, width=25, device="svg")
ggsave("/home/hl46161/publish_living_mulch/mega_NMDS.tiff", height=12, width=25, device="tiff")

######################################

beta_diversity_graph <- function(input_pcoa,metadata,month){
  
  #extract cordinates from qza file 
  input_pcoa_cordination <- input_pcoa$data$Vectors
  input_pcoa_cordination
  #choose first two pc 
  input_pcoa_cordination <- input_pcoa_cordination[,c("SampleID","PC1","PC2")]
  input_pcoa_cordination_metadata <- dplyr::inner_join(input_pcoa_cordination,metadata,by = "SampleID")
  input_pcoa_cordination_metadata
  #exctract % of explaination 
  pc1 <- as.character(round(input_pcoa$data$ProportionExplained[1,1],3)*100)
  pc2 <- as.character(round(input_pcoa$data$ProportionExplained[1,2],3)*100)
  # generate x y axis labels 
  pc1_label = paste0("PC1","(",pc1,"%)")
  pc2_label = paste0("PC2","(",pc2, "%)")
  month = as.character(month)
  input_pcoa_cordination_metadata$Treatment <- factor(input_pcoa_cordination_metadata$Treatment, levels = c("NoCover","CerealRye","CrimsonClover","LivingMulch"))
  input_pcoa_plot <- ggplot(input_pcoa_cordination_metadata,aes(x=PC1, y=PC2, color=Treatment)) + 
    geom_point(aes(colour = Treatment), size = 6, alpha = 0.5) + ggtitle(month) +
    scale_colour_manual(values = c("black","darkgreen", "steelblue","red")) + xlab(pc1_label) + 
    ylab(pc2_label) + 
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), legend.key = element_blank(), 
          legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
          legend.text = element_text(size = 30, colour = "grey30",face = "bold"),
          axis.title.x = element_text(size=30, face="bold"),
          axis.title.y = element_text(size=30, face="bold"),
          plot.title = element_text(size=30, face="bold"),
          axis.text = element_blank()
    )
  return(input_pcoa_plot)
  }


test_gg <- beta_diversity_graph(mid_weighted_unifrac,metadata,"May")
test_gg


metadata <- read.table("living_mulch_data_metafile_massalin.tsv", sep="\t",header = TRUE)
metadata <- metadata[,c(1:3)]
mid_weighted_unifrac  <-read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/weighted_unifrac_pcoa_results.qza")
mid_unweighted_unifrac <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/unweighted_unifrac_pcoa_results.qza")
mid_bray_curtis <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/bray_curtis_pcoa_results.qza")
mid_aitchison <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000/June_aitchison-no-rarefraction-ordination.qza")


mid_weighted_unifrac_plot <- beta_diversity_graph(mid_weighted_unifrac,metadata,"June")
mid_weighted_unifrac_plot
ggsave(plot= mid_weighted_unifrac_plot,"/home/hl46161/publish_living_mulch/mid_weighted_unifrac_plot.png", height=8, width=8, device="png")
ggsave(plot= mid_weighted_unifrac_plot,"/home/hl46161/publish_living_mulch/mid_weighted_unifrac_plot.pdf", height=8, width=8, device="pdf")


mid_unweighted_unifrac_plot <- beta_diversity_graph(mid_unweighted_unifrac,metadata,"June")
ggsave(plot= mid_unweighted_unifrac_plot,"/home/hl46161/publish_living_mulch/mid_unweighted_unifrac_plot.png", height=8, width=8, device="png")
ggsave(plot= mid_unweighted_unifrac_plot,"/home/hl46161/publish_living_mulch/mid_unweighted_unifrac_plot.pdf", height=8, width=8, device="pdf")

mid_bray_curtis_plot <- beta_diversity_graph(mid_bray_curtis,metadata,"June")
ggsave(plot= mid_bray_curtis_plot,"/home/hl46161/publish_living_mulch/mid_bray_curtis_plot.png", height=8, width=8, device="png")
ggsave(plot= mid_bray_curtis_plot,"/home/hl46161/publish_living_mulch/mid_bray_curtis_plot.pdf", height=8, width=8, device="pdf")

mid_aitchison_plot <- beta_diversity_graph(mid_aitchison,metadata,"June")
ggsave(plot= mid_aitchison_plot,"/home/hl46161/publish_living_mulch/mid_aitchison_plot.png", height=8, width=8, device="png")
ggsave(plot= mid_aitchison_plot,"/home/hl46161/publish_living_mulch/mid_aitchison_plot.pdf", height=8, width=8, device="pdf")


##############################################################################

early_weighted_unifrac  <-read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/weighted_unifrac_pcoa_results.qza")
early_unweighted_unifrac <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/unweighted_unifrac_pcoa_results.qza")
early_bray_curtis <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/bray_curtis_pcoa_results.qza")
early_aitchison <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000/May_aitchison-no-rarefraction-ordination.qza")


early_weighted_unifrac_plot <- beta_diversity_graph(early_weighted_unifrac,metadata,"May")
early_weighted_unifrac_plot
ggsave(plot= early_weighted_unifrac_plot,"/home/hl46161/publish_living_mulch/early_weighted_unifrac_plot.png", height=8, width=8, device="png")
ggsave(plot= early_weighted_unifrac_plot,"/home/hl46161/publish_living_mulch/early_weighted_unifrac_plot.pdf", height=8, width=8, device="pdf")

early_weighted_unifrac_pcoa_cordination <- early_weighted_unifrac$data$Vectors
early_weighted_unifrac$data$Vectors

early_unweighted_unifrac_plot <- beta_diversity_graph(early_unweighted_unifrac,metadata,"May")
ggsave(plot= early_unweighted_unifrac_plot,"/home/hl46161/publish_living_mulch/early_unweighted_unifrac_plot.png", height=8, width=8, device="png")
ggsave(plot= early_unweighted_unifrac_plot,"/home/hl46161/publish_living_mulch/early_unweighted_unifrac_plot.pdf", height=8, width=8, device="pdf")

early_bray_curtis_plot <- beta_diversity_graph(early_bray_curtis,metadata,"May")
ggsave(plot= early_bray_curtis_plot,"/home/hl46161/publish_living_mulch/early_bray_curtis_plot.png", height=8, width=8, device="png")
ggsave(plot= early_bray_curtis_plot,"/home/hl46161/publish_living_mulch/early_bray_curtis_plot.pdf", height=8, width=8, device="pdf")

early_aitchison_plot <- beta_diversity_graph(early_aitchison,metadata,"May")
ggsave(plot= early_aitchison_plot,"/home/hl46161/publish_living_mulch/early_aitchison_plot.png", height=8, width=8, device="png")
ggsave(plot= early_aitchison_plot,"/home/hl46161/publish_living_mulch/early_aitchison_plot.pdf", height=8, width=8, device="pdf")


###############################################3

late_weighted_unifrac  <-read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/weighted_unifrac_pcoa_results.qza")
late_unweighted_unifrac <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/unweighted_unifrac_pcoa_results.qza")
late_bray_curtis <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/bray_curtis_pcoa_results.qza")
late_aitchison <- read_qza("/home/hl46161/publish_living_mulch/core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000/August_aitchison-no-rarefraction-ordination.qza")


late_weighted_unifrac_plot <- beta_diversity_graph(late_weighted_unifrac,metadata,"August")
late_weighted_unifrac_plot
ggsave(plot= late_weighted_unifrac_plot,"/home/hl46161/publish_living_mulch/late_weighted_unifrac_plot.png", height=8, width=8, device="png")
ggsave(plot= late_weighted_unifrac_plot,"/home/hl46161/publish_living_mulch/late_weighted_unifrac_plot.pdf", height=8, width=8, device="pdf")

late_unweighted_unifrac_plot <- beta_diversity_graph(late_unweighted_unifrac,metadata,"August")
ggsave(plot= late_unweighted_unifrac_plot,"/home/hl46161/publish_living_mulch/late_unweighted_unifrac_plot.png", height=8, width=8, device="png")
ggsave(plot= late_unweighted_unifrac_plot,"/home/hl46161/publish_living_mulch/late_unweighted_unifrac_plot.pdf", height=8, width=8, device="pdf")

late_bray_curtis_plot <- beta_diversity_graph(late_bray_curtis,metadata,"August")
ggsave(plot= late_bray_curtis_plot,"/home/hl46161/publish_living_mulch/late_bray_curtis_plot.png", height=8, width=8, device="png")
ggsave(plot= late_bray_curtis_plot,"/home/hl46161/publish_living_mulch/late_bray_curtis_plot.pdf", height=8, width=8, device="pdf")

late_aitchison_plot <- beta_diversity_graph(late_aitchison,metadata,"August")
ggsave(plot= late_aitchison_plot,"/home/hl46161/publish_living_mulch/late_aitchison_plot.png", height=8, width=8, device="png")
ggsave(plot= late_aitchison_plot,"/home/hl46161/publish_living_mulch/late_aitchison_plot.pdf", height=8, width=8, device="pdf")


#########################################################################

beta_diversity_plot = ggarrange(early_weighted_unifrac_plot,early_unweighted_unifrac_plot,early_bray_curtis_plot,early_aitchison_plot,
                      mid_weighted_unifrac_plot,mid_unweighted_unifrac_plot,mid_bray_curtis_plot,mid_aitchison_plot,
                      late_weighted_unifrac_plot,late_unweighted_unifrac_plot,late_bray_curtis_plot,late_aitchison_plot,
                      ncol = 4,
                      nrow = 4,
                      font.label = list(size = 30,face='bold'),
                      common.legend = TRUE)
beta_diversity_plot
ggsave("/home/hl46161/publish_living_mulch/beta_diversity_plot.pdf",beta_diversity_plot, height=20, width=20, device="pdf",dpi = 600)
ggsave("/home/hl46161/publish_living_mulch/beta_diversity_plot.png",beta_diversity_plot, height=20, width=20, device="png",dpi = 600)
ggsave("/home/hl46161/publish_living_mulch/beta_diversity_plot.svg",beta_diversity_plot, height=20, width=20, device="svg",dpi = 600)

##########################################################################
#differential abudance test 
library(DESeq2)

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
#rarecurve(t(otu_table(ps)), step=100, cex=0.5)
## save the rarefraction graph if necessary 
ggsave("rarefraction")

###################################################

ps.rarefied = rarefy_even_depth(ps,sample.size=19000, replace=F,rngseed=1)

#phylum_taxa_plot <- plot_bar(ps.rarefied,x="Treatment",fill="Phylum")
#phylum_taxa_plot

top20_taxa_list <- names(sort(taxa_sums(ps.rarefied), decreasing=TRUE)[1:30])
top20_taxa_list #shows 20 results

dat.aglo = tax_glom(ps.rarefied, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))

prune.dat.two = prune_taxa(top20_taxa_list, dat.trans)
dat.dataframe = psmelt(prune.dat.two)
dat.agr = aggregate(Abundance ~ Treatment + Date_taken + Phylum, data=dat.dataframe, FUN=mean)

aggregate


gg <- ggplot(dat.agr, aes(x=Treatment, y=Abundance, fill=Phylum)) + geom_bar(stat="identity") + 
  facet_grid(~Date_taken, scale="free") + 
  theme(
    legend.text = element_text(color = "black", size = 15),       #set the text size to be 20 and bold the text
    legend.title = element_text(color = "black", size = 15),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20, face="bold"),
    axis.text.x = element_text(face="bold"),
    strip.text.x = element_text(size = 10, face="bold"),
    panel.background = element_blank()
  ) + scale_x_discrete(labels=c("CR", "CC", "LM","NC"))

gg 


ggsave("/home/hl46161/publish_living_mulch/phylum_taxa_plot.pdf", height=10, width=12,dpi=600,device="pdf")
ggsave("/home/hl46161/publish_living_mulch/phylum_taxa_plot.svg", height=10, width=12,dpi=600,device="svg")


theme(
  legend.text = element_text(color = "black", size = 20),       #set the text size to be 20 and bold the text
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text.y = element_text(size=20, face="bold")
) + 
  theme(
    panel.grid.major = element_blank(),   ##remove grid line and background
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 20)) +
  geom_line() + 





####################################################


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


##first calculate overall differential abundant OTU between CrimsonClover and CerealRye 
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


##################################################


FamilyFiltered <- subset_taxa(ps, Family != "NA")
ps.family = tax_glom(FamilyFiltered, "Family")
otu_table(FamilyFiltered)
otu_table(ps.family)
family_otu_table = as(otu_table(ps.family), "matrix")
family_otu_table  = as.data.frame(family_otu_table)
#write.table(family_otu_table,"/home/hl46161/publish_replicate_new_living_mulch/family_otu_table.tsv",sep="\t")

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


write.csv(lm_vs_nc_filtered_1_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_VS_nc_differential_OTU.csv")


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


write.csv(lm_vs_cc_2_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_VS_cc_differential_OTU.csv")


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

write.csv(lm_vs_cr_filtered_4_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_VS_cr_differential_OTU.csv")


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

write.csv(cr_vs_nc__filtered_5_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_cr_VS_nc_differential_OTU.csv")


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

write.csv(cc_vs_nc__filtered_6_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_cc_VS_nc_differential_OTU.csv")


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

write.csv(cc_vs_cr__filtered_7_ps.family,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_cc_VS_cr_differential_OTU.csv")

##############################################################################################
library(dplyr)

#summarize differential abudant otu. First try to find out which dfOtu is shared between three treatment 
lm_vs_nc_filtered_1_ps.family$OTU_ID <- rownames(lm_vs_nc_filtered_1_ps.family)
cr_vs_nc__filtered_5_ps.family$OTU_ID <- rownames(cr_vs_nc__filtered_5_ps.family)
cc_vs_nc__filtered_6_ps.family$OTU_ID <- rownames(cc_vs_nc__filtered_6_ps.family)

lm_cr_cc_otu <- dplyr::inner_join(lm_vs_nc_filtered_1_ps.family,cc_vs_nc__filtered_6_ps.family,by = "OTU_ID")
lm_cr_cc_otu <- lm_cr_cc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,log2FoldChange.y,baseMean.x,padj.x)
lm_cr_cc_otu <- dplyr::inner_join(lm_cr_cc_otu, cr_vs_nc__filtered_5_ps.family,by = "OTU_ID")


lm_cr_cc_otu <- lm_cr_cc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,log2FoldChange.y,log2FoldChange,baseMean.x,padj.x)
colnames(lm_cr_cc_otu) <-c("Kingdom","Phylum","Class","Order","Family","OTU_ID","log2FoldChange.lm","log2FoldChange.cc","log2FoldChange.cr","baseMean","P_adjusted")

#round the log2fold change to two decimal point 
lm_cr_cc_otu$log2FoldChange.lm <- round(lm_cr_cc_otu$log2FoldChange.lm, digits = 2)
lm_cr_cc_otu$log2FoldChange.cc <- round(lm_cr_cc_otu$log2FoldChange.cc, digits = 2)
lm_cr_cc_otu$log2FoldChange.cr <- round(lm_cr_cc_otu$log2FoldChange.cr, digits = 2)
#lm_cr_cc_otu$P_adjusted <- round(lm_cr_cc_otu$P_adjusted, digits = 2)

#mannualy fill out the the missing taxonomic rank with taxaonomic information from upper taxonomic level 
lm_cr_cc_otu$Taxa <- as.character(lm_cr_cc_otu$Family)
lm_cr_cc_otu$Taxa[c(2,3,18,19,22,25,26,27,28,29,30,31,32,34,36,37,38,43)] <- c("c__Ktedonobacteria","o__Ktedonobacterales","p__Chloroflexi","o__Tepidisphaerales",
                                                                               " o__Rokubacteriales"," c__Planctomycetes"," c__Ktedonobacteria"," p__WPS-2","c__Phycisphaerae	",
                                                                               "c__Ktedonobacteria","p__Acidobacteriota"," c__Gammaproteobacteria","o__Acidobacteriales",
                                                                               "o__Gaiellales","c__Vicinamibacteria","p__Acidobacteriota","o__Elsterales",
                                                                               "c__Anaerolineae")

#generate the abudance label to indicate the abudance of otu 
lm_cr_cc_otu_table <- lm_cr_cc_otu %>%
  mutate(Abudance =
           case_when(baseMean <= 100 ~ "*",
                     baseMean > 1000 ~ "***",
                     baseMean >100 ~ "**",
           ))

lm_cr_cc_otu_table  <- lm_cr_cc_otu_table  %>% dplyr::select(log2FoldChange.lm:log2FoldChange.cr,Taxa,Abudance)




#output the shared OTU between three cover crop treatment
write.csv(lm_cr_cc_otu,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_lm_cc_cr_VS_nc_differential_OTU.csv",row.names = FALSE)

write.csv(lm_cr_cc_otu_table,"/home/hl46161/publish_living_mulch/exported-feature-table-for-deseq2/family_table_lm_cr_cc_vs_nc.csv",row.names = FALSE)

#create presence and absence matrix for differential abundant otu in each treatment combination 
lm_cr_cc_otu$CerealRye <- 1
lm_cr_cc_otu$LivingMulch <- 1
lm_cr_cc_otu$CrimsonClover <- 1

lm_cc_otu <- dplyr::inner_join(lm_vs_nc_filtered_1_ps.family, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
lm_cc_otu <- lm_cc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cc_otu <- dplyr::left_join(lm_cc_otu,cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
lm_cc_otu <- lm_cc_otu[rowSums(is.na(lm_cc_otu)) > 4,]
lm_cc_otu <- lm_cc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cc_otu$CerealRye <- 0
lm_cc_otu$LivingMulch <- 1
lm_cc_otu$CrimsonClover <- 1


lm_cr_otu <- dplyr::inner_join(lm_vs_nc_filtered_1_ps.family, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
lm_cr_otu <- lm_cr_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cr_otu <- dplyr::left_join(lm_cr_otu, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
lm_cr_otu <- lm_cr_otu[rowSums(is.na(lm_cr_otu)) > 4,]
lm_cr_otu <- lm_cr_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_cr_otu$CerealRye <- 1
lm_cr_otu$LivingMulch <- 1
lm_cr_otu$CrimsonClover <- 0

cc_cr_otu <- dplyr::inner_join(cc_vs_nc__filtered_6_ps.family, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
cc_cr_otu <- cc_cr_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cc_cr_otu <- dplyr::left_join(cc_cr_otu,lm_vs_nc_filtered_1_ps.family, by = "OTU_ID")
cc_cr_otu <- cc_cr_otu[rowSums(is.na(cc_cr_otu)) > 4,]
#cc_cr_otu <- cc_cr_otu %>% dplyr::select(Domain.x:Damily.x,OTU_ID)
cc_cr_otu <- cc_cr_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cc_cr_otu$CerealRye <- 1
cc_cr_otu$LivingMulch <- 0
cc_cr_otu$CrimsonClover <- 1

lm_nc_otu <- dplyr::left_join(lm_vs_nc_filtered_1_ps.family, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
lm_nc_otu <- lm_nc_otu[rowSums(is.na(lm_nc_otu)) > 7,]
lm_nc_otu <- lm_nc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
lm_nc_otu <- dplyr::left_join(lm_nc_otu, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
lm_nc_otu <- lm_nc_otu[rowSums(is.na(lm_nc_otu)) > 7,]
lm_nc_otu <- lm_nc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#lm_nc_otu <- lm_nc_otu [,c(7:ncol(lm_nc_otu))]
lm_nc_otu$CerealRye <- 0
lm_nc_otu$LivingMulch <- 1
lm_nc_otu$CrimsonClover <- 0

cc_nc_otu <- dplyr::left_join(cc_vs_nc__filtered_6_ps.family, lm_vs_nc_filtered_1_ps.family, by = "OTU_ID")
cc_nc_otu <- cc_nc_otu[rowSums(is.na(cc_nc_otu)) > 7,]
cc_nc_otu <- cc_nc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cc_nc_otu <- dplyr::left_join(cc_nc_otu, cr_vs_nc__filtered_5_ps.family, by = "OTU_ID")
cc_nc_otu <- cc_nc_otu[rowSums(is.na(cc_nc_otu)) > 7,]
cc_nc_otu <- cc_nc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#cc_nc_otu <- cc_nc_otu [,c(7:ncol(cc_nc_otu))]
cc_nc_otu$CerealRye <- 0
cc_nc_otu$LivingMulch <- 0
cc_nc_otu$CrimsonClover <- 1

cr_nc_otu <- dplyr::left_join(cr_vs_nc__filtered_5_ps.family, lm_vs_nc_filtered_1_ps.family, by = "OTU_ID")
cr_nc_otu <- cr_nc_otu[rowSums(is.na(cr_nc_otu)) > 7,]
cr_nc_otu <- cr_nc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
cr_nc_otu <- dplyr::left_join(cr_nc_otu, cc_vs_nc__filtered_6_ps.family, by = "OTU_ID")
cr_nc_otu <- cr_nc_otu[rowSums(is.na(cr_nc_otu)) > 7,]
cr_nc_otu <- cr_nc_otu %>% dplyr::select(Domain.x:Family.x,OTU_ID,log2FoldChange.x,baseMean.x)
#cr_nc_otu <- cr_nc_otu [,c(7:ncol(cr_nc_otu))]
cr_nc_otu$CerealRye <- 1
cr_nc_otu$LivingMulch <- 0
cr_nc_otu$CrimsonClover <- 0

#summary all differential otu patterin each treatment comparison 
otu_summary_matrix <- dplyr::bind_rows(lm_nc_otu,cc_nc_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,cr_nc_otu)

otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,lm_cr_cc_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,lm_cr_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,lm_cc_otu)
otu_summary_matrix <- dplyr::bind_rows(otu_summary_matrix,cc_cr_otu)
otu_summary_matrix$Confidence.x <- NULL

upset_data <- otu_summary_matrix %>% dplyr::select(CerealRye:CrimsonClover)

gg <- upset(upset_data,
      point.size = 3, line.size = 2,
      mainbar.y.label = "Distrubution of differential family", sets.x.label = "Total # of differential family",
      text.scale =c(2, 3, 2, 2, 2, 3)
)

gg

ggplot()
ggsave()

ggsave("differential_otu_upset_plot.svg", ggplotify::as.ggplot(gg),width = 15,height = 15,device = "svg",dpi = 600)

ggsave("differential_otu_upset_plot.pdf",width = 15,height = 15,device = "pdf",dpi = 600)
ggsave("differential_otu_upset_plot.png",width = 15,height = 15,device = "png",dpi = 600)
ggsave("differential_otu_upset_plot.svg",width = 15,height = 15,device = "svg",dpi = 600)


