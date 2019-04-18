#dependencies
library(ggplot2)
library(reshape2)
library(EnvStats)
library(ggpubr)
library(ggsignif)
library(psych)
library(aod)
library(survey)
library(lmtest)
library(pheatmap)
library(RColorBrewer)

#functions
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
removeITH <- function(x) {x$HM_status[x$HM_status == "ITH"] <- "NA"}

#data import
dirPath <- '~/Dropbox/Postdoc/MGMT Project/BSAS Results/BSAS Analysis R Project'
data <- read.csv(file.path(dirPath, "MGMT Project Data With Patient Annotations and Recurrences Added BSASmCG ITHmarked.csv"))
allCGs <- read.csv(file.path(dirPath, "allCGs.csv"), check.names=FALSE)
patientannotations <- read.csv(file.path(dirPath, "Matt_Patient_Annotations.csv"))
RNA <- read.csv(file.path(dirPath, "RNA_from_Saumya.csv"))


#data organization
data_identifiers <- data[c("Sample", "Patient", "HM_status", "HM_by_sample", "Tumor", "mol_type")]
allCGs <- merge(data_identifiers, allCGs, by="Sample")
BySample <- allCGs
BySample$PercentMethylation <- rowMeans(BySample[,7:33], na.rm = TRUE)
BySample_Primaries <- BySample[which(BySample$Tumor == "Primary"),]
BySample_Recurrences <- BySample[which(BySample$Tumor == "Recurrence"),]

##############################################################################################################################
#CpG Landscape Analysis
##############################################################################################################################

allCGs_reshape <- melt(allCGs, id.vars = c("Sample", "Patient", "HM_status", "HM_by_sample", "Tumor", "mol_type"))
colnames(allCGs_reshape)[which(names(allCGs_reshape) == "variable")] <- "CpG"
colnames(allCGs_reshape)[which(names(allCGs_reshape) == "value")] <- "Percent_Methylation"
allCGs_primaries <- allCGs_reshape[which(allCGs_reshape$Tumor == "Primary"),]
allCGs_recurrences <- allCGs_reshape[which(allCGs_reshape$Tumor == "Recurrence"),]

#Plot primaries by sample
ggplot(data=allCGs_primaries, aes(x=CpG, y=Percent_Methylation, group=Sample, color=HM_status)) +
  geom_line()+
  geom_point()+
  labs(x= "CpG sites", y = "Methylation (%)", title = "Percent Methylation for Primaries (by HM status at Recurrence)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#Plot recurrences by sample
ggplot(data=allCGs_recurrences, aes(x=CpG, y=Percent_Methylation, group=Sample, color=HM_by_sample)) +
  geom_line()+
  geom_point()+
  labs(x= "CpG sites", y = "Methylation (%)", title = "Percent Methylation for Recurrences (by HM status per sample)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#Heatmap for primaries by sample
ggplot(data=allCGs_primaries, aes(x=CpG, y=as.factor(Sample), fill=Percent_Methylation))+
  geom_tile()+
  scale_fill_gradient(low = "steelblue2", high = "indianred2")+
  facet_grid(HM_status~., scales="free", space="free")+  
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x= "CpG sites", y = "Sample", title = "Primaries")

#add annotation
ggplot(data=allCGs_primaries, aes(x=0, y=as.factor(Sample), fill=HM_status, group=Patient))+
  geom_tile()+
  geom_text(aes(label=Patient))+
  facet_grid(mol_type~., labeller = labeller(mol_type = c(A = "Astros", O = "Oligos")), scales="free", space="free")+  
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  labs(x= "HM", y = "Patient", title = "Primaries")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#Heatmap for recurrences by sample (HM by sample)
ggplot(data=allCGs_recurrences, aes(x=CpG, y=as.factor(Sample), fill=Percent_Methylation, group=Patient))+
  geom_tile()+
  scale_fill_gradient(low = "steelblue2", high = "indianred2")+
  facet_grid(HM_by_sample~., scales="free", space="free")+  
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x= "CpG sites", y = "Sample", title = "Recurrences")

#Heatmap for recurrences by sample (HM by patient)
ggplot(data=allCGs_recurrences, aes(x=CpG, y=as.factor(Sample), fill=Percent_Methylation, group=Patient))+
  geom_tile()+
  scale_fill_gradient(low = "steelblue2", high = "indianred2")+
  facet_grid(HM_status~., scales="free", space="free")+  
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x= "CpG sites", y = "Sample", title = "Recurrences")

#Aggregate sample data
agdata_allCGs <- aggregate(allCGs_reshape$Percent_Methylation, by=list(Patient=allCGs_reshape$Patient, Tumor=allCGs_reshape$Tumor, HM_status=allCGs_reshape$HM_status, mol_type = allCGs_reshape$mol_type, CpG = allCGs_reshape$CpG), mean, na.rm= TRUE)
colnames(agdata_allCGs)[which(names(agdata_allCGs) == "x")] <- "Percent_Methylation"
agdata_allCGs_primaries <- agdata_allCGs[which(agdata_allCGs$Tumor == "Primary"),]
agdata_allCGs_recurrences <- agdata_allCGs[which(agdata_allCGs$Tumor == "Recurrence"),]

#Plot primaries and Recurrences by patient
ggplot(data=agdata_allCGs, aes(x=CpG, y=Percent_Methylation, group=Patient, color=HM_status)) +
  geom_line()+
  geom_point()+
  labs(x= "CpG sites", y = "Methylation (%)", title = "Aggregated by Patient") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(Tumor~.)

#Heatmap for primaries by patient
ggplot(data=agdata_allCGs_primaries, aes(x=CpG, y=as.factor(Patient), fill=Percent_Methylation))+
  geom_tile()+
  scale_fill_gradient(low = "steelblue2", high = "indianred2")+
  facet_grid(mol_type~., labeller = labeller(mol_type = c(A = "Astros", O = "Oligos")), scales="free", space="free")+  
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x= "CpG sites (by distance from TSS)", y = "Patient", title = "Primaries")

#add annotation
ggplot(data=agdata_allCGs_primaries, aes(x=0, y=as.factor(Patient), fill=HM_status, group=Patient))+
  geom_tile()+
  facet_grid(mol_type~., labeller = labeller(mol_type = c(A = "Astros", O = "Oligos")), scales="free", space="free")+  
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  labs(x= "HM", y = "Patient", title = "Primaries")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#Heatmap for recurrences by patient
ggplot(data=agdata_allCGs_recurrences, aes(x=CpG, y=as.factor(Patient), fill=Percent_Methylation))+
  geom_tile()+
  scale_fill_gradient(low = "steelblue2", high = "indianred2")+
  facet_grid(mol_type~., labeller = labeller(mol_type = c(A = "Astros", O = "Oligos")), scales="free", space="free")+  
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x= "CpG sites (by distance from TSS)", y = "Patient", title = "Recurrences")

#Aggregate patient data to plot by HM status
agdata2_allCGs <- aggregate(agdata_allCGs$Percent_Methylation, by=list(Tumor=agdata_allCGs$Tumor, HM_status=agdata_allCGs$HM_status, CpG = agdata_allCGs$CpG), mean, na.rm= TRUE)
colnames(agdata2_allCGs)[which(names(agdata2_allCGs) == "x")] <- "Percent_Methylation"
agdata2_allCGs <- agdata2_allCGs[!(agdata2_allCGs$HM_status == "ITH"),]
agdata2_allCGs$TumorType = paste(agdata2_allCGs$Tumor, agdata2_allCGs$HM_status)

#Plot aggregate data by HM status
ggplot(data=agdata2_allCGs, aes(x=as.numeric.factor(CpG), y=Percent_Methylation, group=TumorType, color=Tumor)) +
  geom_line(aes(linetype=HM_status))+
  geom_point()+
  labs(x= NULL, y = "Methylation (%)", title = "Percent Methylation Across CpGs") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")+
  theme(legend.position="top")+ 
  theme(panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_continuous(limits= c(73,345), breaks = c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308), labels=c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308)) 
#facet_grid(Tumor~.)
#facet_grid(HM_status~.)

#Aggregate patient data to plot by mol_type
agdata3_allCGs <- aggregate(agdata_allCGs$Percent_Methylation, by=list(Tumor=agdata_allCGs$Tumor, mol_type=agdata_allCGs$mol_type, CpG = agdata_allCGs$CpG), mean, na.rm= TRUE)
colnames(agdata3_allCGs)[which(names(agdata3_allCGs) == "x")] <- "Percent_Methylation"
agdata3_allCGs$TumorType = paste(agdata3_allCGs$Tumor, agdata3_allCGs$mol_type)

#Plot aggregate data by mol_type
ggplot(data=agdata3_allCGs, aes(x=as.numeric.factor(CpG), y=Percent_Methylation, group=TumorType, color=mol_type)) +
  geom_line(aes(linetype=Tumor))+
  geom_point()+
  labs(x= NULL, y = "Methylation (%)", title = "Percent Methylation Across CpGs") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set2")+
  theme(legend.position="top")+ 
  theme(panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_continuous(limits= c(73,345), breaks = c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308), labels=c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308)) 

##############################################################################################################################
# By CpG Analysis
##############################################################################################################################

#Reshape ByCpG
ByCpG <- agdata2_allCGs
ByCpG$Tumor=NULL
ByCpG$HM_status=NULL
ByCpG <- reshape(ByCpG, idvar="CpG", timevar="TumorType", direction="wide")
colnames(ByCpG)[which(names(ByCpG)=="Percent_Methylation.Primary HM")] <- "Primary_HM"
colnames(ByCpG)[which(names(ByCpG)=="Percent_Methylation.Recurrence HM")] <- "Recurrence_HM"
colnames(ByCpG)[which(names(ByCpG)=="Percent_Methylation.Primary nonHM")] <- "Primary_nonHM"
colnames(ByCpG)[which(names(ByCpG)=="Percent_Methylation.Recurrence nonHM")] <- "Recurrence_nonHM"

AcrossCpG <- melt(ByCpG, id.vars="CpG")
colnames(AcrossCpG)[which(names(AcrossCpG) == "variable")] <- "TumorType"
colnames(AcrossCpG)[which(names(AcrossCpG) == "value")] <- "Percent_Methylation"

ggplot(data=AcrossCpG, aes(x=TumorType, y=Percent_Methylation, group=TumorType))+
  geom_boxplot(aes(color=TumorType))+
  geom_jitter(aes(color=TumorType))+
  ylim(0, 100)+
  stat_mean_sd_text(digits=2,y.pos=100)+
  labs(x=NULL)+
  geom_signif(test=wilcox.test, y_position=c(88,80,75,75), comparisons=list(c("Primary_HM", "Primary_nonHM"),c("Recurrence_HM", "Recurrence_nonHM"),c("Primary_HM", "Recurrence_HM"),c("Primary_nonHM", "Recurrence_nonHM")))

ggplot(data=AcrossCpG, aes(x=CpG, y=Percent_Methylation, group=TumorType, color=TumorType))+
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

wilcox.test(ByCpG$Primary_HM, ByCpG$Primary_nonHM, paired=TRUE, alternative = "two.sided")
wilcox.test(ByCpG$Primary_HM, ByCpG$Recurrence_HM, paired=TRUE, alternative = "two.sided")
wilcox.test(ByCpG$Primary_nonHM, ByCpG$Recurrence_nonHM, paired=TRUE, alternative = "two.sided")
wilcox.test(ByCpG$Recurrence_HM, ByCpG$Recurrence_nonHM, paired=TRUE, alternative = "two.sided")

ByCpG$DiffPrimaries <- ByCpG$Primary_HM - ByCpG$Primary_nonHM
ByCpG$DiffRecurrences <- ByCpG$Recurrence_HM - ByCpG$Recurrence_nonHM
ByCpG$DiffHM <- ByCpG$Recurrence_HM - ByCpG$Primary_HM
ByCpG$DiffnonHM <- ByCpG$Recurrence_nonHM - ByCpG$Primary_nonHM
CpG_Diffs <- ByCpG 
CpG_Diffs$Primary_HM = NULL
CpG_Diffs$Primary_nonHM = NULL
CpG_Diffs$Recurrence_HM = NULL
CpG_Diffs$Recurrence_nonHM = NULL
CpG_TumorDiffs <- CpG_Diffs
CpG_TumorDiffs$DiffHM = NULL
CpG_TumorDiffs$DiffnonHM = NULL
CpG_HMStatusDiffs <- CpG_Diffs
CpG_HMStatusDiffs$DiffPrimaries <- NULL
CpG_HMStatusDiffs$DiffRecurrences <- NULL
CpG_TumorDiffs <- melt(CpG_TumorDiffs, id.vars="CpG")
colnames(CpG_TumorDiffs)[which(names(CpG_TumorDiffs) == "variable")] <- "Tumor"
colnames(CpG_TumorDiffs)[which(names(CpG_TumorDiffs) == "value")] <- "Percent_Methylation"
CpG_HMStatusDiffs <- melt(CpG_HMStatusDiffs, id.vars="CpG")
colnames(CpG_HMStatusDiffs)[which(names(CpG_HMStatusDiffs) == "variable")] <- "HM_Status"
colnames(CpG_HMStatusDiffs)[which(names(CpG_HMStatusDiffs) == "value")] <- "Percent_Methylation"
CpG_AllDiffs <- melt(CpG_Diffs, id.vars = "CpG")

#Plot Differences in Percent Methylation (across tumor type)
ggplot(data=CpG_TumorDiffs, aes(x=CpG, y=Percent_Methylation, group=Tumor)) +
  geom_line(aes(linetype=Tumor))+
  geom_point()+
  labs(x= "CpG sites", y = "Methylation (%)", title = "Differences in Percent Methylation (HM - nonHM)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")

#Plot Differences in Percent Methylation (across HM_status)
ggplot(data=CpG_HMStatusDiffs, aes(x=CpG, y=Percent_Methylation, group=HM_Status, color=HM_Status)) +
  geom_line()+
  geom_point()+
  labs(x= "CpG sites", y = "Methylation (%)", title = "Differences in Percent Methylation (Recurrence - Primary)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")

#Plot All Differences in Percent Methylation
ggplot(data=CpG_AllDiffs, aes(x=CpG, y=abs(value), group=variable, color=variable)) +
  geom_line()+
  geom_point()+
  labs(x= "CpG sites", y = "Methylation (%)", title = "Differences in Percent Methylation (Absolute Values)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")


##############################################################################################################################
# Intratumoral heterogeneity
##############################################################################################################################

# Boxplot of all samples showing ITH
ggplot(BySample, aes(x=factor(Patient), y=PercentMethylation)) +
  geom_boxplot(aes(colour=HM_status)) +
  facet_grid(Tumor ~ .) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Methylation (%)", title = "Methylation Across Patients") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(BySample[which(BySample$HM_status == "ITH"),], aes(x=Tumor, y=PercentMethylation, color=Tumor)) +
  geom_jitter(aes(shape=HM_by_sample), size=3)+ 
  scale_shape_manual(values=c(1,13))+
  facet_wrap(Patient~.)+
  labs(x= NULL, y = "Methylation (%)", title = "Recurrences with Intratumoral Heteogeneity")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")

##############################################################################################################################
# By Patient Analysis
##############################################################################################################################

ByPatient <- aggregate(agdata_allCGs$Percent_Methylation, by=list(Patient = agdata_allCGs$Patient, Tumor=agdata_allCGs$Tumor, mol_type=agdata_allCGs$mol_type, HM_status=agdata_allCGs$HM_status), mean, na.rm= TRUE)
colnames(ByPatient)[which(names(ByPatient) == "x")] <- "Percent_Methylation"
ByPatient$TumorType <- paste0(ByPatient$Tumor, "_", ByPatient$HM_status)

#Boxplots by tumor type 
ggplot(data=ByPatient, aes(x=TumorType, y=Percent_Methylation, group=TumorType))+
  geom_boxplot(aes(color=Tumor))+
  geom_jitter(aes(color=Tumor, shape=mol_type))+
  ylim(0, 100)+
  stat_mean_sd_text(digits=2,y.pos=100)+
  labs(x=NULL)+
  scale_color_brewer(palette="Set1")+
  geom_signif(test=wilcox.test, y_position=c(80,80,88), comparisons=list(c("Primary_HM", "Primary_nonHM"),c("Recurrence_HM", "Recurrence_nonHM"),c("Primary_nonHM", "Recurrence_nonHM")))

#Boxplots by tumor type (no ITH)
ByPatientNoITH <- ByPatient
ByPatientNoITH$HM_status[ByPatientNoITH$HM_status == "ITH"] <- "HM"
ByPatientNoITH$TumorType <- paste0(ByPatientNoITH$Tumor, "_", ByPatientNoITH$HM_status)
ggplot(data=ByPatientNoITH, aes(x=TumorType, y=Percent_Methylation, group=TumorType))+
  geom_boxplot(aes(color=Tumor))+
  scale_color_brewer(palette="Set1")+
  geom_jitter(aes(color=Tumor, shape=mol_type))+
  ylim(0, 100)+
  stat_mean_sd_text(digits=2,y.pos=100)+
  stat_n_text()+
  labs(x=NULL)+
  geom_signif(test=wilcox.test, y_position=c(78,78,83,90), comparisons=list(c("Primary_HM", "Primary_nonHM"),c("Recurrence_HM", "Recurrence_nonHM"),c("Primary_HM", "Recurrence_HM"),c("Primary_nonHM", "Recurrence_nonHM")))

#Boxplots by molecular type 
ggplot(data=ByPatient, aes(x=mol_type, y=Percent_Methylation))+
  geom_boxplot()+
  geom_jitter(aes(color=Tumor, shape=mol_type))+
  ylim(0, 100)+
  stat_mean_sd_text(digits=2,y.pos=0)+
  labs(x=NULL)+
  #facet_grid(.~Tumor)+
  stat_compare_means()+
  scale_color_brewer(palette="Set1")

#Boxplots by Tumor
ggplot(data=ByPatient, aes(x=Tumor, y=Percent_Methylation))+
  geom_boxplot()+
  geom_jitter(aes(color=Tumor, shape=mol_type))+
  ylim(0, 100)+
  stat_mean_sd_text(digits=2,y.pos=0)+
  labs(x=NULL)+
  stat_compare_means()+
  scale_color_brewer(palette="Set1")


##############################################################################################################################
# Patient Cohort
##############################################################################################################################
Primaries <- ByPatient[which(ByPatient$Tumor == 'Primary'),]
Recurrences <- ByPatient[which(ByPatient$Tumor == 'Recurrence'),]

# Patient cohort
ggplot(data=ByPatient, aes(x=mol_type))+
  geom_bar(stat="count", aes(fill=HM_status))+
  facet_grid(.~Tumor, labeller = labeller(mol_type = c(A = "Astros", O = "Oligos")))+
  scale_color_brewer(palette="Dark2")+
  geom_text(stat="count", aes(label = ..count..), vjust=-0.5)

# Primaries
ggplot(data=Primaries, aes(x=mol_type))+
  geom_bar(stat="count", aes(fill=HM_status))+
  scale_color_brewer(palette="Dark2")+
  geom_text(stat="count", aes(label = ..count..), vjust=-0.5)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24))+
  labs(x="Molecular Subtype", y="Patients")

# Samples included in methylation analysis
ggplot(data=allCGs, aes(x=as.factor(Patient)))+
  geom_bar(stat="count", aes(fill=Tumor))+
  scale_fill_brewer(palette="Set1")+
  geom_text(stat="count", aes(label = ..count..), vjust=-0.5)+
  labs(x="Patient ID", y="Samples For Methylation Analysis")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11))+
  facet_grid(.~mol_type, scales="free", space="free", labeller = labeller(mol_type = c(A = "Astros", O = "Oligos")))+
  theme(panel.grid.minor = element_blank())

# Determination of HM status
ggplot(data=BySample_Recurrences, aes(x=as.factor(Patient)))+
  geom_bar(stat="count", aes(fill=HM_by_sample))+
  scale_fill_brewer(palette="Set2")+
  geom_text(stat="count", aes(label = ..count..), vjust=-0.5)+
  labs(x="Patient ID", y="Recurrent Samples")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11))+
  facet_grid(.~mol_type, scales="free", space="free", labeller = labeller(mol_type = c(A = "Astros", O = "Oligos")))+
  theme(panel.grid.minor = element_blank())

##############################################################################################################################
# Paired data
##############################################################################################################################

PairedData <- merge(Primaries, Recurrences, by = c("Patient", "HM_status", "mol_type"))
PairedData$Tumor.x <- NULL
PairedData$Tumor.y <- NULL
PairedData$TumorType.x <- NULL
PairedData$TumorType.y <- NULL
colnames(PairedData)[which(names(PairedData) == "Percent_Methylation.x")] <- "Primary"
colnames(PairedData)[which(names(PairedData) == "Percent_Methylation.y")] <- "Recurrence"

#Paired plots
ggpaired(data=PairedData, cond1="Primary", cond2="Recurrence", id= "Patient", fill="condition", ylab= "Percent Methylation") 
ggpaired(data=PairedData, cond1="Primary", cond2="Recurrence", id= "Patient", fill="condition", facet.by = "HM_status", ylab= "Percent Methylation") 
ggpaired(data=PairedData, cond1="Primary", cond2="Recurrence", id= "Patient", fill="condition", facet.by = "mol_type", ylab= "Percent Methylation") 

ggpaired(data=PairedData[which(PairedData$mol_type == "A"),], cond1="Primary", cond2="Recurrence", id= "Patient", fill="condition", ylab= "Percent Methylation") 
ggpaired(data=PairedData[which(PairedData$mol_type == "A"),], cond1="Primary", cond2="Recurrence", id= "Patient", fill="condition", facet.by = "HM_status", ylab= "Percent Methylation") 
ggpaired(data=PairedData[which(PairedData$mol_type == "O"),], cond1="Primary", cond2="Recurrence", id= "Patient", fill="condition", facet.by = "HM_status", ylab= "Percent Methylation") 

##############################################################################################################################
# Add patient annotations
##############################################################################################################################

AnnotatedPatients <- merge(ByPatient, patientannotations, by="Patient", all.x=TRUE)
AnnotatedPatientsNoITH <- AnnotatedPatients
AnnotatedPatientsNoITH$HM_status[AnnotatedPatientsNoITH$HM_status == "ITH"] <- "HM"
AnnotatedPatientsNoITH$TumorType <- paste0(AnnotatedPatientsNoITH$Tumor, "_", AnnotatedPatientsNoITH$HM_status)

# Association of HM with grade of recurrence
ggplot(data=AnnotatedPatients[which(AnnotatedPatients$Tumor == "Primary"),], aes(x=as.factor(Grade.at.recurrence)))+
  geom_bar(stat="count", aes(fill=HM_status))+
  scale_color_brewer(palette="Dark2")+
  geom_text(stat="count", aes(label = ..count..), vjust=-0.5)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24))+
  labs(x="Grade at recurrence", y="Patients")+
  facet_grid(.~mol_type, scales="free", space="free", labeller = labeller(mol_type = c(A = "Astros", O = "Oligos")))+
  theme(panel.grid.minor = element_blank())

ggplot(data=AnnotatedPatientsNoITH, aes(x=TumorType, y=Percent_Methylation, group=TumorType))+
  geom_boxplot(aes(color=TumorType))+
  geom_jitter(aes(color=TumorType, shape=mol_type))+
  ylim(0, 100)+
  stat_mean_sd_text(digits=2,y.pos=100)+
  stat_n_text()+
  labs(x=NULL)+
  facet_wrap(Grade.at.recurrence~.)+
  geom_signif(test=wilcox.test, y_position=c(78,78,83,90), comparisons=list(c("Primary_HM", "Primary_nonHM"),c("Recurrence_HM", "Recurrence_nonHM"),c("Primary_HM", "Recurrence_HM"),c("Primary_nonHM", "Recurrence_nonHM")))

ggplot(data=AnnotatedPatientsNoITH, aes(x=TumorType, y=Percent_Methylation, group=TumorType))+
  geom_boxplot(aes(color=TumorType))+
  geom_jitter(aes(color=TumorType, shape=mol_type))+
  ylim(0, 100)+
  stat_mean_sd_text(digits=2,y.pos=100)+
  stat_n_text()+
  labs(x=NULL)+
  #facet_wrap(mol_type~.)+
  geom_signif(test=wilcox.test, y_position=c(78,78,83,90), comparisons=list(c("Primary_HM", "Primary_nonHM"),c("Recurrence_HM", "Recurrence_nonHM"),c("Primary_HM", "Recurrence_HM"),c("Primary_nonHM", "Recurrence_nonHM")))

ggplot(data=ByPatientNoITH, aes(x=TumorType, y=Percent_Methylation, group=TumorType))+
  geom_boxplot(aes(color=TumorType))+
  geom_jitter(aes(color=TumorType, shape=mol_type))+
  ylim(0, 100)+
  stat_mean_sd_text(digits=2,y.pos=100)+
  stat_n_text()+
  labs(x=NULL)+
  #facet_wrap(mol_type~.)+
  geom_signif(test=wilcox.test, y_position=c(78,78,83,90), comparisons=list(c("Primary_HM", "Primary_nonHM"),c("Recurrence_HM", "Recurrence_nonHM"),c("Primary_HM", "Recurrence_HM"),c("Primary_nonHM", "Recurrence_nonHM")))

##############################################################################################################################
# Direct Comparison with Hinke
##############################################################################################################################

Hinke=list(1,5,10,18,21,4,6,11,17)
HinkePatients <- AnnotatedPatientsNoITH[which(AnnotatedPatientsNoITH$Patient %in% Hinke),]
ggplot(data=HinkePatients, aes(x=TumorType, y=Percent_Methylation, group=TumorType))+
  geom_boxplot(aes(color=TumorType))+
  geom_jitter(aes(color=TumorType, shape=mol_type))+
  ylim(0, 100)+
  stat_mean_sd_text(digits=2,y.pos=100)+
  labs(x=NULL)+
  geom_signif(test=wilcox.test, y_position=c(80,80,88), comparisons=list(c("Primary_HM", "Primary_nonHM"),c("Recurrence_HM", "Recurrence_nonHM"),c("Primary_nonHM", "Recurrence_nonHM")))

# Grade of recurrence plot as in Hinke's paper
ggplot(data=HinkePatients[which(HinkePatients$Tumor == "Primary"),], aes(x=HM_status))+
  geom_bar(stat="count", aes(fill=as.factor(Grade.at.recurrence)))+
  scale_color_brewer(palette="Dark2")+
  geom_text(stat="count", aes(label = ..count..), vjust=-0.5)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24))+
  labs(x="Grade at recurrence", y="Patients")+
  facet_grid(.~mol_type, scales="free", space="free", labeller = labeller(mol_type = c(A = "Astros", O = "Oligos")))+
  theme(panel.grid.minor = element_blank())

##############################################################################################################################
# Comparison with Infinium
##############################################################################################################################

ggplot(data, aes(x=MGMT_BSAS_cg12981137, y=MGMT_Infinium_cg12981137, color=Tumor)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "M Value (Infinium Array)", title = "Infinium v. BSAS at CpG site 174") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")
cor.test(data$MGMT_BSAS_cg12981137, data$MGMT_Infinium_cg12981137, method = "pearson", conf.level = 0.95)
cor.test(data$MGMT_BSAS_cg12981137, data$MGMT_Infinium_cg12981137, method = "spearman", conf.level = 0.95)
                   
##############################################################################################################################
# Comparison with RNA data
##############################################################################################################################
RNAdata <- merge(BySample, RNA, by="Sample", all.x=TRUE)

#Correlation matrix
RNACorr <- corr.test(RNAdata[,7:34], RNAdata[,35], use="pairwise", method="pearson")
RNACorrMatrix <- as.data.frame(RNACorr$r)
RNACorrMatrix$CpG <- rownames(RNACorrMatrix)
RNACorrMatrix$p <- RNACorr$p
colnames(RNACorrMatrix)[which(names(RNACorrMatrix) == "V1")] <- "PearsonCorr"
RNACorrCpG <- RNACorrMatrix[-c(28),]

#Plot Pearson Corr
ggplot(data=RNACorrCpG, aes(x=as.numeric(CpG), y=PearsonCorr)) +
  geom_line()+
  geom_point()+
  labs(x= NULL, y = "Pearson Corr", title = "Correlation with MGMT RNA expression") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")+
  theme(legend.position="top")+ 
  theme(panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_continuous(limits= c(73,345), breaks = c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308), labels=c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308)) 

#Plot p-values
ggplot(data=RNACorrCpG, aes(x=as.numeric(CpG), y=p)) +
  geom_line()+
  geom_point()+
  labs(x= NULL, y = "Adjusted P values", title = "P values adjusted for multiple tests") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")+
  theme(legend.position="top")+ 
  theme(panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_continuous(limits= c(73,345), breaks = c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308), labels=c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308)) 

#plot average percent methylation across CpGs
ggplot(RNAdata, aes(x=PercentMethylation, y=MGMT_RNA2, color=Tumor)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "RNA log count", title = "Across CpGs") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  #geom_smooth(method=lm)+
  scale_color_brewer(palette="Set1")
cor.test(RNAdata$PercentMethylation, RNAdata$MGMT_RNA2, method = "pearson", conf.level = 0.95)

ggplot(RNAdata, aes(x=RNAdata[,9], y=MGMT_RNA2, color=Tumor)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "RNA log count", title = "CpG site 106") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  #geom_smooth(method=lm)+
  scale_color_brewer(palette="Set1")
cor.test(RNAdata[,9], RNAdata$MGMT_RNA2, method = "pearson", conf.level = 0.95)

ggplot(RNAdata, aes(x=RNAdata[,10], y=MGMT_RNA2, color=Tumor)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "RNA log count", title = "CpG site 113") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  #geom_smooth(method=lm)+
  scale_color_brewer(palette="Set1")
cor.test(RNAdata[,10], RNAdata$MGMT_RNA2, method = "pearson", conf.level = 0.95)

ggplot(RNAdata, aes(x=RNAdata[,16], y=MGMT_RNA2, color=Tumor)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "RNA log count", title = "CpG site 142") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  #geom_smooth(method=lm)+
  scale_color_brewer(palette="Set1")
cor.test(RNAdata[,16], RNAdata$MGMT_RNA2, method = "pearson", conf.level = 0.95)

##############################################################################################################################
# Logistic Regression Model
##############################################################################################################################
BySample_PrimariesNoITH <- BySample_Primaries
BySample_PrimariesNoITH$HM_status[BySample_PrimariesNoITH$HM_status == "ITH"] <- "HM"
BySample_PrimariesNoITHAnnotated <- merge(BySample_PrimariesNoITH, patientannotations, by= "Patient", all.x=TRUE)
BySample_PrimariesNoITHAnnotated$HM_status <- relevel(BySample_PrimariesNoITHAnnotated$HM_status, ref = 'nonHM')
Primaries$HM_status <- relevel(Primaries$HM_status, ref = 'nonHM')

#Average Percent Methylation
mylogit <- glm(HM_status ~ PercentMethylation + mol_type + factor(Grade.at.recurrence), data = BySample_PrimariesNoITHAnnotated, family = "binomial")
summary(mylogit)
exp(coef(mylogit))

wald.test(b=coef(mylogit), Sigma=vcov(mylogit), Terms = 4:5)
regterm

#Average Percent Methylation by Patient
mylogit <- glm(HM_status ~ Percent_Methylation + mol_type, data = Primaries, family = "binomial")
summary(mylogit)
exp(coef(mylogit))

wald.test(b=coef(mylogit), Sigma=vcov(mylogit), Terms = 4:5)
regterm

#delete
wilcox.test(BySample_PrimariesNoITH$`142`~BySample_PrimariesNoITH$HM_status)
wilcox.test(BySample_PrimariesNoITH$`113`~BySample_PrimariesNoITH$HM_status)
wilcox.test(BySample_PrimariesNoITH$`106`~BySample_PrimariesNoITH$HM_status)
wilcox.test(BySample_PrimariesNoITH$PercentMethylation~BySample_PrimariesNoITH$HM_status)

ByPatient_PercentMethylation <- aggregate(BySample$PercentMethylation, by=list(Patient = BySample$Patient,Tumor=BySample$Tumor, HM_status=BySample$HM_status, mol_type=BySample$mol_type), mean, na.rm= TRUE)
colnames(ByPatient_PercentMethylation)[which(names(ByPatient_PercentMethylation) == "x")] <- "PercentMethylation"
ByPatient_106 <- aggregate(BySample$`106`, by=list(Patient = BySample$Patient,Tumor=BySample$Tumor, HM_status=BySample$HM_status, mol_type=BySample$mol_type), mean, na.rm= TRUE)
colnames(ByPatient_PercentMethylation)[which(names(ByPatient_PercentMethylation) == "x")] <- "CpG106"
ByPatient_113 <- aggregate(BySample$`113`, by=list(Patient = BySample$Patient,Tumor=BySample$Tumor, HM_status=BySample$HM_status, mol_type=BySample$mol_type), mean, na.rm= TRUE)
colnames(ByPatient_PercentMethylation)[which(names(ByPatient_PercentMethylation) == "x")] <- "CpG113"
ByPatient_142 <- aggregate(BySample$`142`, by=list(Patient = BySample$Patient,Tumor=BySample$Tumor, HM_status=BySample$HM_status, mol_type=BySample$mol_type), mean, na.rm= TRUE)
colnames(ByPatient_PercentMethylation)[which(names(ByPatient_PercentMethylation) == "x")] <- "CpG142"

##############################################################################################################################
# 'Pretty' heatmaps
##############################################################################################################################

# By Sample | Recurrences
annotation_row = BySample_Recurrences[c("mol_type", "HM_by_sample")]
colnames(annotation_row)[which(names(annotation_row) == "mol_type")] <- "Molecular Subtype"
colnames(annotation_row)[which(names(annotation_row) == "HM_by_sample")] <- "Hypermutation Status"
labels_row = as.factor(BySample_Recurrences$Patient)
pheatmap(BySample_Recurrences[,7:33], cluster_cols=FALSE, annotation_row = annotation_row, labels_row = labels_row, annotation_names_row=F, main = "All Recurrence Samples")

BySample_Recurrences_concordant <- BySample_Recurrences[c("106", "113", "142")]
pheatmap(BySample_Recurrences_concordant, cluster_cols=FALSE, annotation_row = annotation_row, labels_row = labels_row, annotation_names_row=F, main = "All Recurrence Samples")

# By Sample | Primaries (with ITH)
annotation_row = BySample_Primaries[c("mol_type", "HM_status")]
colnames(annotation_row)[which(names(annotation_row) == "mol_type")] <- "Molecular Subtype"
colnames(annotation_row)[which(names(annotation_row) == "HM_by_sample")] <- "Hypermutation Status At Recurrence"
labels_row = as.factor(BySample_Primaries$Patient)
pheatmap(BySample_Primaries[,7:33], cluster_cols=FALSE, annotation_row = annotation_row, labels_row = labels_row, annotation_names_row=F, main = "All Primary Samples")

BySample_Primaries_concordant <- BySample_Primaries[c("106", "113", "142")]
pheatmap(BySample_Primaries_concordant, cluster_cols=FALSE, annotation_row = annotation_row, labels_row = labels_row, annotation_names_row=F, main = "All Primary Samples")

# By Sample | Primaries (no ITH)
annotation_row = BySample_PrimariesNoITH[c("mol_type", "HM_status")]
colnames(annotation_row)[which(names(annotation_row) == "mol_type")] <- "Molecular Subtype"
colnames(annotation_row)[which(names(annotation_row) == "HM_by_sample")] <- "Hypermutation Status At Recurrence"
labels_row = as.factor(BySample_PrimariesNoITH$Patient)
pheatmap(BySample_PrimariesNoITH[,7:33], cluster_cols=FALSE, annotation_row = annotation_row, labels_row = labels_row, annotation_names_row=F, main = "All Primary Samples")

BySample_PrimariesNoITH_concordant <- BySample_PrimariesNoITH[c("106", "113", "142")]
pheatmap(BySample_PrimariesNoITH_concordant, cluster_cols=FALSE, annotation_row = annotation_row, labels_row = labels_row, annotation_names_row=F, main = "All Primary Samples")

