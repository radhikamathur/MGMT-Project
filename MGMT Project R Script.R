#dependencies
library(ggplot2)
library(reshape2)

#data import
dirPath <- '~/Dropbox/Postdoc/MGMT Project/BSAS Results/BSAS Analysis R Project'
data <- read.csv(file.path(dirPath, "MGMT_Infinium_BSAS_RNA.csv"))

#define subsets
data_primary <- data[which(data$Tumor == 'Primary'),]
data_recurrence <- data[which(data$Tumor == 'Recurrence'),]
data_HM <- data[which(data$HM_status == 'HM'),]
data_nonHM <- data[which(data$HM_status == 'nonHM'),]

#scatterplot Infinium v.BSAS 
ggplot(data, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "All individual samples: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data$MGMT_BSAS, data$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data$MGMT_BSAS, data$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot Infinium v. RNA
ggplot(data, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count", title = "All individual samples: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data$MGMT_RNA, data$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data$MGMT_RNA, data$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot Infinium v. RNA, facet by Tumor
ggplot(data, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count",  title = "All individual samples: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  facet_grid(. ~ Tumor, scales = "free", space = "free")
cor.test(data_primary$MGMT_RNA, data_primary$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_primary$MGMT_RNA, data_primary$MGMT_Infinium, method = "spearman", conf.level = 0.95)
cor.test(data_recurrence$MGMT_RNA, data_recurrence$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_recurrence$MGMT_RNA, data_recurrence$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot BSAS v. RNA
ggplot(data, aes(x=MGMT_BSAS, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Percent CpG Methylation (BSAS)", y = "RNA log count", title = "All individual samples: BSAS v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data$MGMT_RNA, data$MGMT_BSAS, method = "pearson", conf.level = 0.95)
cor.test(data$MGMT_RNA, data$MGMT_BSAS, method = "spearman", conf.level = 0.95)

#aggregate data
agdata_Infinium <- aggregate(data$MGMT_Infinium, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean, na.rm= TRUE)
agdata_BSAS <-aggregate(data$MGMT_BSAS, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean, na.rm= TRUE)
agdata_RNA <-aggregate(data$MGMT_RNA, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean, na.rm= TRUE)
agdata <- merge(agdata_Infinium, agdata_BSAS, by = c("Patient", "Tumor", "HM_status"))
agdata <- merge(agdata, agdata_RNA, by = c("Patient", "Tumor", "HM_status"))
colnames(agdata)[which(names(agdata) == "x.x")] <- "MGMT_Infinium"
colnames(agdata)[which(names(agdata) == "x.y")] <- "MGMT_BSAS"
colnames(agdata)[which(names(agdata) == "x")] <- "MGMT_RNA"

#define subsets
agdata_primary <- agdata[which(agdata$Tumor == 'Primary'),]
agdata_recurrence <- agdata[which(agdata$Tumor == 'Recurrence'),]
agdata_HM <- agdata[which(agdata$HM_status == 'HM'),]
agdata_nonHM <- agdata[which(agdata$HM_status == 'nonHM'),]

#scatterplot of aggregate data - Infinium v. BSAS
ggplot(agdata, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "All patients: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata$MGMT_BSAS, agdata$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata$MGMT_BSAS, agdata$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. RNA
ggplot(agdata, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count", title = "All Patients: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata$MGMT_RNA, agdata$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata$MGMT_RNA, agdata$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. RNA, facet by tumor
ggplot(agdata, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count", title = "All Patients: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  facet_grid(. ~ Tumor, scales = "free", space = "free")
cor.test(agdata_primary$MGMT_RNA, agdata_primary$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_primary$MGMT_RNA, agdata_primary$MGMT_Infinium, method = "spearman", conf.level = 0.95)
cor.test(agdata_recurrence$MGMT_RNA, agdata_recurrence$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_recurrence$MGMT_RNA, agdata_recurrence$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - BSAS v. RNA
ggplot(agdata, aes(x=MGMT_BSAS, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Percent CpG Methylation (BSAS)", y = "RNA log count", title = "All Patients: BSAS v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata$MGMT_RNA, agdata$MGMT_BSAS, method = "pearson", conf.level = 0.95)
cor.test(agdata$MGMT_RNA, agdata$MGMT_BSAS, method = "spearman", conf.level = 0.95)

#boxplot - Infinium Array data for all patients
ggplot(data, aes(x=factor(Patient), y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot - BSAS data for all patients
ggplot(data, aes(x=factor(Patient), y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot - RNA data for all patients
ggplot(data, aes(x=factor(Patient), y=MGMT_RNA)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log count", title = "MGMT RNA Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot - Infinium Array data for each tumor type
ggplot(agdata, aes(x=Tumor, y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#boxplot - BSAS data for each tumor type
ggplot(agdata, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#boxplot - RNA data for each tumor type
ggplot(agdata, aes(x=Tumor, y=MGMT_RNA)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log counts", title = "MGMT RNA Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#Wilcox tests comparing across primaries or across recurrences
wilcox.test(agdata_primary$MGMT_Infinium~agdata_primary$HM_status)
wilcox.test(agdata_recurrence$MGMT_Infinium~agdata_recurrence$HM_status)
wilcox.test(agdata_primary$MGMT_BSAS~agdata_primary$HM_status)
wilcox.test(agdata_primary$MGMT_RNA~agdata_primary$HM_status)
wilcox.test(agdata_recurrence$MGMT_RNA~agdata_recurrence$HM_status)

#Wilcox tests comparing between primaries and recurrences
wilcox.test(agdata$MGMT_Infinium~agdata$Tumor)
wilcox.test(agdata_HM$MGMT_Infinium~agdata_HM$Tumor)
wilcox.test(agdata_nonHM$MGMT_Infinium~agdata_nonHM$Tumor)
wilcox.test(agdata$MGMT_RNA~agdata$Tumor)
wilcox.test(agdata_HM$MGMT_RNA~agdata_HM$Tumor)
wilcox.test(agdata_nonHM$MGMT_RNA~agdata_nonHM$Tumor)

#pair data
pairdata <- merge(agdata_primary, agdata_recurrence, by = c("Patient", "HM_status"))
pairdata$Tumor.x <- NULL
pairdata$Tumor.y <- NULL
colnames(pairdata)[which(names(pairdata) == "MGMT_Infinium.x")] <- "MGMT_Infinium.Primary"
colnames(pairdata)[which(names(pairdata) == "MGMT_Infinium.y")] <- "MGMT_Infinium.Recurrence"
colnames(pairdata)[which(names(pairdata) == "MGMT_BSAS.x")] <- "MGMT_BSAS.Primary"
colnames(pairdata)[which(names(pairdata) == "MGMT_BSAS.y")] <- "MGMT_BSAS.Recurrence"
colnames(pairdata)[which(names(pairdata) == "MGMT_RNA.x")] <- "MGMT_RNA.Primary"
colnames(pairdata)[which(names(pairdata) == "MGMT_RNA.y")] <- "MGMT_RNA.Recurrence"
pairdata_infinium <- pairdata
pairdata_infinium$MGMT_BSAS.Primary <-NULL
pairdata_infinium$MGMT_BSAS.Recurrence <-NULL
pairdata_infinium$MGMT_RNA.Primary <-NULL
pairdata_infinium$MGMT_RNA.Recurrence <-NULL
pairdata_infinium <- na.omit(pairdata_infinium)
pairdata_BSAS <- pairdata
pairdata_BSAS$MGMT_Infinium.Primary <-NULL
pairdata_BSAS$MGMT_Infinium.Recurrence <-NULL
pairdata_infinium$MGMT_RNA.Primary <-NULL
pairdata_infinium$MGMT_RNA.Recurrence <-NULL
pairdata_BSAS <- na.omit(pairdata_BSAS)
pairdata_RNA <- pairdata
pairdata_RNA$MGMT_Infinium.Primary <-NULL
pairdata_RNA$MGMT_Infinium.Recurrence <-NULL
pairdata_RNA$MGMT_BSAS.Primary <-NULL
pairdata_RNA$MGMT_BSAS.Recurrence <-NULL
pairdata_RNA <- na.omit(pairdata_RNA)

#scatterplot of paired data - Infinium
ggplot(pairdata, aes(x=MGMT_Infinium.Primary, y=MGMT_Infinium.Recurrence, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Primary (Probability of Methylation)", y = "Recurrence (Probability of Methylation)", title = "MGMT Methylation Infinium Array Data For Paired Samples") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#scatterplot of paired data - RNA
ggplot(pairdata, aes(x=MGMT_RNA.Primary, y=MGMT_RNA.Recurrence, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Primary (RNA log count)", y = "Recurrence (RNA log count)", title = "MGMT RNA Data For Paired Samples") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#bar plot of paired data - Infinium
pairdata_infinium_reshape <- melt(pairdata_infinium, id.vars = c("Patient", "HM_status"))
ggplot(pairdata_infinium_reshape, aes(x=factor(Patient), y=value, fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data For Paired Samples") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_brewer(palette="Paired") +
  theme(legend.position="bottom")

#bar plot of paired data - RNA
pairdata_RNA_reshape <- melt(pairdata_RNA, id.vars = c("Patient", "HM_status"))
ggplot(pairdata_RNA_reshape, aes(x=factor(Patient), y=value, fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log count", title = "MGMT RNA Data For Paired Samples") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_brewer(palette="Paired") +
  theme(legend.position="bottom")

#box plot of paired data - Infinium
ggplot(pairdata_infinium_reshape, aes(x=variable, y=value)) +
  geom_boxplot(aes(colour=variable)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data For Paired Samples") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#box plot of paired data - RNA
ggplot(pairdata_RNA_reshape, aes(x=variable, y=value)) +
  geom_boxplot(aes(colour=variable)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log count", title = "MGMT RNA Data For Paired Samples") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#Wilcox test for paired Infinium data
pairdata_infinium_HM <- pairdata_infinium[which(pairdata_infinium$HM_status == 'HM'),]
pairdata_infinium_nonHM <- pairdata_infinium[which(pairdata_infinium$HM_status == 'nonHM'),]
wilcox.test(pairdata_infinium_HM$MGMT_Infinium.Primary, pairdata_infinium_HM$MGMT_Infinium.Recurrence, paired=TRUE, alternative = "two.sided")
wilcox.test(pairdata_infinium_nonHM$MGMT_Infinium.Primary, pairdata_infinium_nonHM$MGMT_Infinium.Recurrence, paired=TRUE, alternative = "two.sided")

#Wilcox test for paired RNA data
pairdata_RNA_HM <- pairdata_RNA[which(pairdata_RNA$HM_status == 'HM'),]
pairdata_RNA_nonHM <- pairdata_RNA[which(pairdata_RNA$HM_status == 'nonHM'),]
wilcox.test(pairdata_RNA_HM$MGMT_RNA.Primary, pairdata_RNA_HM$MGMT_RNA.Recurrence, paired=TRUE, alternative = "two.sided")
wilcox.test(pairdata_RNA_nonHM$MGMT_RNA.Primary, pairdata_RNA_nonHM$MGMT_RNA.Recurrence, paired=TRUE, alternative = "two.sided")

#patient cohort
table(agdata$HM_status)
table(agdata$Tumor)
cohort_patient <- aggregate(agdata$Patient, by =list(Tumor=agdata$Tumor, HM_status=agdata$HM_status), length)
cohort_sample <- aggregate(data$Sample, by =list(Tumor=data$Tumor, HM_status=data$HM_status), length)
cohort <- merge(cohort_patient, cohort_sample, by= c("Tumor", "HM_status"))
colnames(cohort)[which(names(cohort) == "x.x")] <- "Patients"
colnames(cohort)[which(names(cohort) == "x.y")] <- "Samples"
cohort_pairs <- aggregate(pairdata$Patient, by=list(HM_status=pairdata$HM_status), length)
names(cohort_pairs)[names(cohort_pairs) == "x"] <- "Patients with paired tumors"
