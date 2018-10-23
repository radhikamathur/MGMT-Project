#dependencies
library(ggplot2)
library(reshape2)

#data import
dirPath <- '~/Dropbox/Postdoc/MGMT Project/BSAS Results/BSAS Analysis R Project'
data <- read.csv(file.path(dirPath, "MGMT_Infinium_BSAS_RNA_FACETS_IDH_cg12981137.csv"))
molecularsubtypes <- read.csv(file.path(dirPath, "all_patients_molecular_subtype_estimation.csv"))
data <- merge(data, molecularsubtypes, by="Patient", all.x = TRUE)

#define subsets
data_primary <- data[which(data$Tumor == 'Primary'),]
data_recurrence <- data[which(data$Tumor == 'Recurrence'),]
data_HM <- data[which(data$HM_status == 'HM'),]
data_nonHM <- data[which(data$HM_status == 'nonHM'),]
data_astros <- data[which(data$mol_type == 'A'),]
data_oligos <- data[which(data$mol_type == 'O'),]
data_pure <- data[which(data$FACET_Purity >=0.695),]
data_pure_primary <- data_pure[which(data_pure$Tumor == 'Primary'),]
data_pure_recurrence <- data_pure[which(data_pure$Tumor == 'Recurrence'),]
data_pure_HM <- data_pure[which(data_pure$HM_status == 'HM'),]
data_pure_nonHM <- data_pure[which(data_pure$HM_status == 'nonHM'),]
data_pure_astros <- data_pure[which(data_pure$mol_type == 'A'),]
data_pure_oligos <- data_pure[which(data_pure$mol_type == 'O'),]

#aggregate data
agdata_Infinium <- aggregate(data$MGMT_Infinium, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean, na.rm= TRUE)
agdata_BSAS <-aggregate(data$MGMT_BSAS, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean, na.rm= TRUE)
agdata_RNA <-aggregate(data$MGMT_RNA, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean, na.rm= TRUE)
agdata <- merge(agdata_Infinium, agdata_BSAS, by = c("Patient", "Tumor", "HM_status"))
agdata <- merge(agdata, agdata_RNA, by = c("Patient", "Tumor", "HM_status"))
colnames(agdata)[which(names(agdata) == "x.x")] <- "MGMT_Infinium"
colnames(agdata)[which(names(agdata) == "x.y")] <- "MGMT_BSAS"
colnames(agdata)[which(names(agdata) == "x")] <- "MGMT_RNA"
agdata <- merge(agdata, molecularsubtypes, by="Patient", all.x = TRUE)
agdata_pure_Infinium <- aggregate(data_pure$MGMT_Infinium, by=list(Patient=data_pure$Patient, Tumor=data_pure$Tumor, HM_status=data_pure$HM_status), mean, na.rm= TRUE)
agdata_pure_BSAS <-aggregate(data_pure$MGMT_BSAS, by=list(Patient=data_pure$Patient, Tumor=data_pure$Tumor, HM_status=data_pure$HM_status), mean, na.rm= TRUE)
agdata_pure_RNA <-aggregate(data_pure$MGMT_RNA, by=list(Patient=data_pure$Patient, Tumor=data_pure$Tumor, HM_status=data_pure$HM_status), mean, na.rm= TRUE)
agdata_pure <- merge(agdata_pure_Infinium, agdata_pure_BSAS, by = c("Patient", "Tumor", "HM_status"))
agdata_pure <- merge(agdata_pure, agdata_pure_RNA, by = c("Patient", "Tumor", "HM_status"))
colnames(agdata_pure)[which(names(agdata_pure) == "x.x")] <- "MGMT_Infinium"
colnames(agdata_pure)[which(names(agdata_pure) == "x.y")] <- "MGMT_BSAS"
colnames(agdata_pure)[which(names(agdata_pure) == "x")] <- "MGMT_RNA"
agdata_pure <- merge(agdata_pure, molecularsubtypes, by="Patient", all.x = TRUE)

#aggregate data for cg12981137
agdata_Infinium_cg12981137 <- aggregate(data$MGMT_Infinium_cg12981137, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean, na.rm= TRUE)
agdata_BSAS_cg12981137 <-aggregate(data$MGMT_BSAS_cg12981137, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean, na.rm= TRUE)
agdata_cg12981137 <- merge(agdata_Infinium_cg12981137, agdata_BSAS_cg12981137, by = c("Patient", "Tumor", "HM_status"))
colnames(agdata_cg12981137)[which(names(agdata_cg12981137) == "x.x")] <- "MGMT_Infinium_cg12981137"
colnames(agdata_cg12981137)[which(names(agdata_cg12981137) == "x.y")] <- "MGMT_BSAS_cg12981137"
agdata_cg12981137 <- merge(agdata_cg12981137, molecularsubtypes, by="Patient", all.x = TRUE)
agdata_pure_Infinium_cg12981137 <- aggregate(data_pure$MGMT_Infinium_cg12981137, by=list(Patient=data_pure$Patient, Tumor=data_pure$Tumor, HM_status=data_pure$HM_status), mean, na.rm= TRUE)
agdata_pure_BSAS_cg12981137 <-aggregate(data_pure$MGMT_BSAS_cg12981137, by=list(Patient=data_pure$Patient, Tumor=data_pure$Tumor, HM_status=data_pure$HM_status), mean, na.rm= TRUE)
agdata_pure_cg12981137 <- merge(agdata_pure_Infinium_cg12981137, agdata_pure_BSAS_cg12981137, by = c("Patient", "Tumor", "HM_status"))
agdata_cg12981137_pure <- merge(agdata_pure_Infinium_cg12981137, agdata_pure_BSAS_cg12981137, by = c("Patient", "Tumor", "HM_status"))
colnames(agdata_pure_cg12981137)[which(names(agdata_pure_cg12981137) == "x.x")] <- "MGMT_Infinium_cg12981137"
colnames(agdata_pure_cg12981137)[which(names(agdata_pure_cg12981137) == "x.y")] <- "MGMT_BSAS_cg12981137"
colnames(agdata_cg12981137_pure)[which(names(agdata_cg12981137_pure) == "x.x")] <- "MGMT_Infinium_cg12981137"
colnames(agdata_cg12981137_pure)[which(names(agdata_cg12981137_pure) == "x.y")] <- "MGMT_BSAS_cg12981137"
agdata_pure_cg12981137 <- merge(agdata_pure_cg12981137, molecularsubtypes, by="Patient", all.x = TRUE)

#define subsets of aggregate data
agdata_primary <- agdata[which(agdata$Tumor == 'Primary'),]
agdata_recurrence <- agdata[which(agdata$Tumor == 'Recurrence'),]
agdata_HM <- agdata[which(agdata$HM_status == 'HM'),]
agdata_nonHM <- agdata[which(agdata$HM_status == 'nonHM'),]
agdata_astros <- agdata[which(agdata$mol_type == 'A'),]
agdata_oligos <- agdata[which(agdata$mol_type == 'O'),]
agdata_primary_astros <- agdata_astros[which(agdata_astros$Tumor == 'Primary'),]
agdata_recurrence_astros <- agdata_astros[which(agdata_astros$Tumor == 'Recurrence'),]
agdata_primary_oligos <- agdata_oligos[which(agdata_oligos$Tumor == 'Primary'),]
agdata_recurrence_oligos <- agdata_oligos[which(agdata_oligos$Tumor == 'Recurrence'),]
agdata_HM_astros <- agdata_astros[which(agdata_astros$HM_status == 'HM'),]
agdata_nonHM_astros <- agdata_astros[which(agdata_astros$HM_status == 'nonHM'),]
agdata_HM_oligos <- agdata_oligos[which(agdata_oligos$HM_status == 'HM'),]
agdata_nonHM_oligos <- agdata_oligos[which(agdata_oligos$HM_status == 'nonHM'),]
agdata_pure_primary <- agdata_pure[which(agdata_pure$Tumor == 'Primary'),]
agdata_pure_recurrence <- agdata_pure[which(agdata_pure$Tumor == 'Recurrence'),]
agdata_pure_HM <- agdata_pure[which(agdata_pure$HM_status == 'HM'),]
agdata_pure_nonHM <- agdata_pure[which(agdata_pure$HM_status == 'nonHM'),]
agdata_pure_astros <- agdata_pure[which(agdata_pure$mol_type == 'A'),]
agdata_pure_oligos <- agdata_pure[which(agdata_pure$mol_type == 'O'),]
agdata_pure_primary_astros <- agdata_pure_astros[which(agdata_pure_astros$Tumor == 'Primary'),]
agdata_pure_recurrence_astros <- agdata_pure_astros[which(agdata_pure_astros$Tumor == 'Recurrence'),]
agdata_pure_primary_oligos <- agdata_pure_oligos[which(agdata_pure_oligos$Tumor == 'Primary'),]
agdata_pure_recurrence_oligos <- agdata_pure_oligos[which(agdata_pure_oligos$Tumor == 'Recurrence'),]
agdata_pure_HM_astros <- agdata_pure_astros[which(agdata_pure_astros$HM_status == 'HM'),]
agdata_pure_nonHM_astros <- agdata_pure_astros[which(agdata_pure_astros$HM_status == 'nonHM'),]
agdata_pure_HM_oligos <- agdata_pure_oligos[which(agdata_pure_oligos$HM_status == 'HM'),]
agdata_pure_nonHM_oligos <- agdata_pure_oligos[which(agdata_pure_oligos$HM_status == 'nonHM'),]

#define subsets for aggregate data for cg12981137
agdata_cg12981137_primary <- agdata_cg12981137[which(agdata_cg12981137$Tumor == 'Primary'),]
agdata_cg12981137_recurrence <- agdata_cg12981137[which(agdata_cg12981137$Tumor == 'Recurrence'),]
agdata_cg12981137_HM <- agdata_cg12981137[which(agdata_cg12981137$HM_status == 'HM'),]
agdata_cg12981137_nonHM <- agdata_cg12981137[which(agdata_cg12981137$HM_status == 'nonHM'),]
agdata_cg12981137_pure_primary <- agdata_cg12981137_pure[which(agdata_cg12981137_pure$Tumor == 'Primary'),]
agdata_cg12981137_pure_recurrence <- agdata_cg12981137_pure[which(agdata_cg12981137_pure$Tumor == 'Recurrence'),]
agdata_cg12981137_pure_HM <- agdata_cg12981137_pure[which(agdata_cg12981137_pure$HM_status == 'HM'),]
agdata_cg12981137_pure_nonHM <- agdata_cg12981137_pure[which(agdata_cg12981137_pure$HM_status == 'nonHM'),]

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

#pair data: High purity samples only
pairdata_pure <- merge(agdata_pure_primary, agdata_pure_recurrence, by = c("Patient", "HM_status"))
pairdata_pure$Tumor.x <- NULL
pairdata_pure$Tumor.y <- NULL
colnames(pairdata_pure)[which(names(pairdata_pure) == "MGMT_Infinium.x")] <- "MGMT_Infinium.Primary"
colnames(pairdata_pure)[which(names(pairdata_pure) == "MGMT_Infinium.y")] <- "MGMT_Infinium.Recurrence"
colnames(pairdata_pure)[which(names(pairdata_pure) == "MGMT_BSAS.x")] <- "MGMT_BSAS.Primary"
colnames(pairdata_pure)[which(names(pairdata_pure) == "MGMT_BSAS.y")] <- "MGMT_BSAS.Recurrence"
colnames(pairdata_pure)[which(names(pairdata_pure) == "MGMT_RNA.x")] <- "MGMT_RNA.Primary"
colnames(pairdata_pure)[which(names(pairdata_pure) == "MGMT_RNA.y")] <- "MGMT_RNA.Recurrence"
pairdata_pure_infinium <- pairdata_pure
pairdata_pure_infinium$MGMT_BSAS.Primary <-NULL
pairdata_pure_infinium$MGMT_BSAS.Recurrence <-NULL
pairdata_pure_infinium$MGMT_RNA.Primary <-NULL
pairdata_pure_infinium$MGMT_RNA.Recurrence <-NULL
pairdata_pure_infinium <- na.omit(pairdata_pure_infinium)
pairdata_pure_BSAS <- pairdata_pure
pairdata_pure_BSAS$MGMT_Infinium.Primary <-NULL
pairdata_pure_BSAS$MGMT_Infinium.Recurrence <-NULL
pairdata_pure_infinium$MGMT_RNA.Primary <-NULL
pairdata_pure_infinium$MGMT_RNA.Recurrence <-NULL
pairdata_pure_BSAS <- na.omit(pairdata_pure_BSAS)
pairdata_pure_RNA <- pairdata_pure
pairdata_pure_RNA$MGMT_Infinium.Primary <-NULL
pairdata_pure_RNA$MGMT_Infinium.Recurrence <-NULL
pairdata_pure_RNA$MGMT_BSAS.Primary <-NULL
pairdata_pure_RNA$MGMT_BSAS.Recurrence <-NULL
pairdata_pure_RNA <- na.omit(pairdata_RNA)


#scatterplot Infinium v.BSAS 
ggplot(data, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "All individual samples: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data$MGMT_BSAS, data$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data$MGMT_BSAS, data$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot Infinium v.BSAS (by molecular subtype)
ggplot(data, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "All individual samples: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) + 
  facet_grid(mol_type ~ .)
cor.test(data_astros$MGMT_BSAS, data_astros$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_astros$MGMT_BSAS, data_astros$MGMT_Infinium, method = "spearman", conf.level = 0.95)
cor.test(data_oligos$MGMT_BSAS, data_oligos$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_oligos$MGMT_BSAS, data_oligos$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot Infinium v.BSAS: High purity samples only
ggplot(data_pure, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "High purity individual samples: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data_pure$MGMT_BSAS, data_pure$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_pure$MGMT_BSAS, data_pure$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot Infinium v.BSAS: High purity samples only (by molecular subtype)
ggplot(data_pure, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "High purity individual samples: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  facet_grid(mol_type ~ .)
cor.test(data_pure_astros$MGMT_BSAS, data_pure_astros$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_pure_astros$MGMT_BSAS, data_pure_astros$MGMT_Infinium, method = "spearman", conf.level = 0.95)
cor.test(data_pure_oligos$MGMT_BSAS, data_pure_oligos$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_pure_oligos$MGMT_BSAS, data_pure_oligos$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot Infinium v. BSAS (cg12981137)
ggplot(data, aes(x=MGMT_BSAS_cg12981137, y=MGMT_Infinium_cg12981137, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS) at cg12981137", y = "M Value at cg12981137 (Infinium Array)", title = "All individual samples: Infinium v. BSAS at cg12981137") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data$MGMT_BSAS_cg12981137, data$MGMT_Infinium_cg12981137, method = "pearson", conf.level = 0.95)
cor.test(data$MGMT_BSAS_cg12981137, data$MGMT_Infinium_cg12981137, method = "spearman", conf.level = 0.95)



#scatterplot Infinium v.BSAS (cg12981137): High purity samples only
ggplot(data_pure, aes(x=MGMT_BSAS_cg12981137, y=MGMT_Infinium_cg12981137, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS) at cg12981137", y = "M Value at cg12981137 (Infinium Array)", title = "High purity individual samples: Infinium v. BSAS at cg12981137") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data_pure$MGMT_BSAS_cg12981137, data_pure$MGMT_Infinium_cg12981137, method = "pearson", conf.level = 0.95)
cor.test(data_pure$MGMT_BSAS_cg12981137, data_pure$MGMT_Infinium_cg12981137, method = "spearman", conf.level = 0.95)

#scatterplot Infinium v. RNA
ggplot(data, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count", title = "All individual samples: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data$MGMT_RNA, data$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data$MGMT_RNA, data$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot Infinium v. RNA: High purity samples only
ggplot(data_pure, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count", title = "High purity individual samples: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data_pure$MGMT_RNA, data_pure$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_pure$MGMT_RNA, data_pure$MGMT_Infinium, method = "spearman", conf.level = 0.95)

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

#scatterplot Infinium v. RNA, facet by Tumor: High purity samples only
ggplot(data_pure, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count",  title = "High purity individual samples: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  facet_grid(. ~ Tumor, scales = "free", space = "free")
cor.test(data_pure_primary$MGMT_RNA, data_pure_primary$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_pure_primary$MGMT_RNA, data_pure_primary$MGMT_Infinium, method = "spearman", conf.level = 0.95)
cor.test(data_pure_recurrence$MGMT_RNA, data_pure_recurrence$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(data_pure_recurrence$MGMT_RNA, data_pure_recurrence$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot BSAS v. RNA
ggplot(data, aes(x=MGMT_BSAS, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Percent CpG Methylation (BSAS)", y = "RNA log count", title = "All individual samples: BSAS v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data$MGMT_RNA, data$MGMT_BSAS, method = "pearson", conf.level = 0.95)
cor.test(data$MGMT_RNA, data$MGMT_BSAS, method = "spearman", conf.level = 0.95)

#scatterplot BSAS v. RNA: High purity samples only
ggplot(data_pure, aes(x=MGMT_BSAS, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Percent CpG Methylation (BSAS)", y = "RNA log count", title = "High purity individual samples: BSAS v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(data_pure$MGMT_RNA, data_pure$MGMT_BSAS, method = "pearson", conf.level = 0.95)
cor.test(data_pure$MGMT_RNA, data_pure$MGMT_BSAS, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. BSAS
ggplot(agdata, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "All patients: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata$MGMT_BSAS, agdata$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata$MGMT_BSAS, agdata$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. BSAS: High purity samples only
ggplot(agdata_pure, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "High purity samples by patient: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata_pure$MGMT_BSAS, agdata_pure$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_pure$MGMT_BSAS, agdata_pure$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. BSAS (by molecular subtype)
ggplot(agdata, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "All patients: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  facet_grid(mol_type ~ .)
cor.test(agdata_astros$MGMT_BSAS, agdata_astros$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_astros$MGMT_BSAS, agdata_astros$MGMT_Infinium, method = "spearman", conf.level = 0.95)
cor.test(agdata_oligos$MGMT_BSAS, agdata_oligos$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_oligos$MGMT_BSAS, agdata_oligos$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. BSAS: High purity samples only (by molecular subtype)
ggplot(agdata_pure, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  geom_hline(yintercept = 0.358, linetype="dotted") +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)", title = "High purity samples by patient: Infinium v. BSAS") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  facet_grid(mol_type ~ .)
cor.test(agdata_pure_astros$MGMT_BSAS, agdata_pure_astros$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_pure_astros$MGMT_BSAS, agdata_pure_astros$MGMT_Infinium, method = "spearman", conf.level = 0.95)
cor.test(agdata_pure_oligos$MGMT_BSAS, agdata_pure_oligos$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_pure_oligos$MGMT_BSAS, agdata_pure_oligos$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. BSAS for cg12981137
ggplot(agdata_cg12981137, aes(x=MGMT_BSAS_cg12981137, y=MGMT_Infinium_cg12981137, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS) at cg12981137", y = "M Value at cg12981137 (Infinium Array)", title = "All patients: Infinium v. BSAS at cg12981137") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata_cg12981137$MGMT_BSAS, agdata_cg12981137$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_cg12981137$MGMT_BSAS, agdata_cg12981137$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. BSAS for cg12981137: High purity samples only
ggplot(agdata_pure_cg12981137, aes(x=MGMT_BSAS_cg12981137, y=MGMT_Infinium_cg12981137, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS) at cg12981137", y = "M Value at cg12981137 (Infinium Array)", title = "High purity samples by patient: Infinium v. BSAS at cg12981137") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata_pure_cg12981137$MGMT_BSAS, agdata_pure_cg12981137$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_pure_cg12981137$MGMT_BSAS, agdata_pure_cg12981137$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. RNA
ggplot(agdata, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count", title = "All Patients: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata$MGMT_RNA, agdata$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata$MGMT_RNA, agdata$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - Infinium v. RNA: High purity samples only
ggplot(agdata_pure, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count", title = "High purity samples by patient: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata_pure$MGMT_RNA, agdata_pure$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_pure$MGMT_RNA, agdata_pure$MGMT_Infinium, method = "spearman", conf.level = 0.95)

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

#scatterplot of aggregate data - Infinium v. RNA, facet by tumor: High purity samples only
ggplot(agdata_pure, aes(x=MGMT_Infinium, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Probability of Methylation (Infinium Array)", y = "RNA log count", title = "High purity samples by patient: Infinium v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  facet_grid(. ~ Tumor, scales = "free", space = "free")
cor.test(agdata_pure_primary$MGMT_RNA, agdata_pure_primary$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_pure_primary$MGMT_RNA, agdata_pure_primary$MGMT_Infinium, method = "spearman", conf.level = 0.95)
cor.test(agdata_pure_recurrence$MGMT_RNA, agdata_pure_recurrence$MGMT_Infinium, method = "pearson", conf.level = 0.95)
cor.test(agdata_pure_recurrence$MGMT_RNA, agdata_pure_recurrence$MGMT_Infinium, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - BSAS v. RNA
ggplot(agdata, aes(x=MGMT_BSAS, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Percent CpG Methylation (BSAS)", y = "RNA log count", title = "All Patients: BSAS v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata$MGMT_RNA, agdata$MGMT_BSAS, method = "pearson", conf.level = 0.95)
cor.test(agdata$MGMT_RNA, agdata$MGMT_BSAS, method = "spearman", conf.level = 0.95)

#scatterplot of aggregate data - BSAS v. RNA: High purity samples only
ggplot(agdata_pure, aes(x=MGMT_BSAS, y=MGMT_RNA, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Percent CpG Methylation (BSAS)", y = "RNA log count", title = "High purity samples by patient: BSAS v. RNA") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
cor.test(agdata_pure$MGMT_RNA, agdata_pure$MGMT_BSAS, method = "pearson", conf.level = 0.95)
cor.test(agdata_pure$MGMT_RNA, agdata_pure$MGMT_BSAS, method = "spearman", conf.level = 0.95)

#boxplot - Infinium Array data for all patients
ggplot(data, aes(x=factor(Patient), y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot - Infinium Array data for all patients: High purity samples only
ggplot(data_pure, aes(x=factor(Patient), y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot - Infinium Array data for all patients (by molecular subtype)
ggplot(data, aes(x=factor(Patient), y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(mol_type ~ HM_status, scales = "free_x", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot - Infinium Array data for all patients: High purity samples only (by molecular subtype)
ggplot(data_pure, aes(x=factor(Patient), y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(mol_type ~ HM_status, scales = "free_x", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status: Pure samples only") +
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

#boxplot - BSAS data for all patients: High purity samples only
ggplot(data_pure, aes(x=factor(Patient), y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot - BSAS data for all patients (by molecular subtype)
ggplot(data, aes(x=factor(Patient), y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(mol_type ~ HM_status, scales = "free_x", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot - BSAS data for all patients: High purity samples only (by molecular subtype)
ggplot(data_pure, aes(x=factor(Patient), y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(mol_type ~ HM_status, scales = "free_x", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status: Pure samples only") +
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

#boxplot - RNA data for all patients: High purity samples only
ggplot(data_pure, aes(x=factor(Patient), y=MGMT_RNA)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log count", title = "MGMT RNA Data By Hypermutation Status: Pure samples only") +
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
wilcox.test(agdata_primary$MGMT_Infinium~agdata_primary$HM_status)
wilcox.test(agdata_recurrence$MGMT_Infinium~agdata_recurrence$HM_status)
wilcox.test(agdata_HM$MGMT_Infinium~agdata_HM$Tumor)
wilcox.test(agdata_nonHM$MGMT_Infinium~agdata_nonHM$Tumor)

#boxplot - Infinium Array data for each tumor type: High purity samples only
ggplot(agdata_pure, aes(x=Tumor, y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_pure_primary$MGMT_Infinium~agdata_pure_primary$HM_status)
wilcox.test(agdata_pure_recurrence$MGMT_Infinium~agdata_pure_recurrence$HM_status)
wilcox.test(agdata_pure_HM$MGMT_Infinium~agdata_pure_HM$Tumor)
wilcox.test(agdata_pure_nonHM$MGMT_Infinium~agdata_pure_nonHM$Tumor)

#boxplot - Infinium Array data for each tumor type (by molecular subtype)
ggplot(agdata, aes(x=Tumor, y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(mol_type ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#boxplot - Infinium Array data for each tumor type: High purity samples only (by molecular subtype)
ggplot(agdata_pure, aes(x=Tumor, y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(mol_type ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#boxplot - Infinium Array data for each tumor type (OLIGOS)
ggplot(agdata_oligos, aes(x=Tumor, y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_primary_oligos$MGMT_Infinium~agdata_primary_oligos$HM_status)
wilcox.test(agdata_recurrence_oligos$MGMT_Infinium~agdata_recurrence_oligos$HM_status)
wilcox.test(agdata_HM_oligos$MGMT_Infinium~agdata_HM_oligos$Tumor)
wilcox.test(agdata_nonHM_oligos$MGMT_Infinium~agdata_nonHM_oligos$Tumor)

#boxplot - Infinium Array data for each tumor type: High purity samples only (OLIGOS)
ggplot(agdata_pure_oligos, aes(x=Tumor, y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_pure_primary_oligos$MGMT_Infinium~agdata_pure_primary_oligos$HM_status)
wilcox.test(agdata_pure_recurrence_oligos$MGMT_Infinium~agdata_pure_recurrence_oligos$HM_status)
wilcox.test(agdata_pure_HM_oligos$MGMT_Infinium~agdata_pure_HM_oligos$Tumor)
wilcox.test(agdata_pure_nonHM_oligos$MGMT_Infinium~agdata_pure_nonHM_oligos$Tumor)

#boxplot - Infinium Array data for each tumor type (ASTROS)
ggplot(agdata_astros, aes(x=Tumor, y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_primary_astros$MGMT_Infinium~agdata_primary_astros$HM_status)
wilcox.test(agdata_recurrence_astros$MGMT_Infinium~agdata_recurrence_astros$HM_status)
wilcox.test(agdata_HM_astros$MGMT_Infinium~agdata_HM_astros$Tumor)
wilcox.test(agdata_nonHM_astros$MGMT_Infinium~agdata_nonHM_astros$Tumor)

#boxplot - Infinium Array data for each tumor type: High purity samples only (ASTROS)
ggplot(agdata_pure_astros, aes(x=Tumor, y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_pure_primary_astros$MGMT_Infinium~agdata_pure_primary_astros$HM_status)
wilcox.test(agdata_pure_recurrence_astros$MGMT_Infinium~agdata_pure_recurrence_astros$HM_status)
wilcox.test(agdata_pure_HM_astros$MGMT_Infinium~agdata_pure_HM_astros$Tumor)
wilcox.test(agdata_pure_nonHM_astros$MGMT_Infinium~agdata_pure_nonHM_astros$Tumor)


#boxplot - Infinium Array data for each tumor type at cg12981137
ggplot(agdata_cg12981137, aes(x=Tumor, y=MGMT_Infinium_cg12981137)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "M-Value at cg12981137", title = "MGMT Methylation Infinium Array Data at cg12981137") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_cg12981137_primary$MGMT_Infinium_cg12981137~agdata_cg12981137_primary$HM_status)
wilcox.test(agdata_cg12981137_recurrence$MGMT_Infinium_cg12981137~agdata_cg12981137_recurrence$HM_status)
wilcox.test(agdata_cg12981137_HM$MGMT_Infinium_cg12981137~agdata_cg12981137_HM$Tumor)
wilcox.test(agdata_cg12981137_nonHM$MGMT_Infinium_cg12981137~agdata_cg12981137_nonHM$Tumor)

#boxplot - Infinium Array data for each tumor type at cg12981137: High purity samples only
ggplot(agdata_pure_cg12981137, aes(x=Tumor, y=MGMT_Infinium_cg12981137)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data at cg12981137: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_cg12981137_pure_primary$MGMT_Infinium_cg12981137~agdata_cg12981137_pure_primary$HM_status)
wilcox.test(agdata_cg12981137_pure_recurrence$MGMT_Infinium_cg12981137~agdata_cg12981137_pure_recurrence$HM_status)
wilcox.test(agdata_cg12981137_pure_HM$MGMT_Infinium_cg12981137~agdata_cg12981137_pure_HM$Tumor)
wilcox.test(agdata_cg12981137_pure_nonHM$MGMT_Infinium_cg12981137~agdata_cg12981137_pure_nonHM$Tumor)

#boxplot - BSAS data for each tumor type
ggplot(agdata, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_primary$MGMT_BSAS~agdata_primary$HM_status)

#boxplot - BSAS data for each tumor type: High purity samples only
ggplot(agdata_pure, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_pure_primary$MGMT_BSAS~agdata_pure_primary$HM_status)

#boxplot - BSAS data for each tumor type (by molecular subtype)
ggplot(agdata, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(mol_type ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#boxplot - BSAS data for each tumor type: High purity samples only (by molecular subtype)
ggplot(agdata_pure, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(mol_type ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#boxplot - BSAS data for each tumor type (OLIGOS)
ggplot(agdata_oligos, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_primary_oligos$MGMT_BSAS~agdata_primary_oligos$HM_status)

#boxplot - BSAS data for each tumor type: High purity samples only (OLIGOS)
ggplot(agdata_pure_oligos, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_pure_primary_oligos$MGMT_BSAS~agdata_pure_primary_oligos$HM_status)

#boxplot - BSAS data for each tumor type (ASTROS)
ggplot(agdata_astros, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_primary_astros$MGMT_BSAS~agdata_primary_astros$HM_status)

#boxplot - BSAS data for each tumor type: High purity samples only (ASTROS)
ggplot(agdata_pure_astros, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_pure_primary_astros$MGMT_BSAS~agdata_pure_primary_astros$HM_status)

#boxplot - BSAS data for each tumor type at cg12981137
ggplot(agdata_cg12981137, aes(x=Tumor, y=MGMT_BSAS_cg12981137)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data at cg12981137") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_cg12981137_primary$MGMT_BSAS_cg12981137~agdata_cg12981137_primary$HM_status)

#boxplot - BSAS data for each tumor type at cg12981137: High purity samples only
ggplot(agdata_pure_cg12981137, aes(x=Tumor, y=MGMT_BSAS_cg12981137)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Percentage CpG Methylation", title = "MGMT Methylation BSAS Data at cg12981137: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_cg12981137_pure_primary$MGMT_BSAS_cg12981137~agdata_cg12981137_pure_primary$HM_status)


#boxplot - RNA data for each tumor type
ggplot(agdata, aes(x=Tumor, y=MGMT_RNA)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log counts", title = "MGMT RNA Data By Hypermutation Status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_primary$MGMT_RNA~agdata_primary$HM_status)
wilcox.test(agdata_recurrence$MGMT_RNA~agdata_recurrence$HM_status)
wilcox.test(agdata_HM$MGMT_RNA~agdata_HM$Tumor)
wilcox.test(agdata_nonHM$MGMT_RNA~agdata_nonHM$Tumor)

#boxplot - RNA data for each tumor type: High purity samples only
ggplot(agdata_pure, aes(x=Tumor, y=MGMT_RNA)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log counts", title = "MGMT RNA Data By Hypermutation Status: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
wilcox.test(agdata_pure_primary$MGMT_RNA~agdata_pure_primary$HM_status)
wilcox.test(agdata_pure_recurrence$MGMT_RNA~agdata_pure_recurrence$HM_status)
wilcox.test(agdata_pure_HM$MGMT_RNA~agdata_pure_HM$Tumor)
wilcox.test(agdata_pure_nonHM$MGMT_RNA~agdata_pure_nonHM$Tumor)

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
pairdata_infinium_HM <- pairdata_infinium[which(pairdata_infinium$HM_status == 'HM'),]
pairdata_infinium_nonHM <- pairdata_infinium[which(pairdata_infinium$HM_status == 'nonHM'),]
wilcox.test(pairdata_infinium_HM$MGMT_Infinium.Primary, pairdata_infinium_HM$MGMT_Infinium.Recurrence, paired=TRUE, alternative = "two.sided")
wilcox.test(pairdata_infinium_nonHM$MGMT_Infinium.Primary, pairdata_infinium_nonHM$MGMT_Infinium.Recurrence, paired=TRUE, alternative = "two.sided")

#box plot of paired data - RNA
ggplot(pairdata_RNA_reshape, aes(x=variable, y=value)) +
  geom_boxplot(aes(colour=variable)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log count", title = "MGMT RNA Data For Paired Samples") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
pairdata_RNA_HM <- pairdata_RNA[which(pairdata_RNA$HM_status == 'HM'),]
pairdata_RNA_nonHM <- pairdata_RNA[which(pairdata_RNA$HM_status == 'nonHM'),]
wilcox.test(pairdata_RNA_HM$MGMT_RNA.Primary, pairdata_RNA_HM$MGMT_RNA.Recurrence, paired=TRUE, alternative = "two.sided")
wilcox.test(pairdata_RNA_nonHM$MGMT_RNA.Primary, pairdata_RNA_nonHM$MGMT_RNA.Recurrence, paired=TRUE, alternative = "two.sided")


#scatterplot of paired data - Infinium: High purity samples only
ggplot(pairdata_pure, aes(x=MGMT_Infinium.Primary, y=MGMT_Infinium.Recurrence, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Primary (Probability of Methylation)", y = "Recurrence (Probability of Methylation)", title = "MGMT Methylation Infinium Array Data For Paired Samples: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#scatterplot of paired data - RNA: High purity samples only
ggplot(pairdata_pure, aes(x=MGMT_RNA.Primary, y=MGMT_RNA.Recurrence, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Primary (RNA log count)", y = "Recurrence (RNA log count)", title = "MGMT RNA Data For Paired Samples: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

#bar plot of paired data - Infinium: High purity samples only
pairdata_pure_infinium_reshape <- melt(pairdata_pure_infinium, id.vars = c("Patient", "HM_status"))
ggplot(pairdata_pure_infinium_reshape, aes(x=factor(Patient), y=value, fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data For Paired Samples: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_brewer(palette="Paired") +
  theme(legend.position="bottom")

#bar plot of paired data - RNA: High purity samples only
pairdata_pure_RNA_reshape <- melt(pairdata_RNA, id.vars = c("Patient", "HM_status"))
ggplot(pairdata_pure_RNA_reshape, aes(x=factor(Patient), y=value, fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_grid(. ~ HM_status, scales = "free", space = "free", labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log count", title = "MGMT RNA Data For Paired Samples: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_brewer(palette="Paired") +
  theme(legend.position="bottom")

#box plot of paired data - Infinium: High purity samples only
ggplot(pairdata_pure_infinium_reshape, aes(x=variable, y=value)) +
  geom_boxplot(aes(colour=variable)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  geom_hline(yintercept = 0.358, col = "red", linetype="dotted") +
  labs(x= "Patient", y= "Probability of Methylation", title = "MGMT Methylation Infinium Array Data For Paired Samples: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
pairdata_pure_infinium_HM <- pairdata_pure_infinium[which(pairdata_pure_infinium$HM_status == 'HM'),]
pairdata_pure_infinium_nonHM <- pairdata_pure_infinium[which(pairdata_pure_infinium$HM_status == 'nonHM'),]
wilcox.test(pairdata_pure_infinium_HM$MGMT_Infinium.Primary, pairdata_pure_infinium_HM$MGMT_Infinium.Recurrence, paired=TRUE, alternative = "two.sided")
wilcox.test(pairdata_pure_infinium_nonHM$MGMT_Infinium.Primary, pairdata_pure_infinium_nonHM$MGMT_Infinium.Recurrence, paired=TRUE, alternative = "two.sided")

#box plot of paired data - RNA: High purity samples only
ggplot(pairdata_pure_RNA_reshape, aes(x=variable, y=value)) +
  geom_boxplot(aes(colour=variable)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated", nonHM = "Non-Hypermutated"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x= "Patient", y= "RNA log count", title = "MGMT RNA Data For Paired Samples: Pure samples only") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
pairdata_pure_RNA_HM <- pairdata_pure_RNA[which(pairdata_pure_RNA$HM_status == 'HM'),]
pairdata_pure_RNA_nonHM <- pairdata_pure_RNA[which(pairdata_pure_RNA$HM_status == 'nonHM'),]
wilcox.test(pairdata_pure_RNA_HM$MGMT_RNA.Primary, pairdata_pure_RNA_HM$MGMT_RNA.Recurrence, paired=TRUE, alternative = "two.sided")
wilcox.test(pairdata_pure_RNA_nonHM$MGMT_RNA.Primary, pairdata_pure_RNA_nonHM$MGMT_RNA.Recurrence, paired=TRUE, alternative = "two.sided")

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

#pure patient cohort
table(agdata_pure$HM_status)
table(agdata_pure$Tumor)
pure_cohort_patient <- aggregate(agdata_pure$Patient, by =list(Tumor=agdata_pure$Tumor, HM_status=agdata_pure$HM_status), length)
pure_cohort_sample <- aggregate(data_pure$Sample, by =list(Tumor=data_pure$Tumor, HM_status=data_pure$HM_status), length)
pure_cohort <- merge(pure_cohort_patient, pure_cohort_sample, by= c("Tumor", "HM_status"))
colnames(pure_cohort)[which(names(pure_cohort) == "x.x")] <- "Patients"
colnames(pure_cohort)[which(names(pure_cohort) == "x.y")] <- "Samples"
pure_cohort_pairs <- aggregate(pairdata_pure$Patient, by=list(HM_status=pairdata_pure$HM_status), length)
names(pure_cohort_pairs)[names(pure_cohort_pairs) == "x"] <- "Patients with paired tumors"

## For Joe's donor meeting
#boxplot - Infinium Array data for each tumor type: High purity samples only - Primary only
ggplot(agdata_pure_primary, aes(x=Tumor, y=MGMT_Infinium)) +
  geom_boxplot(aes(colour=Tumor)) +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated at recurrence", nonHM = "Non-Hypermutated at recurrence"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(y= "MGMT promoter methylation in primary tumor (Probability)", title = "Method 1: Infinium Array") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(legend.position="none")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

#boxplot - BSAS data for each tumor type: High purity samples only - Primary only
ggplot(agdata_pure_primary, aes(x=Tumor, y=MGMT_BSAS)) +
  geom_boxplot(aes(colour=Tumor)) +
  scale_color_brewer(palette="Set2") +
  facet_grid(. ~ HM_status, labeller = labeller(HM_status = c(HM = "Hypermutated at recurrence", nonHM = "Non-Hypermutated at recurrence"))) +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(y= "MGMT promoter methylation in primary tumor (Percent)", title = "Method 2: Bisulfite Amplicon Sequencing") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  theme(legend.position="none")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

# FACETS v. IDH
ggplot(data, aes(x=FACET_Purity, y=IDH_Purity, color = mol_type)) +
  geom_point() 
cor.test(data$FACET_Purity, data$IDH_Purity, method = "pearson", conf.level = 0.95)
cor.test(data$FACET_Purity, data$IDH_Purity, method = "spearman", conf.level = 0.95)

# FACETS v. IDH
ggplot(data, aes(x=FACET_Purity, y=IDH_Purity)) +
  geom_point() +
  facet_grid(.~mol_type)
cor.test(data_astros$FACET_Purity, data_astros$IDH_Purity, method = "pearson", conf.level = 0.95)
cor.test(data_astros$FACET_Purity, data_astros$IDH_Purity, method = "spearman", conf.level = 0.95)
cor.test(data_oligos$FACET_Purity, data_oligos$IDH_Purity, method = "pearson", conf.level = 0.95)
cor.test(data_oligos$FACET_Purity, data_oligos$IDH_Purity, method = "spearman", conf.level = 0.95)
