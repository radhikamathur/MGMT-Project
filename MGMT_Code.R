#dependencies

library(pacman)
p_load(ggplot2, reshape2, EnvStats, ggpubr, ggsignif, psych, pheatmap, RColorBrewer, dplyr, plyr, table1, cowplot, magick, writexl, tidyverse, data.table, caret, plotROC, lmtest, ROCR, aod, survey, partDSA, rpart)


#functions
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#data import
dirPath <- '~/Dropbox/Postdoc/Papers/MGMT Paper/Code'
data_identifiers <- read.csv(file.path(dirPath, "Sample_info.csv"))
allCGs <- read.csv(file.path(dirPath, "BSAS_data.csv"), check.names=FALSE)
patientannotations <- read.csv(file.path(dirPath, "Patient_Annotations.csv"))
Infinium <- read.csv(file.path(dirPath, "Infinium.csv"))
RNA <- read.csv(file.path(dirPath, "RNA.csv"))

#data organization
BySample <- merge(data_identifiers, allCGs, by = "Sample")
BySample$PercentMethylation <- rowMeans(BySample[,7:33], na.rm = TRUE)
ByPatient <- aggregate(BySample$PercentMethylation, by = list(Patient = BySample$Patient, Tumor=BySample$Tumor, mol_type=BySample$mol_type, HM_status=BySample$HM_status), mean, na.rm= TRUE)
colnames(ByPatient)[which(names(ByPatient) == "x")] <- "Percent_Methylation"
ByPatient$TumorType <- paste0(ByPatient$Tumor, "_", ByPatient$HM_status)
AnnotatedPatients <- merge(ByPatient, patientannotations, by="Patient", all.x=TRUE)

##########################################
# TABLE 1 (Description of patient cohort)
##########################################

Dem <- AnnotatedPatients[which(AnnotatedPatients$Tumor == "Initial"),]
Dem <- select(Dem, mol_type, Age.at.diagnosis, Sex, Grade.at.recurrence, TMZ.cycles.prior.to.recurrence,HM_status)
Dem$Grade.at.recurrence <- as.factor(Dem$Grade.at.recurrence)
Dem$mol_type <- revalue(Dem$mol_type, c("A" = "Astrocytoma"))
Dem$mol_type <- revalue(Dem$mol_type, c("O" = "Oligodendroglioma"))

table1::label(Dem$mol_type) <- "Molecular Subtype"
table1::label(Dem$Age.at.diagnosis) <- "Age at diagnosis"
table1::label(Dem$Grade.at.recurrence) <- "Grade at recurrence"
table1::label(Dem$TMZ.cycles.prior.to.recurrence) <- "TMZ cycles"
table1::label(Dem$HM_status) <- "HM status at recurrence"

table1::table1(~Age.at.diagnosis + Sex + TMZ.cycles.prior.to.recurrence + Grade.at.recurrence + HM_status | mol_type, data = Dem)

##########################################
# FIGURE 1
##########################################

# Plot by grade at recurrence (Figure 1a)
Fig1a <- ggplot(data=Dem, aes(x=Grade.at.recurrence))+
  geom_bar(color="black", stat="count", aes(fill=HM_status))+
  geom_text(stat="count", aes(label = ..count..), vjust=-0.5)+
  labs(x="Grade at recurrence", y="Patients", fill = "Hypermutation (HM) status")+
  facet_grid(.~mol_type, scales="free", space="free")+
  theme(panel.grid.minor = element_blank())+
  scale_fill_manual(labels = c("HM", "Mixed", "nonHM"), values = c("#F8766D", "#E7B800","#00BFC4"))+ 
  theme_bw(base_size=11)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(size=11), strip.background = element_rect(colour="white", fill="white"), axis.text=element_text(size=11))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24))+
  theme(legend.position="bottom", legend.text=element_text(size=11))

# Remove mixed samples
DemNoMixed <- Dem[!(Dem$HM_status == "Mixed"),]

# Plot by age (Figure 1b)
Fig1b <- ggplot(data=DemNoMixed, aes(x=HM_status, y=Age.at.diagnosis))+
  geom_boxplot(aes(fill=HM_status))+
  ylim(15,55)+
  scale_color_brewer(palette="Dark2")+
  labs(x=NULL)+
  stat_compare_means()+
  facet_grid(.~mol_type, scales="free", space="free")+
  labs(x="Hypermutation (HM) status", y="Age at diagnosis (years)")+ 
  theme_bw(base_size=11)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(size=11), strip.background = element_rect(colour="white", fill="white"), axis.text=element_text(size=11))+
  theme(legend.position="none")

# Plot by TMZ cycles (Figure 1c)
Fig1c <- ggplot(data=DemNoMixed, aes(x=HM_status, y=TMZ.cycles.prior.to.recurrence))+
  geom_boxplot(aes(fill=HM_status))+
  ylim(0,55)+
  scale_color_brewer(palette="Dark2")+
  labs(x=NULL)+
  stat_compare_means()+
  facet_grid(.~mol_type, scales="free", space="free")+
  labs(x="Hypermutation (HM) status", y="TMZ cycles received")+ 
  theme_bw(base_size=11)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(size=11), strip.background = element_rect(colour="white", fill="white"), axis.text=element_text(size=11))+
  theme(legend.position="none")

Figure1 <- plot_grid(Fig1a, Fig1b, Fig1c, labels = c("a", "b", "c"), nrow = 1,  rel_widths = c(1.2,0.9,0.9))
save_plot('~/Dropbox/Postdoc/Papers/MGMT Paper/Figures/Figure1.pdf', Figure1, base_width = 12, base_height =5.5, dpi = 300)

#############################################################################################################################
# FIGURE 2
##############################################################################################################################

# Comparison with Infinium (Fig 2a)
InfiniumPlot <- merge(BySample, Infinium, by="Sample")
Fig2a <- ggplot(InfiniumPlot, aes(x=`174`, y=MGMT_Infinium_cg12981137, color=Tumor)) + 
  geom_point() + 
  labs(x= "BSAS (% methylation)", y = "Infinium Array (M value)") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")+
  facet_grid(Tumor~.)+ 
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", label.x = 3, label.y = 4, size =3)+
  theme_bw(base_size=8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(size=8), strip.background = element_rect(colour="white", fill="white"), axis.text=element_text(size=8))+
  theme(legend.position="none")

#Correlation matrix
RNAdata <- merge(BySample, RNA, by="Sample", all.x=TRUE)
RNACorr <- corr.test(RNAdata[,7:34], RNAdata[,35], use="pairwise", method="pearson")
RNACorrMatrix <- as.data.frame(RNACorr$r)
RNACorrMatrix$CpG <- rownames(RNACorrMatrix)
RNACorrMatrix$p <- RNACorr$p
colnames(RNACorrMatrix)[which(names(RNACorrMatrix) == "V1")] <- "PearsonCorr"
RNACorrCpG <- RNACorrMatrix[-c(28),]

#Plot Pearson Corr (Fig 2b - upper)
Fig2bTop <- ggplot(data=RNACorrCpG, aes(x=as.numeric(CpG), y=PearsonCorr)) +
  geom_line()+
  geom_point()+
  labs(x= NULL, y = "Pearson's r")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")+
  theme_bw(base_size=8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position="top")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7), axis.text.y = element_text(size = 8))+
  scale_x_continuous(limits= c(87,331), breaks = c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308), labels=c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308)) 

#Plot p-values (Fig 2b - lower)
Fig2bBottom <- ggplot(data=RNACorrCpG, aes(x=as.numeric(CpG), y=p)) +
  geom_line()+
  geom_point()+
  labs(x= NULL, y = "Adjusted P values") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")+
  theme_bw(base_size=8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position="top")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7), axis.text.y = element_text(size = 8))+
  geom_hline(yintercept = 0.05, color = "red", linetype =2)+
  annotate("text", x = 325, y = 0.15, label = "P = 0.05", color = "red", size = 3)+
  scale_x_continuous(limits= c(87,331), breaks = c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308), labels=c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308)) 

Fig2bAnnotation <- ggdraw()+
  draw_image(file.path(dirPath,"SequenceAnnotation.png"))

Fig2b <- plot_grid(Fig2bTop, Fig2bBottom, Fig2bAnnotation, nrow = 3, rel_heights = c(0.6,0.6,0.4))
Figure2 <- plot_grid(Fig2a, Fig2b, labels = c("a", "b"), label_size = 10, nrow = 1, rel_widths = c(1, 2.5))
save_plot('~/Dropbox/Postdoc/Papers/MGMT Paper/Figures/Figure2.pdf', Figure2, base_width = 8.5, base_height =4, dpi = 300)

##############################################################################################################################
# FIGURE 3
##############################################################################################################################

# Heatmaps of recurrences by sample (FIGURE 3a)
BySample_Recurrences <- BySample[which(BySample$Tumor == "Recurrence"),]
BySample_Recurrences$HM_by_sample <- relevel(BySample_Recurrences$HM_by_sample, ref = 'Yes')

annotation_row = BySample_Recurrences[c("mol_type", "HM_by_sample")]
colnames(annotation_row)[which(names(annotation_row) == "mol_type")] <- "Subtype"
colnames(annotation_row)[which(names(annotation_row) == "HM_by_sample")] <- "Hypermutation"
labels_row = as.factor(BySample_Recurrences$Patient)
pheatmap(BySample_Recurrences[,7:33], cluster_cols=FALSE, annotation_row = annotation_row, labels_row = labels_row, annotation_names_row=F)
Fig3aLeft  <- ggdraw()+
  draw_image(file.path(dirPath,"Heatmap_all27.png"))

BySample_Recurrences_concordant <- BySample_Recurrences[c("106", "113", "142")]
pheatmap(BySample_Recurrences_concordant, cluster_cols=FALSE, annotation_row = annotation_row, labels_row = labels_row, annotation_names_row=F)
Fig3aRight  <- ggdraw()+
  draw_image(file.path(dirPath,"Heatmap_3.png"))

#Intratumoral analysis (FIGURE 3b)
Mixed <- BySample[which(BySample$HM_status == "Mixed"),]
ITHplot <- Mixed[-8,] 
colnames(ITHplot)[which(names(ITHplot) == "HM_by_sample")] <- "HM"
patientlabels <- c(`1` = "Patient 1", `33` = "Patient 33", `94` = "Patient 94")
Fig3b <- ggplot(ITHplot, aes(x=Tumor, y=PercentMethylation, color=Tumor)) +
  geom_jitter(aes(shape=HM), size=3)+ 
  scale_shape_manual(values=c(1,13))+
  facet_wrap(Patient~., labeller=labeller(Patient = patientlabels))+
  labs(x= NULL, y = "Average methylation (%)", shape = "Hypermutation")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")+
  theme_bw(base_size = 11)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="white", fill="white"))+
  theme(legend.position="bottom", legend.title=element_text(size=10))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  guides(color =guide_legend(nrow=3, title.position = "top"), shape =guide_legend(nrow=3, title.position = "top"))

Figure3 <- plot_grid(Fig3aLeft, Fig3aRight, Fig3b, labels = c("a", "", "b"), nrow = 1, rel_widths = c(1,0.55,0.55))
save_plot('~/Dropbox/Postdoc/Papers/MGMT Paper/Figures/Figure3.pdf', Figure3, base_width = 11, base_height =5, dpi = 300)

##############################################################################################################################
# FIGURE 4
##############################################################################################################################

# Methylation landscape by HM status and tumor type (FIGURE 4a)
allCGs <- merge(data_identifiers, allCGs, by="Sample")
allCGs_reshape <- melt(allCGs, id.vars = c("Sample", "Patient", "HM_status", "HM_by_sample", "Tumor", "mol_type"))
colnames(allCGs_reshape)[which(names(allCGs_reshape) == "variable")] <- "CpG"
colnames(allCGs_reshape)[which(names(allCGs_reshape) == "value")] <- "Percent_Methylation"

agdata_allCGs <- aggregate(allCGs_reshape$Percent_Methylation, by=list(Patient=allCGs_reshape$Patient, Tumor=allCGs_reshape$Tumor, HM_status=allCGs_reshape$HM_status, mol_type = allCGs_reshape$mol_type, CpG = allCGs_reshape$CpG), mean, na.rm= TRUE)
colnames(agdata_allCGs)[which(names(agdata_allCGs) == "x")] <- "Percent_Methylation"

agdata2_allCGs <- aggregate(agdata_allCGs$Percent_Methylation, by=list(Tumor=agdata_allCGs$Tumor, HM_status=agdata_allCGs$HM_status, CpG = agdata_allCGs$CpG), mean, na.rm= TRUE)
colnames(agdata2_allCGs)[which(names(agdata2_allCGs) == "x")] <- "Percent_Methylation"
agdata2_allCGs <- agdata2_allCGs[!(agdata2_allCGs$HM_status == "Mixed"),]
agdata2_allCGs$TumorType = paste(agdata2_allCGs$Tumor, agdata2_allCGs$HM_status)

Fig4a <- ggplot(data=agdata2_allCGs, aes(x=as.numeric.factor(CpG), y=Percent_Methylation, group=TumorType, color=Tumor)) +
  geom_line(aes(linetype=HM_status))+
  geom_point()+
  labs(x= "CpG sites (+bp from TSS)", y = "Methylation (%)", linetype="Hypermutation (HM) status") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  scale_color_brewer(palette="Set1")+
  theme_bw(base_size=14)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position="bottom")+ 
  #theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_continuous(limits= c(85,325), breaks = c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308), labels=c(93,95,106,113,118,121,125,135,137,142,147,153,174,179,185,195,208,213,225,241,252,255,258,270,291,295,308)) 

Figure <- ByPatient
Figure$mol_type <- revalue(Figure$mol_type, c("A" = "Astrocytoma"))
Figure$mol_type <- revalue(Figure$mol_type, c("O" = "Oligodendroglioma"))
colnames(Figure)[which(names(Figure) == "mol_type")] <- "Subtype"

#Boxplots by HM status  (FIGURE 4b)
Fig4b <- ggplot(data=Figure, aes(x=TumorType, y=Percent_Methylation, group=TumorType))+
  geom_boxplot(aes(color=Tumor))+
  geom_jitter(aes(color=Tumor, shape=Subtype))+
  ylim(0, 100)+
  #stat_mean_sd_text(digits=2,y.pos=100, size = 3.5)+
  labs(x=NULL, y = "Average methylation (%)")+
  scale_color_brewer(palette="Set1")+
  geom_signif(test=wilcox.test, textsize=4, y_position=c(80,80,88), comparisons=list(c("Initial_HM", "Initial_nonHM"),c("Recurrence_HM", "Recurrence_nonHM"),c("Initial_nonHM", "Recurrence_nonHM")))+
  theme_bw(base_size=14)+
  scale_shape_manual(values=c(17, 19))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_discrete(labels=c("Initial_HM" = "HM", "Initial_Mixed" = "Mixed",
                            "Initial_nonHM" = "nonHM", "Recurrence_HM" = "HM", "Recurrence_Mixed" = "Mixed",
                            "Recurrence_nonHM" = "nonHM"))

Fig4c <- ggdraw()+
  draw_image(file.path(dirPath,"PartDSA.png"))

bottomrow <- plot_grid(Fig4b, Fig4c, labels = c("b", "c"), nrow = 1, rel_widths = c(2.5,1))
Figure4 <- plot_grid(Fig4a, bottomrow, labels = c("a", ""), nrow = 2, rel_widths = c(1,1))
save_plot('~/Dropbox/Postdoc/Papers/MGMT Paper/Figures/Figure4.pdf', Figure4, base_width = 12, base_height =7.5, dpi = 300)

##############################################################################################################################
# SUPPLEMENTARY FIGURE 2
##############################################################################################################################

# Samples included in methylation analysis
FigureS2 <- ggplot(data=BySample, aes(x=as.factor(Patient)))+
  geom_bar(stat="count", aes(fill=Tumor))+
  scale_fill_brewer(palette="Set1")+
  geom_text(stat="count", aes(label = ..count..), vjust=-0.5)+
  labs(x="Patient ID", y="Samples For Methylation Analysis")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11))+
  facet_grid(.~mol_type, scales="free", space="free", labeller = labeller(mol_type = c(A = "Astrocytoma", O = "Oligodendroglioma")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="white", fill="white"))
save_plot('~/Dropbox/Postdoc/Papers/MGMT Paper/Figures/FigureS2.pdf', FigureS2, base_width = 11, base_height =5, dpi = 300)

##############################################################################################################################
# SUPPLEMENTARY FIGURE 3
##############################################################################################################################

colnames(RNAdata)[which(names(RNAdata) == "mol_type")] <- "Subtype"

# RNA correlation plot (average percent methylation across CpG) 
FigS3Legend <- ggplot(RNAdata, aes(x=PercentMethylation, y=MGMT_RNA2, color=Tumor, shape=Subtype)) +
  scale_shape_manual(values=c(17, 19))+
  geom_point() + 
  #geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "RNA log count", title = "Across CpGs (R = -0.43)") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  #geom_smooth(method=lm)+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold"))
cor.test(RNAdata$PercentMethylation, RNAdata$MGMT_RNA2, method = "pearson", conf.level = 0.95)

FigS3a <- ggplot(RNAdata, aes(x=PercentMethylation, y=MGMT_RNA2, color=Tumor, shape=Subtype)) +
  scale_shape_manual(values=c(17, 19))+
  geom_point() + 
  #geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "RNA log count", title = "Across CpGs (R = -0.43)") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  #geom_smooth(method=lm)+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold"))+ 
  theme(legend.position = "none")
cor.test(RNAdata$PercentMethylation, RNAdata$MGMT_RNA2, method = "pearson", conf.level = 0.95)

FigS3b <- ggplot(RNAdata, aes(x=RNAdata[,9], y=MGMT_RNA2, color=Tumor, shape=Subtype)) +
  scale_shape_manual(values=c(17, 19))+
  geom_point() + 
  #geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "RNA log count", title = "CpG site 106 (R = -0.48)") +
  theme(plot.title = element_text(face="bold"))+
  #geom_smooth(method=lm)+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold"))+ 
  theme(legend.position = "none")
cor.test(RNAdata[,9], RNAdata$MGMT_RNA2, method = "pearson", conf.level = 0.95)

FigS3c <- ggplot(RNAdata, aes(x=RNAdata[,10], y=MGMT_RNA2, color=Tumor, shape=Subtype)) +
  scale_shape_manual(values=c(17, 19))+ 
  geom_point() + 
  #geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "RNA log count", title = "CpG site 113 (R = -0.50)") +
  theme(plot.title = element_text(face="bold"))+
  #geom_smooth(method=lm)+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold"))+ 
  theme(legend.position = "none")
cor.test(RNAdata[,10], RNAdata$MGMT_RNA2, method = "pearson", conf.level = 0.95)

FigS3d <-ggplot(RNAdata, aes(x=RNAdata[,16], y=MGMT_RNA2, color=Tumor, shape=Subtype)) +
  scale_shape_manual(values=c(17, 19))+ 
  geom_point() + 
  #geom_text(aes(label=Patient), nudge_x=2.5) +
  labs(x= "Percent Methylation (BSAS)", y = "RNA log count", title = "CpG site 142 (R = -0.52)") +
  theme(plot.title = element_text(face="bold"))+
  #geom_smooth(method=lm)+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face="bold"))+ 
  theme(legend.position = "none")
cor.test(RNAdata[,16], RNAdata$MGMT_RNA2, method = "pearson", conf.level = 0.95)

legendS3 <- get_legend(FigS3Legend)
FigureS3 <- plot_grid(FigS3a, FigS3b, legendS3, FigS3c, FigS3d, nrow = 2, rel_widths = c(1,1,0.25,1,1))
save_plot('~/Dropbox/Postdoc/Papers/MGMT Paper/Figures/FigureS3.pdf', FigureS3, base_width = 10, base_height =8, dpi = 300)

##############################################################################################################################
# SUPPLEMENTARY FIGURE 3
##############################################################################################################################

#Boxplots by Tumor (FIGURE S4a)
FigS4a <- ggplot(data=Figure, aes(x=Tumor, y=Percent_Methylation))+
  geom_boxplot()+
  geom_jitter(aes(color=Tumor, shape=Subtype))+
  ylim(0, 100)+
  #stat_mean_sd_text(digits=2,y.pos=0, size = 3.5)+
  labs(x=NULL, y = "Average methylation (%)")+
  stat_compare_means()+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
  scale_shape_manual(values=c(17, 19))

#Boxplots by molecular type (FIGURE S4b)
FigS4b <- ggplot(data=Figure, aes(x=Subtype, y=Percent_Methylation))+
  geom_boxplot()+
  geom_jitter(aes(color=Tumor, shape=Subtype))+
  ylim(0, 100)+
  #stat_mean_sd_text(digits=2,y.pos=0, size = 3.5)+
  labs(x=NULL, y = "Average methylation (%)")+
  facet_grid(.~Tumor)+
  stat_compare_means()+
  scale_color_brewer(palette="Set1")+
  theme_bw(base_size = 11)+
  scale_shape_manual(values=c(17, 19))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(size=11), strip.background = element_rect(colour="white", fill="white"))

FigureS4 <- plot_grid(FigS4a, FigS4b, labels = c("a", "b"), nrow = 1, rel_widths = c(1,2))
save_plot('~/Dropbox/Postdoc/Papers/MGMT Paper/Figures/FigureS4.pdf', FigureS4, base_width = 10, base_height =4, dpi = 300)


##############################################################################################################################
# LOGISTIC REGRESSION
##############################################################################################################################

ForRegression <- dcast(agdata_allCGs, Patient + Tumor + HM_status + mol_type ~ CpG)
ForRegression$PercentMethylation <- rowMeans(ForRegression[,5:31], na.rm = TRUE)
ForRegression <- ForRegression[which(ForRegression$Tumor == "Initial"),]
ForRegression <- merge(ForRegression, patientannotations, by= "Patient", all.x=TRUE)
ForRegression$HM_status[ForRegression$HM_status == "Mixed"] <- "HM"
MGMT <- ForRegression %>% mutate(HM_status = relevel(HM_status, ref = 'nonHM'))
MGMT$Age.at.diagnosis <- as.numeric(as.character(MGMT$Age.at.diagnosis))

#Individual CpG sites
cols = names(MGMT[5:31])
model_summary = NULL
for(i in seq_along(cols)){
  model_fit <- glm(HM_status ~ MGMT[[cols[i]]],data = MGMT, family = "binomial")
  OR_est <- exp(cbind(OR=coef(model_fit),confint(model_fit)))
  p_value <- as.numeric(coef(summary(model_fit))[,4][2])
  out <- c(OR_est[2,], p_value)
  model_summary = rbind(model_summary, out)
}
rownames(model_summary) <- cols
colnames(model_summary)[4] <- 'p.value'
model_summary

#PartDSA
varNames = c("mol_type","PercentMethylation","Age.at.diagnosis","Sex","Extent.of.resection","TERTp.status","Initial.tumor.grade","TMZ.cycles.prior.to.recurrence")
DSA.control(vfold=10, minsplit=10, minbuck=round(10/3), cut.off.growth=10, MPD=0.1,
            missing="default", loss.function="default")
model<-partDSA(x=MGMT[,varNames],y=MGMT$HM_status,
               control=DSA.control(missing="impute.at.split") )
#showDSA(model)

#Univariate model
mylogit <- glm(HM_status~PercentMethylation,data = MGMT, family = "binomial")
summary(mylogit)
mylogit_OR = exp(cbind(OR=coef(mylogit),confint(mylogit)))
mylogit_OR
## That is, we are 95% confident that the odds of hypermutation are 1.01 to 1.13 times greater for each additional increase in percentMethylation.

#Multivariate models
mylogit <- glm(HM_status~PercentMethylation+mol_type,data = MGMT, family = "binomial")
summary(mylogit)
mylogit_OR = exp(cbind(OR=coef(mylogit),confint(mylogit)))
mylogit_OR

mylogit <- glm(HM_status~PercentMethylation+Age.at.diagnosis,data = MGMT, family = "binomial")
summary(mylogit)
mylogit_OR = exp(cbind(OR=coef(mylogit),confint(mylogit)))
mylogit_OR

mylogit <- glm(HM_status~PercentMethylation+Age.at.diagnosis+mol_type,data = MGMT, family = "binomial")
summary(mylogit)
mylogit_OR = exp(cbind(OR=coef(mylogit),confint(mylogit)))
mylogit_OR

##############################################################################################################################
# SUPPLEMENTARY TABLE
##############################################################################################################################

PatientList <- as.data.frame(ForRegression$Patient)
colnames(PatientList)[which(names(PatientList) == "ForRegression$Patient")] <- "Patient"

Johnson2014 <- list("1", "2", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "15", "16", "17", "18", "21", "24", "25", "26", "27", "28", "29")
vanThuijl2015 <- list(1,5,10,18,21,4,6,11,17,2,7,12,13,14,16,22)
Mazor2015_Exome <- list(1,2,3,4,7,8,10,11,12,13,14,16,17,18,22,36,37,38,49,68,90)
Mazor2015_RNASeq <- list(1,3,4,7,10,12,13,14,16,17,22,37,38)
Mazor2015_Infinium <- list(1,2,3,4,7,8,10,11,12,13,14,16,17,18,22,36,38,49,68,90)
Mazor2017_Exome <- list(14,169,68,17,21,27)
Mazor2017_RNASeq <- list(14,169,17)
Mazor2017_Infinium <- list(14,169,17,68,21)

PatientList$BSAS <- "New"
PatientList$ExomeSeq <- "New"
PatientList$ExomeSeq[PatientList$Patient %in% Johnson2014 | PatientList$Patient %in% vanThuijl2015 | PatientList$Patient %in% Mazor2015_Exome | PatientList$Patient %in% Mazor2017_Exome] <- "Published"

RNAdata <- na.omit(RNAdata)
PatientList$RNASeq[PatientList$Patient %in% RNAdata$Patient] <- "New"
PatientList$RNASeq[PatientList$Patient %in% Mazor2015_RNASeq | PatientList$Patient %in% Mazor2017_RNASeq] <- "Published"

InfiniumPlot <- na.omit(InfiniumPlot)
PatientList$Infinium[PatientList$Patient %in% InfiniumPlot$Patient] <- "New"
PatientList$Infinium[PatientList$Patient %in% Mazor2015_Infinium | PatientList$Patient %in% Mazor2017_Infinium] <- "Published"

PatientList$Johnson2014[PatientList$Patient %in% Johnson2014] <- "Exome-Seq"
PatientList$vanThuijl2015[PatientList$Patient %in% vanThuijl2015] <- "Exome-Seq"
PatientList$Mazor2015[PatientList$Patient %in% Mazor2015_Exome] <- "Exome-Seq"
PatientList$Mazor2015[PatientList$Patient %in% Mazor2015_RNASeq] <- "RNA-Seq"
PatientList$Mazor2015[PatientList$Patient %in% Mazor2015_Infinium] <- "Infinium"
PatientList$Mazor2015[PatientList$Patient %in% Mazor2015_Exome & PatientList$Patient %in% Mazor2015_RNASeq] <- "Exome-Seq, RNA-Seq"
PatientList$Mazor2015[PatientList$Patient %in% Mazor2015_Exome & PatientList$Patient %in% Mazor2015_Infinium] <- "Exome-Seq, Infinium"
PatientList$Mazor2015[PatientList$Patient %in% Mazor2015_RNASeq & PatientList$Patient %in% Mazor2015_Infinium] <- "RNA-Seq, Infinium"
PatientList$Mazor2015[PatientList$Patient %in% Mazor2015_Exome & PatientList$Patient %in% Mazor2015_RNASeq & PatientList$Patient %in% Mazor2015_Infinium] <- "Exome-Seq, RNA-Seq, Infinium"
PatientList$Mazor2017[PatientList$Patient %in% Mazor2017_Exome] <- "Exome-Seq"
PatientList$Mazor2017[PatientList$Patient %in% Mazor2017_RNASeq] <- "RNA-Seq"
PatientList$Mazor2017[PatientList$Patient %in% Mazor2017_Infinium] <- "Infinium"
PatientList$Mazor2017[PatientList$Patient %in% Mazor2017_Exome & PatientList$Patient %in% Mazor2017_RNASeq] <- "Exome-Seq, RNA-Seq"
PatientList$Mazor2017[PatientList$Patient %in% Mazor2017_Exome & PatientList$Patient %in% Mazor2017_Infinium] <- "Exome-Seq, Infinium"
PatientList$Mazor2017[PatientList$Patient %in% Mazor2017_RNASeq & PatientList$Patient %in% Mazor2017_Infinium] <- "RNA-Seq, Infinium"
PatientList$Mazor2017[PatientList$Patient %in% Mazor2017_Exome & PatientList$Patient %in% Mazor2017_RNASeq & PatientList$Patient %in% Mazor2017_Infinium] <- "Exome-Seq, RNA-Seq, Infinium"

#SuppTable <- list(Patient_List = PatientList, Sample_info = data_identifiers, Patient_annotations = patientannotations, BSAS_data = allCGs, Infinium_data = Infinium, RNA_data = RNA)
#write_xlsx(SuppTable, "~/Dropbox/Postdoc/Papers/MGMT Paper/Supplementary Table.xlsx", col_names = TRUE, format_headers = TRUE)

