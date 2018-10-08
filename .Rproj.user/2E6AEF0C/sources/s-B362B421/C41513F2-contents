library(ggplot2)

dirPath <- '~/Dropbox/Postdoc/MGMT Project/BSAS Results/BSAS Analysis R Project'
data <- read.csv(file.path(dirPath, "MGMT_Infinium_BSAS.csv"))

#scatterplot of all data 
ggplot(data, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)")

#aggregate data
agdata_Infinium <- aggregate(data$MGMT_Infinium, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean)
agdata_BSAS <-aggregate(data$MGMT_BSAS, by=list(Patient=data$Patient, Tumor=data$Tumor, HM_status=data$HM_status), mean)
agdata <- merge(agdata_Infinium, agdata_BSAS, by = c("Patient", "Tumor", "HM_status"))
colnames(agdata)[which(names(agdata) == "x.x")] <- "MGMT_Infinium"
colnames(agdata)[which(names(agdata) == "x.y")] <- "MGMT_BSAS"
agdata_primary <- agdata[which(agdata$Tumor == 'Primary'),]
agdata_recurrence <- agdata[which(agdata$Tumor == 'Recurrence'),]

#scatterplot of aggregate data
ggplot(agdata, aes(x=MGMT_BSAS, y=MGMT_Infinium, color=HM_status)) + 
  geom_point() + 
  geom_text(aes(label=Patient), nudge_x=0.025) +
  labs(x= "Percent CpG Methylation (BSAS)", y = "Probability of Methylation (Infinium Array)")

#boxplots
boxplot(agdata$MGMT_Infinium~agdata$Tumor*agdata$HM_status, col=(c("gold","darkgreen")))
boxplot(agdata$MGMT_BSAS~agdata$Tumor*agdata$HM_status, col=(c("gold","darkgreen")))
wilcox.test(agdata_primary$MGMT_Infinium~agdata_primary$HM_status)
wilcox.test(agdata_recurrence$MGMT_Infinium~agdata_recurrence$HM_status)
wilcox.test(agdata_primary$MGMT_BSAS~agdata_primary$HM_status)
