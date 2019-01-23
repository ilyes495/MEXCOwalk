#library(maxstat)
library(survival)
library(survminer)

args=(commandArgs(TRUE))
cluster_filename=args[1]
clinical_filename=args[2]
plot_filename=args[3]
pvalue_filename = args[4]


survival.table <- read.table(clinical_filename, header = TRUE, sep = '\t')
if(data_type == 'old_data'){
  needed.cols = c("PatientID","Survival","Death") #, "age", "gender")  
  survival.table2 <- survival.table[,needed.cols]
  rownames(survival.table2) <- survival.table2[,'PatientID']
  surv_table <- survival.table2[,2:3]
  
  surv_table$GROUP = readLines(cluster_filename)
  num_clusters = length(unique(surv_table$GROUP ))
  
  kmsurvival <- survfit(Surv(surv_table$Survival,surv_table$Death) ~ strata(surv_table$GROUP))
  
  test.diff <- survdiff(Surv(surv_table$Survival,surv_table$Death) ~  surv_table$GROUP, rho=0)
  
} else{
  needed.cols = c("patient_id","time_to_event","event") #, "age", "gender")
  survival.table2 <- survival.table[,needed.cols]
  rownames(survival.table2) <- survival.table2[,'patient_id']
  surv_table <- survival.table2[,2:3]
  
  surv_table$GROUP = readLines(cluster_filename)
  num_clusters = length(unique(surv_table$GROUP ))
  
  kmsurvival <- survfit(Surv(surv_table$time_to_event,surv_table$event) ~ strata(surv_table$GROUP))
  
  test.diff <- survdiff(Surv(surv_table$time_to_event,surv_table$event) ~  surv_table$GROUP, rho=0)
}

p.value <- round(1-pchisq(test.diff$chi,df=num_clusters-1),10)
print(p.value)
write(p.value, file = pvalue_filename)
pdf(paste("~/Dropbox/Tunde/analysis/new_results/plots/", plot_filename, sep=""),width=7,height=5)
#plot(kmsurvival, main = 'Breast Survival Data', xlab = 'Time (Days)', ylab = 'Survival Probability', col = c('black', 'gray', 'blue'), lwd = 3, cex.lab=1.2, cex.axis=1.2)
#text(1000,0.05,paste(" p.value=", p.value), cex = 1.2)
#legend(x="topright", legend = c('cluster 1', 'cluster 2', 'cluster 3'), col=c("black","gray","blue"), pch=1, bty = "n", cex = 1.2)
#dev.off()

ggsurvplot(kmsurvival)
dev.off()
#ggsurvplot(kmsurvival, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
#           legend.labs=c("cluster1", "cluster2", "cluster3"), legend.title="Groups",  
#           palette=c("dodgerblue2", "orchid2"), 
#           main="Kaplan-Meier Curve for Colon Cancer Survival", 


#coxph(Surv(surv_table$time_to_event,surv_table$event) ~   GROUP + gender + age, surv_table)


