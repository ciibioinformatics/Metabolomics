library(randomForest)
library(xgboost)
met_data_train<-read.csv("/Users/bl2493/Desktop/Work/Metabolomics/RSP/metadata_a1_v2.csv",header=TRUE)
met_data_test<-read.csv("/Users/bl2493/Desktop/Work/Metabolomics/RSP/extension_data/ess_metadata_test_v2.csv",header=TRUE)


MyMWtest <- function(data,npar=TRUE,print=TRUE){
  last_col<-ncol(data)
  mw_test_results<-lapply(data[,3:last_col], function(x) wilcox.test(x~data$mecfs_status))
  #mw_test_results<-lapply(data[,5:last_col], function(x) wilcox.test(x~data$mecfs_status))
  #mw_test_results<-lapply(data[,11:last_col], function(x) wilcox.test(x~data$mecfs_status))
  
  mw_test_pvals<- vapply(mw_test_results,"[[",0, i = "p.value")
  print_mw_test_pvals<-data.frame(mw_test_pvals)
  print_mw_test_adjusted_pvals<-data.frame(p.adjust(mw_test_pvals, method="BH"))
  
  print_mw_test_all<-merge(print_mw_test_pvals,print_mw_test_adjusted_pvals, by=0)
  colnames(print_mw_test_all)<-c("metabolite","mw.test.p.value","mw.test.adjusted.p.value")
  
  return(print_mw_test_all)
}

#combined_mw_output<-MyMWtest(met_data_train)
combined_mw_output<-MyMWtest(met_data_test)


set.seed(1)
mecfs.rf <- randomForest(met_data_train[,3:19], met_data_train[,2], prox=TRUE, importance=TRUE, keep.forest = TRUE) 
importance(mecfs.rf)

rf_prediction<-predict(mecfs.rf, newdata=met_data_test)
table(rf_prediction, as.vector(met_data_test$mecfs_status))


set.seed(1)
mecfs.xg<-xgboost(data=as.matrix(met_data_train[,c(3:19)]), label=(met_data_train$mecfs_status), nrounds=1000 )  # TUBB / TUBA1A
importance_matrix<-xgb.importance(model=mecfs.xg)
print(importance_matrix)

xg_pred=predict(mecfs.xg, as.matrix(met_data_test[,3:19]))
#xg_prediction<-matrix(xg_pred, ncol=3,byrow=TRUE)
table(round(xg_pred), as.vector(met_data_test$mecfs_status))


#write.table(combined_mw_output, "/Users/bl2493/Desktop/Work/Metabolomics/RSP/metadata_a1_MW_20190528.csv",sep=",",row.names=FALSE)
write.table(combined_mw_output,"/Users/bl2493/Desktop/Work/Metabolomics/RSP/extension_data/ess_metadata_MW_20190528.csv",sep=",",row.names=FALSE)

par(new=True)
par(cex.axis=0.9)
par(cex.lab=1) # is for y-axis
par(cex.main=1.2) # is for y-axis
par(mfrow=c(2,4))
boxplot(met_data$sf_36_vitality~met_data$mecfs_status, main="sf_36_vitality", xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$sf_36_phys_funct~met_data$mecfs_status, main="sf_36_phys_funct",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$sf_36_phys_limit~met_data$mecfs_status, main="sf_36_phys_limit",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$sf_36_emo_limit~met_data$mecfs_status, main="sf_36_emo_limit",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$sf_36_emo_wellbeing~met_data$mecfs_status, main="sf_36_emo_wellbeing",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$sf_36_soc_fuct~met_data$mecfs_status, main="sf_36_soc_fuct",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$sf_36_pain~met_data$mecfs_status,main="sf_36_pain", xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$sf_36_general~met_data$mecfs_status, main="sf_36_general",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)


par(new=True)
par(cex.axis=0.9)
par(cex.lab=1) # is for y-axis
par(cex.main=1.2) # is for y-axis
par(mfrow=c(1,5))
boxplot(met_data$mfi_gen_fatigue_1~met_data$mecfs_status, main="mfi_gen_fatigue", xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$mfi_phys_fatigue_1~met_data$mecfs_status, main="mfi_phys_fatigue",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$mfi_mental_fatigue_1~met_data$mecfs_status, main="mfi_mental_fatigue",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$mfi_red_activity_1~met_data$mecfs_status, main="mfi_red_activity",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)
boxplot(met_data$mfi_red_motiv_1~met_data$mecfs_status, main="mfi_red_motiv",xlab="me/cfs status ( case = 1)", ylab="score", col=c("light blue","salmon"), font.main =  0.1)


library(ggfortify)
library(devtools)
library(ggbiplot)

pca_metadata<-prcomp(met_data[,c(13:26)])
ggbiplot(pca_metadata,ellipse=TRUE,choices=c(5,6))
