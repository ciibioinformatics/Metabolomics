# libraries
library(randomForest)
library(ggplot2)
library("ggpubr")
library(data.table)

library(ggplot2)
library(ggsignif)

#, check.names=FALSE
data<-read.csv("/Users/bl2493/Desktop/Work/Metabolomics/zika/zika_master_merged_datasets_final_w_meta.csv",header=TRUE)
#data<-read.csv("/Users/leebo/Desktop/Work/Metabolomics/zika/zika_master_merged_datasets_final_w_meta.csv",header=TRUE)

#################################################################
# My Functions 
MyLogTransformation <- function(x, min_factor, npar=TRUE, print=TRUE){
  last_col<-ncol(x)
  processed_data<-cbind(x[,1:5], log10(x[,6:last_col]*min_factor+10))
  return(processed_data)
}

MyScaling <- function (data, npar=TRUE, print=TRUE){
  last_col<-ncol(data)
  
  # calculate per group stdev & divide data with control's stdv.
  print_stdv<-aggregate(data[,6:last_col], list(data$SampleType), sd)
  index<-ifelse((grepl("control", print_stdv[1,1]) || grepl("Control", print_stdv[1,1])),1,2)
  for (i in 6:last_col){ data[,i]<-data[,i]/print_stdv[index,i-4]}

  return(data)
}


MyMWtest <- function(data,npar=TRUE,print=TRUE){
  last_col<-ncol(data)
  compounds<-as.data.frame(colnames(data[,c(6:last_col)]))
  mw_test_pvals<-sapply(6:ncol(data), function(i){ wilcox.test(data[,i]~data$SampleType, exact=FALSE )$p.value})
  mw_test_adj_pvals<-data.frame(p.adjust(mw_test_pvals, method="BH"))
  mw_test_w<-sapply(6:ncol(data), function(i){as.numeric(wilcox.test(data[,i]~data$SampleType, exact=FALSE )$statistic)})
  
  print_mw_test_all<-cbind(compounds, as.data.frame(mw_test_pvals),
                        as.data.frame(mw_test_adj_pvals),
                        as.data.frame(mw_test_w)
                        )
  colnames(print_mw_test_all)<-c("metabolite","mw.test.p.value","mw.test.adjusted.p.value","mw.w.statistic")
  
  return(print_mw_test_all)
}

MyLNRegression<-function(data, npar=TRUE,print=TRUE){
  
  last_col<-ncol(data)
  lnr.out<-lapply(data[,6:last_col], function (x) gamlss(x~data$SampleType, family=LOGNO2, data=data))
  lnr.coef<-data.frame(lapply(lnr.out[], function(f) summary(f,save=TRUE)$coef.table[2,c(1,4)]))
  #lnr.pval<-data.frame(lapply(lnr.out[], function(f) summary(f,save=TRUE)$pvalue[2]))
  lnr_adjusted_pvals<-data.frame(p.adjust(lnr.coef[2,], method="BH"))
  
  lnr_deviance<-data.frame(lapply(lnr.out[], function(f) f$G.deviance))
  lnr_aic<-data.frame(lapply(lnr.out[], function(f) f$aic))
  
  print_lnr_coef<-t(lnr.coef)
  print_lnr_deviance<-t(lnr_deviance)
  print_lnr_aic<-t(lnr_aic)

  print_lnr_pvals<-merge(print_lnr_coef,lnr_adjusted_pvals, by=0)
  colnames(print_lnr_pvals)<-c("metabolite", "LNR.coef(mu)", "LNR.pval(mu)", "LNR.adj.pval(mu)")
  print_lnr_stats<-merge(print_lnr_deviance,print_lnr_aic,by=0)
  colnames(print_lnr_stats)<-c("metabolite", "LNR.deviance", "LNR.AIC")
  
  print_lnr_all<-merge(print_lnr_pvals,print_lnr_stats, by="metabolite")
  
  return(print_lnr_all)
}

MyRegression<-function(data, npar=TRUE,print=TRUE){
  
  last_col<-ncol(data)
  #glm.out<-lapply(data[,6:last_col], function (x) glm(x~data$SampleType, family=gaussian(identity)))
  glm.out<-lapply(data[,6:last_col], function (x) glm(x~data$SampleType, family=Gamma("log")))
  
  glm_coef<-data.frame(lapply(glm.out[], function(f) summary(f)$coefficients[2,c(1,4)]))
  glm_deviance<-data.frame(lapply(glm.out[], function(f) summary(f)$deviance))
  glm_aic<-data.frame(lapply(glm.out[], function(f) summary(f)$aic))
  glm_adjusted_pvals<-data.frame(p.adjust(glm_coef[2,], method="BH"))
  
  print_glm_coef<-t(glm_coef)
  print_glm_deviance<-t(glm_deviance)
  print_glm_aic<-t(glm_aic)
  
  print_glm_pvals<-merge(print_glm_coef,glm_adjusted_pvals, by=0)
  colnames(print_glm_pvals)<-c("metabolite", "LR.coef", "LR.pval", "LR.adj.pval")
  print_glm_stats<-merge(print_glm_deviance,print_glm_aic,by=0)
  colnames(print_glm_stats)<-c("metabolite", "LR.deviance", "LR.AIC")
  
  print_glm_all<-merge(print_glm_pvals,print_glm_stats, by="metabolite")
  
  return(print_glm_all)
}

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

#################################################################
# Step 0 : descriptive stats
library(psych)
desc_stats_summary<-describeBy(data,group=data$SampleType)
desc_stats_print<-do.call("rbind",desc_stats_summary)
write.table(desc_stats_print, "/Users/bl2493/Desktop/Work/Metabolomics/zika/analysis_20190829/descriptive_stats.csv",sep=",",row.names=TRUE)

# Step 1 : log-transformation
data<-MyLogTransformation(data,1)
#desc_stats_summary<-describeBy(data,group=data$SampleType)
#desc_stats_print<-do.call("rbind",desc_stats_summary)
#write.table(desc_stats_print, "/Users/bl2493/Desktop/Work/Metabolomics/zika/analysis_20190829/descriptive_stats_log_transformed.csv",sep=",",row.names=TRUE)

# Step 2 : rescaling  (Divide by control's stdv)
data<-MyScaling(data)
#desc_stats_summary<-describeBy(data,group=data$SampleType)
#desc_stats_print<-do.call("rbind",desc_stats_summary)
#write.table(desc_stats_print, "/Users/bl2493/Desktop/Work/Metabolomics/zika/analysis_20190829/descriptive_stats_log_transformed_rescaled.csv",sep=",",row.names=TRUE)

# Step 3 : subset data for wilcox 
#data<-data[which(data$SampleType=='MCS-ZIKV+'|data$SampleType == 'Control'),]  # regarding zika
data<-data[which(data$SampleType=='MCS+ZIKV'|data$SampleType == 'Control'),]  # regarding zika
#data<-data[which(data$SampleType=='MCS-ZIKV+'|data$SampleType == 'MCS+ZIKV'),]
data$SampleType<-factor(data$SampleType)

# Step 4 : run MW analysis
combined_mw_output<-MyMWtest(data)
write.table(combined_mw_output, "/Users/bl2493/Desktop/Work/Metabolomics/zika/analysis_20190829/mcs_vs_control_mw_w_w_stats_20190910.csv",sep=",",row.names=FALSE)

# Step 5: correlation analysis to identify predictors

library(corrplot)
library(RColorBrewer)
library("Hmisc")
M<-rcorr(as.matrix(data[,c(6:1043)]))
corrplot(M$r, type="upper", order="hclust", 
         p.mat = M$P, sig.level = 0.01, insig = "blank")
corr_out<-flattenCorrMatrix(M$r, M$P)
corr_out<-corr_out[which(corr_out$p<0.05),]
write.table(corr_out, "/Users/bl2493/Desktop/Work/Metabolomics/zika/analysis_20190829/zika_cotrol_corr_out_20190911.csv",sep=",",row.names=FALSE)

# Step 6 : run linear regression analysis
#model <- lm(CPD00001 ~ SampleType, data = data)
#summary(model)
combined_lr_output<-MyRegression(data)
write.table(combined_lr_output, "/Users/bl2493/Desktop/Work/Metabolomics/zika/analysis_20190829/mcs_vs_control_gamma_regression_20190911.csv",sep=",",row.names=FALSE)

combined_lnr_output<-MyLNRegression(data)
write.table(combined_lnr_output, "/Users/bl2493/Desktop/Work/Metabolomics/zika/analysis_20190829/mcs_vs_control_log_norm_regression_20190911.csv",sep=",",row.names=FALSE)

# biomarkers (zika vs control)
CPD00073
CPD00039
CPD00074
CPD00646
CPD00048
CPD00291

zika_markers_data<-data[,c(2,78,44,79,651,53,296)]


#Step 1. Logistic Regression
library(caret)

train_control <- trainControl(method="boot", number=100) # cv 
ckd.lr <- glm(classification~.,family=binomial(link='logit'),data=data_train, control=glm.control(maxit=1) )
summary(ckd.lr)
ckd.lr.pred<-predict(ckd.lr)
ckd.lr.result <- confusionMatrix(as.factor(ckd.lr.pred), as.factor(zika_markers_data$SampleType))

library('ROCR')
zika_pred_model<-glm(SampleType~., family=binomial(logit), data=zika_markers_data)
zika_p <- predict(zika_pred_model, newdata=zika_markers_data, type="response")
zika_pr <- prediction(zika_p, zika_markers_data$SampleType)
zika_prf <- performance ( zika_pr, measure = "tpr", x.measure="fpr")
plot(zika_prf)

auc <- performance(zika_pr, measure = "auc")
auc <-auc@y.values[[1]]
auc

# Step 5: visualization

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
par(new=TRUE)
par(cex.axis=0.9)
par(cex.lab=1) # is for y-axis
par(cex.main=1.2) # is for y-axis
par(mfrow=c(1,1))
ggplot(data, aes(x=SampleType, y=CPD00073)) + 
  geom_boxplot(fill = c("grey","light blue","salmon")) +
  stat_compare_means(comparisons = list(c("Control", "MCS-ZIKV+"),c("Control","MCS+ZIKV"),c("MCS-ZIKV+","MCS+ZIKV")), symnum.args=symnum.args)#+
  #stat_compare_means(label.y = 50) # Add pairwise comparisons p-value
