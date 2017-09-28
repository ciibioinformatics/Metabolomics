# libraries
library(randomForest)
library(dunn.test)
library(ggplot2)
library("ggpubr")

# Case vs. Control data

#, check.names=FALSE
combined_data<-read.csv("/Users/andrewschultz/Desktop/Projects/MECFS_spinal_fluid/MECFS_CSF_w_cytokine.csv",header=TRUE)
#combined_data<-read.csv("/Volumes/KINGSTON/MECFS_spinal/data/MECFS_CSF_case_vs_normal_raw_data.csv",header=TRUE)
#combined_data<-read.csv("F://MECFS_spinal/data/MECFS_CSF_case_vs_normal_raw_data.csv",header=TRUE)

#combined_data<-read.csv("/Volumes/KINGSTON/MECFS_spinal/data/MECFS_CSF_case_vs_MS_raw_data.csv",header=TRUE)
#combined_data<-read.csv("/Volumes/KINGSTON/MECFS_spinal/data/MECFS_CSF_MS_vs_normal_raw_data.csv",header=TRUE)

#meta_data<-read.csv("/Volumes/KINGSTON/PIA/PIA_metadata.csv", header=TRUE)
meta_data<-read.csv("/Users/andrewschultz/Desktop/Projects/PIA/data/PIA_metadata.csv", header=TRUE)

#################################################################
# My Functions 
MyDescStats <- function(data, npar=TRUE,print=TRUE){
  last_col<-ncol(data)
  metNames<-colnames(data)
  desc_stats<-data.frame()
  
  median_output<-lapply(data[,3:last_col], function (x) aggregate(formula=x~SampleType, data=combined_data,FUN=median))
  print_desc_stats<-lapply(data[,3:last_col], function (y) median_output$y$x)
  
  return(print_desc_stats)
}

MyImputation <- function(x, cutoff, npar=TRUE,print=TRUE) {
  imputed<-x[, -which(colMeans(is.na(x)) > cutoff )]
 # imputed<-x[, which(colMeans(is.na(x)) > cutoff )]
  return(imputed)
}

MyFillinNA <- function(data, npar=TRUE, print=TRUE) {
  last_col<-ncol(data)
  for ( i in 3:last_col){
    min_val<-min(na.omit(data[,i]))/2
    data[,i]<-ifelse(is.na(data[,i]),min_val,data[,i])
  }

  return(data)
}

MyNormalization <- function(x, npar=TRUE, print=TRUE) {
    last_col<-ncol(x)
    per_sample_norm<-cbind(x[,1:2], x[,3:last_col]/rowSums(x[,3:last_col]))
    return(per_sample_norm)
}

MyFindMinFactor <- function(x, npar=TRUE, print=TRUE){
  last_col<-ncol(x)
  min_val<-min(x[,3:last_col])
  return(min_val)
}

MyLogTransformation <- function(x, min_factor, npar=TRUE, print=TRUE){
  last_col<-ncol(x)
  processed_data<-cbind(x[,1:2], log10(x[,3:last_col]*min_factor))
  return(processed_data)
}

MyTtest <- function(data,npar=TRUE,print=TRUE){
  last_col<-ncol(data)
  t_test_results<-lapply(data[,3:last_col], function(x) t.test(x~data$SampleType))
  t_test_pvals<- vapply(t_test_results,"[[",0, i = "p.value")
  print_t_test_pvals<-data.frame(t_test_pvals)
  print_t_test_adjusted_pvals<-data.frame(p.adjust(t_test_pvals, method="BH"))
  
  print_t_test_all<-merge(print_t_test_pvals,print_t_test_adjusted_pvals, by=0)
  colnames(print_t_test_all)<-c("metabolite","t.test.p.value","t.test.adjusted.p.value")
  
  return(print_t_test_all)
}

MyMWtest <- function(data,npar=TRUE,print=TRUE){
  last_col<-ncol(data)
  mw_test_results<-lapply(data[,3:last_col], function(x) wilcox.test(x~data$SampleType))
  mw_test_pvals<- vapply(mw_test_results,"[[",0, i = "p.value")
  print_mw_test_pvals<-data.frame(mw_test_pvals)
  print_mw_test_adjusted_pvals<-data.frame(p.adjust(mw_test_pvals, method="BH"))
  
  print_mw_test_all<-merge(print_mw_test_pvals,print_mw_test_adjusted_pvals, by=0)
  colnames(print_mw_test_all)<-c("metabolite","mw.test.p.value","mw.test.adjusted.p.value")
  
  return(print_mw_test_all)
}

MyKWtest <- function(data,npar=TRUE,print=TRUE){
  last_col<-ncol(data)
  metNames<-colnames(data)
  kw_pvals<-data.frame()
  
  for (i in 3:ncol(data)){
    curr_met<-metNames[i]
    kw_test_results<-kruskal.test(data[,i]~data$SampleType)
    pvals<-kw_test_results$p.value
  
    post_hoc<-dunn.test(data[,i],g=data$SampleType, method="bh")
    kw_pvals<-rbind(kw_pvals,cbind(curr_met,pvals,post_hoc$P.adjusted[1],post_hoc$P.adjusted[2],post_hoc$P.adjusted[3]))
    
  }  
  colnames(kw_pvals)<-c("metabolite","kw.p.value",post_hoc$comparisons[1],post_hoc$comparisons[2],post_hoc$comparisons[3])
  return(kw_pvals)  
}

MyANOVAtest <- function(data,npar=TRUE,print=TRUE){

  metNames<-colnames(data)
  anova_pvals<-data.frame()
  
  for (i in 3:ncol(data)){
    curr_met<-metNames[i]
    anova_result<-aov(data[,i]~as.factor(data$SampleType), data=data) 
    pval<-anova(aov(data[,i]~as.factor(data$SampleType)))[['Pr(>F)']][1]
    #tukey_result<-TukeyHSD(x=anova_result, 'as.factor(data$SampleType)', conf.level=0.95)
    posthoc_out<-pairwise.t.test(data[,i], as.factor(data$SampleType), p.adjust="none", pool.sd = T) 
    anova_pvals<-rbind(anova_pvals,cbind(curr_met,pval,posthoc_out$p.value[1,1], posthoc_out$p.value[2,1], posthoc_out$p.value[2,2]))
  }
  
  adjusted_pvals<-as.numeric(levels(anova_pvals$pval))[anova_pvals$pval] 
  
  print_anova_adjusted_pval<-p.adjust(adjusted_pvals, method="BH")
  print_anova_all<-merge(anova_pvals,print_anova_adjusted_pval, by=0)
  colnames(print_anova_all)<-c("rownum","metabolite","anova.p.value","anova.post.hoc.1","anova.post.hoc.2","anova.post.hoc.3","anova.adjusted.p.value")
  print_anova_all<-subset(print_anova_all, select = -c(rownum) )
  
return(print_anova_all)
  #return(posthoc_out)
}

MyRandomForest<-function(primary_data, npar=TRUE, print=TRUE){
  last_col<-ncol(primary_data)

  # calculate RF importnace
  rf.out<-lapply(primary_data[,3:last_col], function (x) randomForest(primary_data$SampleType~x, ntree=10000, importance=TRUE))
  rf_importance<-data.frame(lapply(rf.out, function(f) importance(f)[1]))
  print_rf_importance<-cbind(rank(-rf_importance),t(rf_importance))

  print_rf_all<-cbind(rownames(print_rf_importance), print_rf_importance)
  colnames(print_rf_all)<-cbind("metabolite","RF.importance.Rank","RF.importance")
  
  return(print_rf_all)
}

MyLogisticRegression<-function(data, npar=TRUE,print=TRUE){

  last_col<-ncol(data)

  # calculate per group stdev for odds ratio adjustment & divide data with control's stdv.
  print_stdv<-aggregate(data[,3:last_col], list(data$SampleType), sd)
  index<-ifelse(((grepl("control", print_stdv[1,1]) || grepl("Control", print_stdv[1,1])) && grepl("Non neuro", print_stdv[1,1])),1,2)
  for (i in 3:last_col){ data[,i]<-data[,i]/print_stdv[index,i-1]}
 
  # merge dataset with metadata to add covariates to the analysis 
  # (sex, age, race/ethnicity, site, season, BMI)  
  merged_data<-merge(data,meta_data,by="SampleName")

  #glm.out<-lapply(merged_data[,3:last_col], function (x) glm(merged_data$SampleType ~ x +factor(merged_data$CC) + factor(merged_data$Sex), family=binomial(logit)))
  glm.out<-lapply(data[,3:last_col], function (x) glm(merged_data$SampleType ~ x, family=binomial(logit)))
  glm_pvals<-data.frame(lapply(glm.out, function(f) summary(f)$coefficients[2,4]))
  print_glm_pvals<-t(glm_pvals)
  colnames(print_glm_pvals)<-c("LR.p.value")
  
  # OR = 1, OR > 1, higher in case, OR < 1 lower in case 
  glm_odd_ratio<-data.frame(lapply(glm.out, function(f) exp(coef(f)[2])))
  print_glm_odd_ratio<-t(glm_odd_ratio)
  colnames(print_glm_odd_ratio)<-c("LR.odd.ratio")
  
  print_lr_all<-merge(print_glm_odd_ratio,print_glm_pvals, by=0)
  colnames(print_lr_all)<-c("metabolite","lr.odd.ratio","lr.p.value")
  print_lr_all$directionality<-ifelse(print_lr_all$lr.odd.ratio >= 1, "Increased","Decreased")
  
  return(print_lr_all)
}


MyMerge<-function(x,y,npar=TRUE,print=TRUE){
  combined_data1<-merge(x, y, by="metabolite")
  #combined_data2<-merge(combined_data1,z,by="metabolite")
  return(combined_data1)
}
#################################################################
# Step 0 : descriptive stats
desc_stats_summary<-MyDescStats(combined_data)

# step 1 : select metabolites with more than cut-off percentages (default:50%)
combined_data<-MyImputation(combined_data, 0.75)
#write.table(combined_data, "/Volumes/KINGSTON/MECFS_spinal/combined_data_w_50%_missing_values.csv",sep=",",row.names=FALSE)

# step 2: fill-in NA 
combined_data<-MyFillinNA(combined_data)

# step 3: per-sample normaliztion 
combined_data<-MyNormalization(combined_data)

# step 4 : find minimum value to use as normalization factor (for AYASDI/LefSe)
MyFindMinFactor(combined_data)

# step 5: Transformation - multiply by min factor to move distribution to pos domain & do log-transformation
combined_data<-MyLogTransformation(combined_data, 1.0e+07)

#write.table(combined_data, "/Volumes/KINGSTON/MECFS_spinal/MECFS_CSF_Case_vs_ND_preprocessed.csv",sep=",",row.names=FALSE)
#write.table(combined_data, "/Volumes/KINGSTON/MECFS_spinal/MECFS_CSF_Case_vs_MS_preprocessed.csv",sep=",",row.names=FALSE)
#write.table(combined_data, "/Volumes/KINGSTON/MECFS_spinal/MECFS_CSF_MS_vs_ND_preprocessed.csv",sep=",",row.names=FALSE)


# step 6: t-test 
met_data<-combined_data[,c(1:226)]
met_data<-met_data[which(met_data$SampleType=='MS control' |
                           met_data$SampleType == 'P1/classic CFS'),]

#met_data<-met_data[which(met_data$SampleType=='Non neuro disease control' |
#                           met_data$SampleType == 'P1/classic CFS'),]
met_data$SampleType <- factor(met_data$SampleType)

combined_mw_output<-MyMWtest(met_data)

#combined_mw_output<-MyMWtest(combined_data)
#write.table(combined_anova_output, "/Volumes/KINGSTON/MECFS_CSF_3grp_anova_output_20170807_2.csv",sep=",",row.names=FALSE)

combined_rf_output<-MyRandomForest(met_data)

# step 7: Logistic Regression
#combined_lr_output<-MyLogisticRegression(combined_data)

# merge output tables for printing
print_combined<-MyMerge(combined_mw_output, combined_rf_output)
#print_combined<-MyMerge(combined_mw_output, combined_rf_output, combined_lr_output)

# Write final output tables  - CaseIBS
#write.table(print_combined, "/Users/andrewschultz/Desktop/Projects/MECFS_spinal_fluid/C_case_vs_ND_mw_rf_manuscript.csv",sep=",",row.names=FALSE)
write.table(print_combined, "/Users/andrewschultz/Desktop/Projects/MECFS_spinal_fluid/C_case_vs_MS_mw_rf_manuscript.csv",sep=",",row.names=FALSE)
#write.table(print_combined, "/Volumes/KINGSTON/MECFS_spinal/analysis/MS_vs_Normal_no_covariate_20170830_validation.csv",sep=",",row.names=FALSE)
#write.table(print_combined, "/Volumes/KINGSTON/MECFS_spinal/analysis/3grp_20170818.csv",sep=",",row.names=FALSE)


# step 7: Random Forest
# MECFS vs ND
combined_rf_data<-combined_data[c(1,2,54,220,19,53,81,6,227,88,200,199,182,74)]
#> colnames(combined_rf_data)
#[1] "SampleName"              "SampleType"              "mannose"                
#[4] "Arachidonyldopamine"     "citramalic.acid"         "maleimide"              
#[7] "threonic.acid"           "X2.hydroxybutanoic.acid" "Acetylcarnitine"        
#[10] "xylitol"                 "Cytidine"                "Cytosine"               
#[13] "L.beta.Homoleucine"      "shikimic.acid"         

#combined_rf_data<-combined_data[c(1,2,227,19,54,31,198,241,215,67,235,150)]
#> colnames(combined_rf_data)
#[1] "SampleName"                                  "SampleType"                                 
#[3] "Acetylcarnitine"                             "citramalic.acid"                            
#[5] "mannose"                                     "glyceric.acid"                              
#[7] "Diethanolamine"                              "X1.Hexadecanoyl.sn.glycero.3.phosphocholine"
#[9] "Betaine"                                     "phosphoethanolamine"                        
#[11] "X3.Hydroxypyridine"                          "FA..10.0...capric.acid."                    


combined_rf_output<-MyRandomForest(combined_rf_data)
