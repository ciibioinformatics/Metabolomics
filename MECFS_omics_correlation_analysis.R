library("Hmisc")

# step 1. read metabolomics dataset
met_data<-read.csv("/Users/bl2493/Desktop/Work/Metabolomics/RSP/raw_data/RSP_merged_normalized_data_w_basic_metadata.csv",header=TRUE)
#met_data<-read.csv("/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/RSP_merged_normalized_data_w_basic_metadata_subset.csv",header=TRUE, check.names=FALSE)
fecal_species_data<-read.csv("/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/braken/RSS/RSS_species_otu_braken_reduced_mw.csv",header=TRUE, check.names=FALSE)
saliva_species_data<-read.csv("/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/braken/RSA/RSA_species_otu_braken_reduced_mw.csv",header=TRUE, check.names=FALSE)
ITS<-read.csv("/Users/bl2493/Desktop/correlation_b_f_test1.csv",header=TRUE, check.names=FALSE)

ITS_cor_out<-rcorr(as.matrix(ITS[,2:76]), as.matrix(ITS[,77:126]), type="spearman")
ITS_cor_print<-flattenCorrMatrix(ITS_cor_out$r, ITS_cor_out$P)
write.table(ITS_cor_print, "/Users/bl2493/Desktop/correlation_b_f_test_output2.csv",sep=",",row.names=FALSE)

RSS_cor_out<-rcorr(as.matrix(fecal_species_data[,3:217]), type="spearman")
RSS_cor_print<-flattenCorrMatrix(RSS_cor_out$r, RSS_cor_out$P)
write.table(RSS_cor_print, "/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/braken/RSS/RSS_interspecies_corr_20190419.csv",sep=",",row.names=FALSE)

RSA_cor_out<-rcorr(as.matrix(fecal_species_data[,3:68]), type="spearman")
RSA_cor_print<-flattenCorrMatrix(RSA_cor_out$r, RSA_cor_out$P)
write.table(RSA_cor_print, "/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/braken/RSA/RSA_interspecies_corr_20190419.csv",sep=",",row.names=FALSE)

RSP_cor_out<-rcorr(as.matrix(met_data[,9:47]), type="spearman")
RSP_cor_print<-flattenCorrMatrix(RSP_cor_out$r, RSP_cor_out$P)
write.table(RSP_cor_print, "/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/RSP_subset_corr_20190419.csv",sep=",",row.names=FALSE)


combined_data<-merge(met_data,fecal_species_data, by=c("SampleName", "SampleType"))
#RSP_RSS_cor_out<-rcorr(as.matrix(combined_data[,9:262]), type="spearman")
#RSP_RSS_cor_print<-flattenCorrMatrix(RSS_cor_out$r, RSS_cor_out$P)
#write.table(RSS_cor_print, "/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/RSP_RSS_subset_corr_20190419.csv",sep=",",row.names=FALSE)

#data<-combined_data[which(combined_data$SampleType=='Case'),]
#data<-combined_data[which(combined_data$SampleType=='Control'),]
#data <- droplevels(data)
#combined_data<-data
#combined_data<-merge(saliva_species_data,fecal_species_data, by=c("SampleName", "SampleType"))
#combined_data<-merge(saliva_genus_data,fecal_genus_data, by=c("SampleName", "SampleType"))
#combined_data<-merge(met_data,humann_data, by=c("SampleName", "SampleType"))


library("psych")
out_data<-data.frame()
for (i in 9:896){
  for(j in 897:1111){
    res <- corr.test(combined_data[,i], combined_data[,j],method="spearman", adjust="BH") 
    temp_out<-cbind(colnames(combined_data)[i], colnames(combined_data[j]), res$r, res$p) 
    out_data<-rbind(out_data,temp_out)
    cat(colnames(combined_data)[i],"\t", colnames(combined_data)[j],"\t",res$r,"\t",res$p,"\n", file="/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/RSP_RSS_all_comp_sub_species_corr_20190419.csv", append=TRUE)
   }
}
#colnames(out_data)<-c("metabolite","species","rho","pval")
#out_data$adj_pval<-data.frame(p.adjust(out_data$pval, method="BH"))
#write.table(out_data, "/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/RSP_RSS_subset_all_species_corr_20190324.csv",sep=",",row.names=FALSE)

# o2pls analysis
library("OmicsPLS")
o2_met_data<-ifelse((grepl("Control", combined_data[,2])),0,1)
o2_fecal_species_data<-ifelse((grepl("Control", fecal_species_data[,2])),0,1)

o2_met_data<-as.data.frame(met_data[,9:896])
o2_fecal_species_data<-as.data.frame(fecal_species_data[,3:217])

########################################################################################################
# with all dataset
#o2_met_data<-cbind(o2_met_data,combined_data[,c(9:896)])
#o2_fecal_species_data<-cbind(o2_fecal_species_data,fecal_species_data[,c(3:217)])

# with subset
#o2_met_data<-cbind(o2_met_data,combined_data[,c(9:47)])
#o2_fecal_species_data<-cbind(o2_fecal_species_data,combined_data[,c(48:123)])

########################################################################################################
CV1<-crossval_o2m_adjR2(o2_met_data, o2_fecal_species_data, 1:3, c(0,1,5,10), c(0,1,5,10),
                        nr_folds = 2, nr_cores = 4)
CV2<-crossval_o2m(o2_met_data, o2_fecal_species_data,1:2, 0:2, 9:11,
                        nr_folds = 10, nr_cores = 4)
fit=o2m(o2_met_data,o2_fecal_species_data,2,2,11)
summary(fit)
#x_loading<-loadings(fit,loading_name = c("Xjoint"))
#y_loading<-loadings(fit,loading_name = c("Yjoint"))
x_loading<-loadings(fit,loading_name = c("XOrth"))
y_loading<-loadings(fit,loading_name = c("YOrth"))
#################################################
##### plotting o2pls 
library(magrittr)
library(ggplot2)
library(gridExtra)

met_names<-colnames(o2_met_data)
met_names[1:71]<-c("OL")
met_names[72:551]<-c("CL")
met_names[552:637]<-c("PM")
met_names[638:888]<-c("BA")
met_names

orig_spec_names<-colnames(o2_fecal_species_data)
spec_names<-ifelse(grepl("[[]",orig_spec_names),as.factor(orig_spec_names),'1')

#################################################
# plot for metabolomics data
#################################################
metabo_p<-plot(fit, loading_name="Xorth",label ="c",i=1, j=2, alpha=0.5) + 
  theme_bw() + coord_fixed(ratio = 1, xlim=c(-.3,.3),ylim=c(-.2,.2)) +
  geom_point( # Set color and size
    aes(col=met_names, shape=met_names, show.legend = TRUE)
  ) +
labs(title = "Metabolite Orth loadings",
     x = "First Orth Loadings", y = "Second Orth Loadings") +
  theme(legend.position="right") +
  theme(plot.title = element_text(face='bold'),
      legend.title=element_text(face='bold')) +
geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

#################################################
# plot for metagenomic data
#################################################
fecal_p<- plot(fit, loading_name="Yorth",label ="c",i=1, j=2, alpha=0.5) + 
  theme_bw() + coord_fixed(ratio = 1, xlim=c(-.3,.3),ylim=c(-.2,.2)) +
  #geom_point(col=spec_names,shape=spec_names)+
  geom_point(col="grey",alpha=0.5, show.legend = TRUE)+
  
  labs(title = "Metagenome Orth loadings",
       x = "First Orth Loadings", y = "Second Orth Loadings") +
  theme(legend.position="right") +
  theme(plot.title = element_text(face='bold'),legend.title=element_text(face='bold')) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

## Finally plot both plots in one figure.
grid.arrange(metabo_p, fecal_p, ncol=2)


####################################
MyLogTransformation <- function(x,npar=TRUE, print=TRUE){
  last_col<-ncol(x)
  processed_data<-cbind(x[,1:2], log10(x[,3:last_col]))
  return(processed_data)
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

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

combined_mw_output<-MyMWtest(saliva_species_data)
combined_mw_output<-MyMWtest(fecal_species_data)
#write.table(combined_mw_output, "/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/braken/RSS/RSS_species_mw_20190419.csv",sep=",",row.names=FALSE)
#write.table(combined_mw_output, "/Users/bl2493/Desktop/Work/Metabolomics/RSP_RSA_RSS_correlation/braken/RSA/RSA_species_mw_20190419.csv",sep=",",row.names=FALSE)

boxplot(fecal_species_data$`[Eubacterium]_rectale`~met_data$SampleType, col=c("light blue","salmon"), font.main =  0.1)

#library(PMA)
require(GGally)
require(CCA)

install.packages("ade4")
library(ade4)
cca_met_data<-as.data.frame(met_data[,9:896])
cca_fecal_species_data<-as.data.frame(fecal_species_data[,3:217])

cc1 <- cc(as.matrix(cca_met_data), as.matrix(cca_fecal_species_data))
dudi1 <- dudi.pca(cca_met_data, scale = TRUE, scan = FALSE, nf = 3)
dudi2 <- dudi.pca(cca_fecal_species_data, scale = FALSE, scan = FALSE, nf = 2)
coin1 <- coinertia(dudi1,dudi2, scan = FALSE, nf = 2)
coin1


if(adegraphicsLoaded()) {
  g1 <- s.arrow(coin1$l1, plab.cex = 0.5)
  g2 <- s.arrow(coin1$c1, plab.cex = 0.5)
  g3 <- s.corcircle(coin1$aX, plot = FALSE)
  g4 <- s.corcircle(coin1$aY, plot = FALSE)
  cbindADEg(g3, g4, plot = TRUE)
  g5 <- plot(coin1)
  
} else {
  s.arrow(coin1$l1, clab = 0.5)
  s.arrow(coin1$c1, clab = 0.5)
  par(mfrow = c(1,2))
  s.corcircle(coin1$aX)
  s.corcircle(coin1$aY)
  par(mfrow = c(1,1))
  plot(coin1)
}

par(mfrow = c(1,2))
s.arrow(coin1$l1, clab = 0.7, xlim=c(-0.1,0.1),ylim=c(-0.1,0.1), boxes=FALSE)
s.arrow(coin1$c1, clab = 0.7, xlim=c(-0.1,0.1),ylim=c(-0.1,0.1), boxes=FALSE)

summary(coin1)



library(ccrepe) 

data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10)
data.rowsum <- apply(data,1,sum)
data.norm <- data/data.rowsum
testdata <- data.norm
dimnames(testdata) <- list(paste("Sample",seq(1,10)),paste("Feature",seq(1,4)))
ccrepe.results  <-ccrepe  (x=testdata, iterations=20, min.subj=10)
ccrepe.results.nc.score <- ccrepe(x=combined_data[,9:896],y=combined_data[,897:1111],iterations=20,min.subj=10,sim.score=nc.score)


ccrepe.results
ccrepe.results.nc.score


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue")

devtools::install_github("borenstein-lab/MIMOSA/mimosa")