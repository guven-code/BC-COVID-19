#############################################################################
# Gene expression analysis of expression data from human bone
# marrow hematopoietic stem cells
# Organism	Homo sapiens
# Dataset: GSE32719 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32719)
# Pang WW, Price EA, Sahoo D, Beerman I et al. Human bone marrow hematopoietic stem cells 
# are increased in frequency and myeloid-biased with age.
# Proc Natl Acad Sci U S A 2011 Dec 13;108(50):20012-7.
# PMID: 22123971 (http://www.ncbi.nlm.nih.gov/pubmed/22123971)
# R code: Emine Guven
##############################################################################

library(gplots)
library(biomaRt)
library(matrixStats)
library(Biobase)
library(GEOquery)
library(rgr)

# # load series and platform data from GEO
# 
# gset <- getGEO("GSE28735", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 
# # set parameters and draw the plot
# 
# dev.new(width=4+dim(gset)[[2]]/5, height=6)
# par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
# title <- paste ("GSE28735", '/', annotation(gset), " selected samples", sep ='')
# boxplot(exprs(gset), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
# 
# GSE28735exprs<-exprs(gset)

#saveRDS(GSE28735exprs, "gset28735.rds")
#GSEdata<-readRDS("gset28735.rds")
GSEdata<-read.csv("GSE156445raw.csv",stringsAsFactors = FALSE)



#samples<-colnames(GSEdata)

geneNames<-GSEdata$probe_ID
samples<-colnames(GSEdata[,-1])

GSEdata.mat<-data.matrix(GSEdata,rownames.force = NA)

DEdata.mat=GSEdata.mat[,-1]
#DEdata.mat[,1]<-GeneNames
row.names(DEdata.mat)<-geneNames
head(DEdata.mat,5)


library(functional)
GSEdata<-DEdata.mat[apply(DEdata.mat, 1, Compose(is.finite, all)),]


# Check the loaded dataset
dim(GSEdata) # Dimension of the dataset
head(GSEdata) # First few rows
tail(GSEdata) # Last few rows

rownames(GSEdata)
colnames(GSEdata)

#raw_data=data.matrix(GSEdata,rownames.force = NA)
#rownames(raw_data)=geneIDs

pdf("plots/histogramGSE.pdf")
hist(GSEdata, col = "gray", main="GSE156445 - Histogram")
dev.off()


# Check the behavior of logdata
pdf("plots/histogramLogGSE.pdf")

title <- paste ("GSE156445", "/", "GPL17586", sep ="")
h<-hist(log2(GSEdata), xlab= "log2(GSE156445)", col = "gray", ylim=c(0,length(log2(GSEdata))/20),
        breaks=80,main="GSE156445 - Histogram")

g=log2(GSEdata)
abline(v = mean(g,na.rm=TRUE), col = 'red', lwd=3, lty=2)
text(x=mean(g), y=1, "mean", cex=1.2,col = "red")


xfit <- seq(min(g), max(g), length = 50) 
yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
yfit <- yfit * diff(h$mids[1:2]) * length(g) 

lines(xfit, yfit, col = "blue", lwd = 2)
legend("topleft", inset=.05, legend=c("data", "Gaussian"),
       fill=c("gray", "blue"),cex=0.8,text.font=4)

dev.off()

data2=log2(GSEdata)


pdf("plots/boxPlotLog2Data.pdf")
# box-and-whisker plot
#par(mar=c(7,4,2,1))
title <- paste ("GSE156445", "/", "GPL17586", sep ="")
boxplot(data2, col=c("blue", "blue" , "blue", 
                     "red","red","red","red",
                     "red","red","red","red","red","red","red",
                     "red","red","red","red"), main="GSE16515 - boxplots", las=2,cex.axis=0.6)

legend("topleft", inset=.05, legend=c("vehicle treated cells", "cells treated with CIPA"),
       fill=c("blue", "red"),cex=0.8,text.font=4)

#boxplot(data2, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2,cex.axis=0.6)
dev.off()


#legend("center", inset=.02, legend=c("tumor", "normal"),
      # fill=c("red", "blue"),cex=0.8,text.font=4)
#dev.off()

# Hierarchical clustering of the "samples" based on
# the correlation coefficients of the expression values
hc = hclust(as.dist(1-cor(data2)), method="complete")


# add horizontal line to illustrate cutting dendrogram

pdf("plots/HclustLog2Data.pdf")
plot(hc, main="GSE156445 - Hierarchical Clustering",cex=0.6)
abline(h = 0.15, col = "red", lwd = 2) 
dev.off()



vehicleTreated=data2[,1:3]
treatment_CIPAug=data2[,4:18]
#treatment10ug=exdata.mat[,7:9]
#treatment100ug=exdata.mat[,10:12]
#treatment500ug=exdata.mat[,13:15]
#treatment1000ug=exdata.mat[,16:18]

# Compute the means of the samples of each condition
vhc.tr.mean = apply(vehicleTreated, 1, mean)
tr_CIPAug.mean = apply(treatment_CIPAug, 1, mean)
#tr10.mean=apply(treatment10ug, 1, mean)
#tr100.mean=apply(treatment100ug, 1, mean)
#tr500.mean=apply(treatment500ug, 1, mean)
#tr1000.mean=apply(treatment1000ug, 1, mean)

limit=max(vhc.tr.mean,tr_CIPAug.mean)#,tr10.mean,tr100.mean,tr500.mean,tr1000.mean)

# Scatter plot
pdf("plots/Log2raw_vhc.tr vs tr_CIPAug Scatter.pdf")
plot(vhc.tr.mean ~ tr_CIPAug.mean, xlab = "mean vehicle treated", ylab = "mean treated with CIPA",
     main = "vhc.tr and tr_CIPAug", xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "red")
dev.off()




#######################################
# Differential expression (DE) analysis
#######################################

# Separate each of the conditions into three smaller data frames
#tumor = data2[,c(1:4,6,8,10,12:14,16,18,20:23,25:29,31,33,35:38,40:42,44:46,48,50:51)]
#normal = data2[,-c(1:4,6,8,10,12:14,16,18,20:23,25:29,31,33,35:38,40:42,44:46,48,50:51)]

DEGs_byTreatments<-read.csv("DEGs_unique_ex.csv")

genes<-DEGs_byTreatments$X
samples<-colnames(DEGs_byTreatments[,-1])

ex.mat<-data.matrix(DEGs_byTreatments,rownames.force = NA)


exdata.mat=ex.mat[,-1]
#DEdata.mat[,1]<-GeneNames
row.names(exdata.mat)<-genes
str(exdata.mat)
dim(exdata.mat)

exdata.mat2 = log2(exdata.mat)

#rownames(exdata.mat2)

# vehicleTreated=rowMeans(exdata.mat[,1:3])
# treatment1ug=rowMeans(exdata.mat[,4:6])
# treatment10ug=rowMeans(exdata.mat[,7:9])
# treatment100ug=rowMeans(exdata.mat[,10:12])
# treatment500ug=rowMeans(exdata.mat[,13:15])
# treatment1000ug=rowMeans(exdata.mat[,16:18])

vehicleTreated=exdata.mat2[,1:3]
treatment_CIPAug=exdata.mat2[,4:18]


#treatment10ug=exdata.mat[,7:9]
#treatment100ug=exdata.mat[,10:12]
#treatment500ug=exdata.mat[,13:15]
#treatment1000ug=exdata.mat[,16:18]

# Compute the means of the samples of each condition
vhc.tr.mean = apply(vehicleTreated, 1, mean)
tr_CIPAug.mean = apply(treatment_CIPAug, 1, mean)
DEGs_control_treatment=data.frame("genes"=genes,"control"=vhc.tr.mean, "treatment"=tr_CIPAug.mean)
write.csv(DEGs_control_treatment,"DEGs_byControlvsTreatment.csv")
#tr10.mean=apply(treatment10ug, 1, mean)
#tr100.mean=apply(treatment100ug, 1, mean)
#tr500.mean=apply(treatment500ug, 1, mean)
#tr1000.mean=apply(treatment1000ug, 1, mean)

limit=max(vhc.tr.mean,tr_CIPAug.mean)#,tr10.mean,tr100.mean,tr500.mean,tr1000.mean)

# Scatter plot
pdf("plots/DEGS_vhc.tr vs tr_CIPAug Scatter.pdf")
plot(vhc.tr.mean ~ tr_CIPAug.mean, xlab = "mean vehicle treated", ylab = "mean treated with CIPA",
     main = "DEGs_vhc.tr and tr_CIPAug", xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "red")
dev.off()




#####by fold change #########

fc_vhc_trCIPAug = vhc.tr.mean-tr_CIPAug.mean


# Histogram of the fold differences
pdf("plots/histFCvhc_otherTrts.pdf")
#par(mfrow=c(3,2))
hist(fc_vhc_trCIPAug, col = "gray")
dev.off()


#dev.off()
# Compute statistical significance (using t-test)
pvalue_vhc_trCIPAug = NULL # Empty list for the p-values
tstat_vhc_trCIPAug = NULL# Empty list of the t test statistics

for(i in 1 : nrow(exdata.mat2)) { # For each gene : 
  x = vehicleTreated[i,] # control of gene number i
  y = treatment_CIPAug[i,] # treated of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue_vhc_trCIPAug[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat_vhc_trCIPAug[i] = t$statistic
}

#head(pvalue_vhc_tr100)
head(pvalue_vhc_trCIPAug)
pdf("plots/vhc_trCIPAug_hist&scatter.pdf")
par(mfrow=c(2,1))
hist(-log10(pvalue_vhc_trCIPAug), breaks = 100, ylim=c(0,120),col = "blue")
plot(fc_vhc_trCIPAug, -log10(pvalue_vhc_trCIPAug), main = "GSE156445 - Volcano")
fold_cutoff = 1 #0.05
pvalue_cutoff = 0.05
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
dev.off()

# Highlighting up-regulated in red and down-regulated in blue
pdf("plots/up&down.pdf")
plot(fc_vhc_trCIPAug, -log10(pvalue_vhc_trCIPAug), main = "GSE156445 - Volcano #3")
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
points (fc_vhc_trCIPAug[filter_combined & fc_vhc_trCIPAug < 0],
        -log10(pvalue_vhc_trCIPAug[filter_combined & fc_vhc_trCIPAug < 0]),
        pch = 16, col = "red")
points (fc_vhc_trCIPAug[filter_combined & fc_vhc_trCIPAug > 0],
        -log10(pvalue_vhc_trCIPAug[filter_combined & fc_vhc_trCIPAug > 0]),
        pch = 16, col = "blue")
legend("topleft",c("NO","down", "up"),cex=.8,
       col=c("black","red","blue"),pch=c(21,19,19),title="Genes")
dev.off()


filter_by_fold = abs(fc_vhc_trCIPAug) >= fold_cutoff
dim(exdata.mat2[filter_by_fold, ])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue_vhc_trCIPAug <= pvalue_cutoff
dim(exdata.mat2[filter_by_pvalue, ])

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue

filtered_vhc_trCIPAug = exdata.mat2[filter_combined,]
dim(filtered_vhc_trCIPAug)

write.csv(filtered_vhc_trCIPAug,"DEGs_vhc_trCIPAug.csv")

upRegulated_DEGs = exdata.mat2[filter_combined & fc_vhc_trCIPAug > 0,]
dim(upRegulated_DEGs)
write.csv(upRegulated_DEGs,"vhc_trCIPAug_upregulatedDEGs.csv")

downRegulated_DEGs = exdata.mat2[filter_combined & fc_vhc_trCIPAug < 0,]
dim(downRegulated_DEGs)
write.csv(downRegulated_DEGs,"vhc_trCIPAug_downregulatedDEGs.csv")








# Cluster the rows (genes) & columns (samples) by correlation
rowvUP = as.dendrogram(hclust(as.dist(1-cor(t(upRegulated_DEGs)))))
colvUP = as.dendrogram(hclust(as.dist(1-cor(upRegulated_DEGs))))


# Generate a heatmap
pdf("plots/upRegulated_heatmap.pdf")
heatmap.2(upRegulated_DEGs, Rowv=rowvUP,Colv=colvUP, 
          col = rev(greenred(20*10)), cexRow=0.6,
          cexCol=0.6,scale="row",notecol="black",
          margins=c(5,5), # ("margin.Y", "margin.X")
          trace='none', 
          symkey=FALSE, 
          symbreaks=FALSE, 
          dendrogram='none',
          density.info='histogram', 
          denscol="black",
          keysize=1, 
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(3.5,0,3,0)))
dev.off()

#expres<-apply(upRegulated_DEGs,1,mean)
#upRegulated_DEGs=data.frame(rownames(upRegulated_DEGs),expres)
#write.csv(upRegulated_DEGs,"upRegulated_DEGs.csv")

#downRegulated_DEGs=  data2[filter_combined & fold < 0,]
#dim(downRegulated_DEGs)

# Cluster the rows (genes) & columns (samples) by correlation
rowvDown = as.dendrogram(hclust(as.dist(1-cor(t(downRegulated_DEGs)))))
colvDown = as.dendrogram(hclust(as.dist(1-cor(downRegulated_DEGs))))



pdf("plots/downRegulated_heatmap.pdf")
 heatmap.2(downRegulated_DEGs, Rowv=rowvDown,Colv=colvDown, 
          col = rev(greenred(20*10)), cexRow=0.6,
          cexCol=0.6,scale="row",notecol="black",
          margins=c(5,5), # ("margin.Y", "margin.X")
          trace='none', 
          symkey=FALSE, 
          symbreaks=FALSE, 
          dendrogram='none',
          density.info='histogram', 
          denscol="black",
          keysize=1,
          key.par=list(mar=c(3.5,0,3,0)))
dev.off()

expres<-apply(downRegulated_DEGs,1,mean)
downRegulated_DEGs=data.frame(rownames(downRegulated_DEGs),expres)
write.csv(downRegulated_DEGs,"downRegulated_DEGs.csv")
# Cluster the rows (genes) & columns (samples) by correlation
#rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
#colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))




# install.packages("gplots")		# Uncomment if not already installed
# install.packages("RColorBrewer")	# Uncomment if not already installed

library(gplots)

# Enhanced heatmap
#pdf("plots/heatmap2.pdf")
#heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
         # col = rev(redblue(256)), scale = "row")
#dev.off()
# Save the heatmap to a PDF file

# Cluster the rows (genes) & columns (samples) by correlation
rowvDEG = as.dendrogram(hclust(as.dist(1-cor(t(filtered_vhc_trCIPAug)))))
colvDEG = as.dendrogram(hclust(as.dist(1-cor(filtered_vhc_trCIPAug))))


pdf ("plots/GSE156445_DE_Heatmap_Genes.pdf")
#pdf("plots/upRegulated_heatmap.pdf")
heatmap.2(filtered_vhc_trCIPAug, Rowv=rowvDEG,Colv=colvDEG, 
          col = rev(greenred(20*10)), cexRow=0.6,
          cexCol=0.6,scale="row",notecol="black",
          # ("margin.Y", "margin.X")
          trace='none', 
          symkey=FALSE, 
          symbreaks=FALSE, 
          dendrogram ='both',
          density.info='histogram', 
          denscol="black",
          keysize=1,
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(3,0,3,0)))
dev.off()



#############################################################################

