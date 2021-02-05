library("data.table")
library("readr")

#ex<-readRDS("ex.RDS") 

ex<-read.csv("GSE156445raw.csv")

ex<-ex[-grep("\\_", ex$X),]
ex<-ex[grep("\\.", ex$X),]
# write.csv(ex,"regex.csv")

probe_ids<-ex$X
samples<-colnames(ex[,-1])

ex.mat<-data.matrix(ex,rownames.force = NA)

exdata.mat=ex.mat[,-1]
#DEdata.mat[,1]<-GeneNames
row.names(exdata.mat)<-probe_ids
str(exdata.mat)
dim(exdata.mat)

wo.DOT=gsub(".1","", rownames(exdata.mat))
wo.DOT<-as.list(wo.DOT)
wo.DOT=setDT(transpose(wo.DOT, fill=0))[]

rownames(exdata.mat)<-unlist(wo.DOT)

DEGs.names<-readRDS("DEGs.RDS")

DEGs.ex<-exdata.mat[DEGs.names$affy_hta_2_0,]
#dim(DEGs.ex)

rownames(DEGs.ex)<-DEGs.names$external_gene_name



DEGs.unique.ex<-DEGs.ex[!duplicated(rownames(DEGs.ex)),]

write.csv(DEGs.unique.ex,"DEGs_unique_ex.csv")

exdata.mat2<-DEGs.unique.ex

rownames(exdata.mat2)

vehicleTreated=exdata.mat2[,1:3]
treatment1ug=exdata.mat2[,4:6]
treatment10ug=exdata.mat2[,7:9]
treatment100ug=exdata.mat2[,10:12]
treatment500ug=exdata.mat2[,13:15]
treatment1000ug=exdata.mat2[,16:18]




tr.df.All= data.frame(vhc_tr=vehicleTreated,tr1=treatment1ug,tr10=treatment10ug,
                      tr100=treatment100ug, tr500=treatment500ug, tr1000=treatment1000ug)

library(data.table)

# get data


 # transpose
 t_tr.df.All <- t(tr.df.All)
 # get row and colnames in order
 colnames(t_tr.df.All) <- rownames(exdata.mat2)
 rownames(t_tr.df.All) <- colnames(exdata.mat2)

vhc_tr <- data.frame(rep("vhc_tr", nrow(t_tr.df.All)/6), t_tr.df.All[1:3,])
colnames(vhc_tr)[1] <- "treatment"
rownames(vhc_tr)=c(1,2,3)

tr1 <- data.frame(rep("tr1", nrow(t_tr.df.All)/6), t_tr.df.All[4:6,])
colnames(tr1)[1] <- "treatment"
rownames(tr1)=c(1,2,3)

tr10 <- data.frame(rep("tr10", nrow(t_tr.df.All)/6), t_tr.df.All[7:9,])
colnames(tr10)[1] <- "treatment"
rownames(tr10)=c(1,2,3)

tr100 <- data.frame(rep("tr100", nrow(t_tr.df.All)/6), t_tr.df.All[10:12,])
colnames(tr100)[1] <- "treatment"
rownames(tr100)=c(1,2,3)

tr500 <- data.frame(rep("tr500", nrow(t_tr.df.All)/6), t_tr.df.All[13:15,])
colnames(tr500)[1] <- "treatment"
rownames(tr500)=c(1,2,3)

tr1000 <- data.frame(rep("tr1000", nrow(t_tr.df.All)/6), t_tr.df.All[16:18,])
colnames(tr1000)[1] <- "treatment"
rownames(tr1000)=c(1,2,3)

merge1 <- rbind(vhc_tr, tr1, tr10, tr100, tr500, tr1000)


p <- vector(mode = "list", length(2:ncol(merge1)))
anova_table<-vector(mode = "list", length(2:ncol(merge1)))
genes<-vector(mode = "list", length(2:ncol(merge1)))
genes_significant<-vector(mode = "list", length(2:ncol(merge1)))
tukeyHSD<-vector(mode = "list", length(2:ncol(merge1)))
#colnames(merge1) <- probe_ids
#rownames(merge) <- samples

baseformula <- " ~ treatment"
for (i in 2:(ncol(merge1))) {
  formula <- paste(colnames(merge1)[i], baseformula, sep="")
  
  p[[i]] <- summary(aov(as.formula(formula), data=merge1))[[1]][["Pr(>F)"]][1]
  #p<0.05
  anova_table[[i]] <- as.data.frame(summary(aov(as.formula(formula), data=merge1))
                                    [[1]][["Pr(>F)"]][1])
  tukeyHSD[[i]]<-TukeyHSD(aov(as.formula(formula), data=merge1))
  
  genes[[i]]=paste(formula, ": p=", p[[i]], sep="")
  genes_significant[[i]]=paste(formula, ": p=", p[[i]]<0.05, sep="")
  #write.csv(print[[i]])
  
}




capture.output(anova_table, file = "anova_results2.txt")

capture.output(tukeyHSD, file = "tukey_results.txt")
capture.output(genes, file = "pvalue_Genes2.txt")
capture.output(genes_significant, file = "significant_Genes2.csv")

significantGenes=read.csv("significant_Genes2.csv")




lst<-t(genes_significant)
DEGsll<-lst[grep("TRUE", lst)]

DEGs=setDT(transpose(DEGsll, fill=0))[]

write.csv(DEGs,"2ndDEGslist.csv")

filteredGenes=read.csv("2ndDEGslist.csv")

filtered=gsub("~","", filteredGenes$V1)
filtered<-as.list(filtered)
filtered=setDT(transpose(filtered, fill=0))[]

write.csv(filtered,"2ndDEGslist.csv")


for (i in (2:length(tukeyHSD))){
  
 # pdf(paste("plotsL/",i, "_regressionFit.pdf", sep=''))
  pdf(paste("plotsTUKEY/",colnames(merge1)[i], "_tukey.pdf", sep=''))
  plot(tukeyHSD[[i]], las=1,
    col="blue",cex.axis=0.6)
  dev.off()
}





