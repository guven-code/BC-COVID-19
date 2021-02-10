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



vehicleTreated=exdata.mat[,1:3]
treatment1ug=exdata.mat[,4:6]
treatment10ug=exdata.mat[,7:9]
treatment100ug=exdata.mat[,10:12]
treatment500ug=exdata.mat[,13:15]
treatment1000ug=exdata.mat[,16:18]




tr.df.All= data.frame(vhc_tr=vehicleTreated,tr1=treatment1ug,tr10=treatment10ug,
                      tr100=treatment100ug, tr500=treatment500ug, tr1000=treatment1000ug)

library(data.table)

# get data


# transpose
t_tr.df.All <- t(tr.df.All)



# get row and colnames in order
colnames(t_tr.df.All) <- probe_ids
rownames(t_tr.df.All) <- samples
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

#merge1$treatment <- as.factor(merge1$treatment, levels=c("vhc_tr","tr1","tr10","tr100", "tr500", "tr1000"))


p <- vector(mode = "list", length(2:ncol(merge1)))
anova_table<-vector(mode = "list", length(2:ncol(merge1)))
pr<-vector(mode = "list", length(2:ncol(merge1)))
pr_significant<-vector(mode = "list", length(2:ncol(merge1)))
#colnames(merge1) <- probe_ids
#rownames(merge) <- samples

baseformula <- " ~ treatment"
for (i in 2:(ncol(merge1))) {
  formula <- paste(colnames(merge1)[i], baseformula, sep="")
  
  p[[i]] <- summary(aov(as.formula(formula), data=merge1))[[1]][["Pr(>F)"]][1]
  #p<0.05
  anova_table[[i]] <- as.data.frame(summary(aov(as.formula(formula), data=merge1))
                                    [[1]][["Pr(>F)"]][1])
  #p<0.05 
  
  pr[[i]]=paste(formula, ": p=", p[[i]], sep="")
  pr_significant[[i]]=paste(formula, ": p=", p[[i]]<0.05, sep="")
  #write.csv(print[[i]])
  
}



#capture.output(p, file = "p_value.txt")
capture.output(anova_table, file = "anova results.txt")
capture.output(pr, file = "pvalue_Genes.txt")
capture.output(pr_significant, file = "psignificant_Genes.csv")

significantGenes=read.csv("psignificant_Genes.csv")


lst<-t(pr_significant)
DEGsll<-lst[grep("TRUE", lst)]

DEGs=setDT(transpose(DEGsll, fill=0))[]

#write.csv(DEGs,"DEGslist.csv")

filteredProbeIDs=read.csv("DEGslist.csv")
filtered=gsub("~","", filteredProbeIDs$probe_ids)
filtered<-as.list(filtered)
filtered=setDT(transpose(filtered, fill=0))[]

write.csv(filtered,"DEGslist.csv")






