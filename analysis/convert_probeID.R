# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

DEGs<-read.csv("DEGslist.csv")

filtered=gsub(".1","", DEGs$V1)
filtered<-as.list(filtered)
filtered=setDT(transpose(filtered, fill=0))[]
write.csv(filtered,"probeIds_dot.csv")


#str(gset)
#ex_cleaned<-read.csv("cleaned_ex.csv")
#require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)


genes <- getBM(attributes = c("affy_hta_2_0",#"ensembl_gene_id", 
                              "external_gene_name"),
               filters = "affy_hta_2_0",
               values = filtered$V1,
               mart = mart,useCache =  FALSE)

saveRDS(genes,"DEGs.RDS")




