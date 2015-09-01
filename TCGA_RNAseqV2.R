library(TCGA-Assembler)
library(reshape2)

# install TCGA-Assembler and use Module_A & Module_B
setwd("TCGA-Assembler")
source("./Module_A.r");
source("./Module_B.r");

# read in TCGA annotation file
# prepare this file using prepTCGA_Annotation.R
TCGA.annot <- read.delim('TCGA.Sept2010.09202010.gene.bed',header=T,check.names=F) 

# use biomaRt or biomart to get all lincRNAs in hg19 Ensembl75
# get EntrezGene.ID, Associated.Gene.Name and Gene.Biotype
# read in all hg19 ensembl75 lincRNA annotation
hg19.ensembl.lincs <- read.delim('all_hg19_ensembl75_lincRNAs.txt',header=T)

# get RNASeq V2 data
# set cancerType parameter
# 8 cancer types but can be expanded to as many as you want
# this can be converted into a function - to do
RNASeqV2RawData.luad <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./QuickStartGuide_Results/luad_RawData/", cancerType = "LUAD",assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results",tissueType = c("TP","NT"), outputFileName ="2minuteExample");
RNASeqV2RawData.lusc <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./QuickStartGuide_Results/luad_RawData/", cancerType = "LUSC",assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results",tissueType = c("TP","NT"), outputFileName ="2minuteExample");
RNASeqV2RawData.lihc <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./QuickStartGuide_Results/luad_RawData/", cancerType = "LIHC",assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results",tissueType = c("TP","NT"), outputFileName ="2minuteExample");
RNASeqV2RawData.hnsc <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./QuickStartGuide_Results/luad_RawData/", cancerType = "HNSC",assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results",tissueType = c("TP","NT"), outputFileName ="2minuteExample");
RNASeqV2RawData.coad <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./QuickStartGuide_Results/luad_RawData/", cancerType = "COAD",assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results",tissueType = c("TP","NT"), outputFileName ="2minuteExample");
RNASeqV2RawData.paad <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./QuickStartGuide_Results/luad_RawData/", cancerType = "PAAD",assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results",tissueType = c("TP","NT"), outputFileName ="2minuteExample");
RNASeqV2RawData.stad <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./QuickStartGuide_Results/luad_RawData/", cancerType = "STAD",assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results",tissueType = c("TP","NT"), outputFileName ="2minuteExample");
RNASeqV2RawData.thca <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./QuickStartGuide_Results/luad_RawData/", cancerType = "THCA",assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results",tissueType = c("TP","NT"), outputFileName ="2minuteExample");

#################################################### processRNASeqV2Data ##################################################### 
# function to process RNASeqV2 data from TCGA
processRNASeqV2Data <- function(TCGA.annot,hg19.ensembl.lincs,RNASeqV2RawData)
{
  rnaseqdata <- unlist(RNASeqV2RawData[[1]])
  colnames(rnaseqdata) <- sub('A-[0-9]{2}R-[A-Z0-9]{4}-07','',RNASeqV2RawData[[1]][1,])
  rnaseqdata <- rnaseqdata[-c(1:2),]
  colnames(rnaseqdata)[1] <- "gene_id"
  
  # separate into cases & controls
  cases <- rnaseqdata[,c(1,grep("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01",colnames(rnaseqdata)))]
  controls <- rnaseqdata[,c(1,grep("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-11",colnames(rnaseqdata)))]
  
  colnames(cases) <- sub("-01","",colnames(cases))
  colnames(controls) <- sub("-11","",colnames(controls))
  
  # use only matched samples
  cases <- cases[,which(!is.na(match(colnames(cases),colnames(controls))))]
  controls <- controls[,which(!is.na(match(colnames(controls),colnames(cases))))]
  
  cases <- as.data.frame(cases)
  cases[,2:ncol(cases)] <- sapply(cases[,2:ncol(cases)],function(x) as.numeric(levels(x)[x]))
  
  controls <- as.data.frame(controls)
  controls[,2:ncol(controls)] <- sapply(controls[,2:ncol(controls)],function(x) as.numeric(levels(x)[x]))
  
  # arrange colnames in cases & controls to match 
  cases <- cases[,order(colnames(cases))]
  controls <- controls[,order(colnames(controls))]
  
  # calculate FC
  fc <- data.frame(id=cases$gene_id, fc = rowMeans(log2(cases[,2:ncol(cases)]+1))-rowMeans(log2(controls[,2:ncol(controls)]+1)))
  
  # calculate paired t-test
  ttest.df <- data.frame(id = cases$gene_id)
  ttest.df <- cbind(ttest.df, t(sapply(1:nrow(cases), function(i) t.test(x=as.numeric(cases[i,2:ncol(cases)]),y=as.numeric(controls[i,2:ncol(controls)]),paired=T,na.action=na.omit)[c(3,1)])))
  
  # merge t-test results & fc
  tcga.lincs <- merge(ttest.df,fc,by="id")
  
  # split entrez id & gene names
  tcga.lincs <- cbind(tcga.lincs,colsplit(tcga.lincs$id,pattern="[|]",names=c("genename","entrez")))
  
  # check which are lincs from allhg19lincs
  tcga.lincs <- merge(tcga.lincs,hg19.ensembl.lincs,by.x="entrez",by.y="EntrezGene.ID")
  
  # p-adjust
  tcga.lincs$padj <- p.adjust(p=tcga.lincs$p.value,method="fdr")
  
  # p-adjust filter = 0.05
  tcga.lincs.sig.005 <- tcga.lincs[which(tcga.lincs$padj < 0.05),]
  if(nrow(tcga.lincs.sig.005)==0)
  {
    stop("No significantly differentially expressed lincRNAs")
  }
  
  tcga.lincs.sig.005.coord <- TCGA.annot[which(TCGA.annot$name %in% tcga.lincs.sig.005$id),] 

  # melt the dataframe of cases & controls
  dat.melted <- rbind(data.frame(melt(controls),treat="Normal"),data.frame(melt(cases),treat="Cancer"))
  dat.melted <- merge(dat.melted,tcga.lincs.sig.005,by.x="gene_id",by.y="id")
  dat.melted$label <- with(dat.melted,paste("Gene:",Associated.Gene.Name,"\n","FoldChange:",round(fc,3),"\n","Pvalue.Adj:",eval(padj)))
  
  resultList <- list("dat" = dat.melted, "tcga.lincs.sig.005.coord" = tcga.lincs.sig.005.coord)
  
  return(resultList)
}

# run the above function for different cancer type
# here we run it for LIHC
resultList <- processRNASeqV2Data(TCGA.annot,hg19.ensembl.lincs,RNASeqV2RawData.lihc)

# lincRNAs that are significantly differentially expressed at pvalue cutoff of 0.05
# t-test results for the lincRNAs
t.test.res <- resultList$dat
# coordinates in bed file format
tcga.lincs.sig.005.coord <- resultList$tcga.lincs.sig.005.coord
