# Same kind of analysis as TCGA_RNASeqV2.R but for RNAseq values 

# example normalized counts downloaded & renamed for ESCA in folder TCGA/ESCA/genecounts
setwd('TCGA/ESCA/')
cases    <- list.files(path="genecounts",pattern="-01",full.names=T)
controls <- list.files(path="genecounts",pattern="-11",full.names=T)

# get normalized values for cases
dat.cases <- data.frame(id=read.delim(cases[1])$gene)
for(i in 1:length(cases))
{
  tmp <- read.delim(cases[i])
  dat.cases <- cbind(dat.cases,x=tmp$median_length_normalized)
  colnames(dat.cases)[i+1] = sub("genecounts/","",cases[i])
}

# get normalized values for controls
dat.control <- data.frame(id=read.delim(controls[1])$gene)
for(i in 1:length(controls))
{
  tmp <- read.delim(controls[i])
  dat.control <- cbind(dat.control,x=tmp$median_length_normalized)
  colnames(dat.control)[i+1] = sub("genecounts/","",controls[i])
}

# fc = data.frame(id=dat.cases$id,fc=rowMeans(dat.cases[,2:ncol(dat.cases)])/rowMeans(dat.control[,2:ncol(dat.control)]))
# fc = fc[-which(fc$fc == "Inf" | is.na(fc$fc)),]
# fc$logfc = log2(fc$fc)

# calculate fold change
fc <- data.frame(id=dat.cases$id, FoldChange = rowMeans(log2(dat.cases[,2:ncol(dat.cases)]+1))-rowMeans(log2(dat.control[,2:ncol(dat.control)]+1)))

# get pvalue & t-statistic
ttest.df <- data.frame(id=dat.cases$id)
ttest.df <- cbind(ttest.df, t(sapply(1:nrow(dat.cases), function(i) t.test(x=as.numeric(dat.cases[i,2:ncol(dat.cases)]),
                                                                           y=as.numeric(dat.control[i,2:ncol(dat.control)]),
                                                                           paired=T,na.action=na.omit)[c(3,1)])))
# merge t-test results with fold change
tcga.lincs <- merge(ttest.df,fc,by="id")
tcga.lincs$id <- sub("_calculated","",tcga.lincs$id) # remove unwanted strings
tcga.lincs <- cbind(tcga.lincs,colsplit(tcga.lincs$id,pattern="[|]",names=c("genename","entrez"))) # split id into symbol & entrezid
tcga.lincs <- merge(tcga.lincs,allhg19lincs,by.x="entrez",by.y="EntrezGene.ID") # get all lncRNAs
tcga.lincs$Pvalue.Adj <- p.adjust(p=tcga.lincs$p.value,method="fdr") # p.adjust only lncRNAs
tcga.lincs.sig.005 <- tcga.lincs[which(tcga.lincs$padj < 0.05),] # significant pvalue

# all the coordinates for TCGA gene data
final <- read.delim('TCGA.Sept2010.09202010.gene.bed',header=T)

# write out coordinates for significant lncRNAs
write.table(x=final[which(final$name %in% tcga.lincs.sig.005$id),],file='hg19_siggenes.bed',quote=F,row.names=F,col.names=F,sep="\t")

#########################################################################################################
# if you want to make plots for expression
# make plot for LOC647979
geneset <- rbind(data.frame(melt(dat.control),treat="Normal"),data.frame(melt(dat.cases),treat="Cancer"))
geneset <- cbind(dat.total.m,colsplit(dat.total.m$id,"[|]",c("gene","entrez")))
geneset <- dat.total.m[grep("LOC647979",dat.total.m$gene),] 
geneset$entrez <- sub('_calculated','',geneset$entrez)
geneset <- merge(geneset,tcga.lincs,by="entrez")

# geneset$FoldChange <- ifelse(geneset$logfc<0,round(1/geneset$FoldChange*(-1),3),round(1/geneset$FoldChange,3))
geneset$label <- with(geneset,paste("Gene:",Associated.Gene.Name,"\n","FoldChange:",FoldChange,"\n","Pvalue.Adj:",Pvalue.Adj))

ggplot(geneset,aes(x=treat,y=log2(value+1),fill=treat)) + geom_boxplot(show_guide=F) + 
  facet_wrap(~label,ncol=4) + ylab("log2(Expression)\n") + 
  ggtitle("TCGA: Boxplot of Differentially Expressed lincRNAs\n") +
  theme(axis.text.x=element_text(size=14,color="black"),
        axis.text.y=element_text(size=14,color="black"),
        strip.text=element_text(size=14,color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14,color="black"),
        plot.title=element_text(size=16,color="black"))
#########################################################################################################
