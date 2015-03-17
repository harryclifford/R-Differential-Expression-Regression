#########################################
# Script to run DESeq2 time series test #
#########################################


### user input

# load DESeq2
library("DESeq2")
library("gplots")
library("RColorBrewer")


# set working directory
setwd("/net/isi-scratch/")

# input counts file
# (including header, first column is gene ids, the rest are each sample's counts)
# (each sample should be named by sample, then an underscore, then by replicate)
# (e.g. P0_R1, P0_R2, P1_R1, P1_R2)
# (because the final plot looks at means of replicates, and will merge by pre-underscore name)
incountfile <- "all_counts.txt"

# sample information data frame
# (column named "condition", rownames match colnames in above counts file)
colData <- data.frame( condition=sort(rep(c(0,4,8,14,21),3)) )
rownames(colData) <- c(
                        "P0_R1","P0_R2","P0_R3",
                        "P4_R1","P4_R2","P4_R3",
                        "P8_R1","P8_R2","P8_R3",
                        "P14_R1","P14_R2","P14_R3",
                        "P21_R1","P21_R2","P21_R3"
                        )

# give FDR threshold for significance
FDR <- 0.05

## for subset heatmap
# if you would like to convert Ensembl gene names, provide biomart file location
# (ensembl biomart, csv, with header, 1st col = ensembl gene name, 2nd col = assoc. gene name)
# (if you would not like this, state biomart_file = NA)
biomart_file <- "mart_export.txt"
# give number of genes (will be doubled, as done for up and down) to include
# (suggested value of 50 - totalling 100 genes)
heatmap_number <- 50


########################
########################
########################


### main script

# reads in counts
# (only takes those included in the colData, and in the order of the colData)
countData <- read.table(incountfile,header=T,row.names=1,stringsAsFactors=F)[,rownames(colData)]

# creates DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
#dds$condition <- factor(dds$condition)

# filter out zero counts
idx <- which(rowSums(counts(dds)) > 0)
dds <- dds[idx,]

# applies DESeq, with likelihood ratio test
dds2 <- DESeq(dds, test="LRT", full= ~ condition, reduced= ~ 1)

# results table
res <- results(dds2)
resOrdered <- res[order(res$padj),]
sig_res <- res[which(res$padj<FDR),]


### Commonly used plots / results / a few extras, written to results folder

dir.create("results",showWarnings=F)

# writes size factors and plots
write.table(sizeFactors(dds2),file="results/size_factors.txt",quote=F,sep="\t",col.names=F)
png(file="results/size_factors.png",height=650,width=650)
barplot(sizeFactors(dds2),main="Size Factors",las=2)
dev.off()

# writes different modifications of counts
write.table(counts(dds2,normalized=T),file="results/normalized_counts.txt",quote=F,sep="\t")
vsd <- assay(DESeq2::varianceStabilizingTransformation(dds2,blind=T))
write.table(vsd,file="results/variance_stabilized_counts.txt",quote=F,sep="\t")

# heatmap of all genes
png(file="results/heatmap_distances_allgenes.png",height=650,width=650)
heatmap.2(as.matrix(dist(t(vsd))), trace="none")
dev.off()

# basic PCA plot
prcomp_output<-prcomp(t(vsd))
png(file="results/PCA_allgenes.png",height=650,width=650)
plot(prcomp_output$x[,1:2],pch=20,
    main="Variance Stabilized Principal Component Analysis",
    xlab=paste("PC1(",round(summary(prcomp_output)$importance[2,1]*100),"% of Variance)",sep=""),
    ylab=paste("PC2(",round(summary(prcomp_output)$importance[2,2]*100),"% of Variance)",sep="")
    )
text(prcomp_output$x[,1:2], rownames(prcomp_output$x[,1:2]), pos=3, cex=0.9)
dev.off()

# main results
write.table(res,file="results/main_results.txt",quote=F,sep="\t")
write.table(sig_res,file="results/significant_results.txt",quote=F,sep="\t",row.names=F)

# results description
write.table(mcols(res)$description,file="results/main_results_description.txt",quote=F,sep="\n",row.names=F,col.names=F)
### NOTE: ensure these state "condition" as a factor
### as opposed to something like "condition 21 vs 0"

# p-value plot
png(file="results/pvalplot.png",height=650,width=650)
hist(res$pvalue,breaks=100,col="skyblue",border="slateblue",main="P-value distribution",xlab="P-values")
dev.off()

# dispersion estimates plot
png("results/dispersion_estimates.png")
plotDispEsts(dds2)
dev.off()



### Plots of significantly differentially expressed genes


## Subset heatmap

# subsets to top differentially expressed genes
vsd_up <- vsd[rownames(resOrdered[which(resOrdered$log2FoldChange > 0)[1:heatmap_number],]),]
vsd_down <- vsd[rownames(resOrdered[which(resOrdered$log2FoldChange < 0)[1:heatmap_number],]),]
vsd_sub <- rbind(
            vsd_up[order(vsd_up[,ncol(vsd_up)],decreasing=T),] ,
            vsd_down[order(vsd_down[,ncol(vsd_down)],decreasing=T),]
            )

# determines gene names
if(!is.na(biomart_file)){
    biomart <- read.csv(biomart_file,stringsAsFactors=F)
    ylabels <- biomart[match(rownames(vsd_sub),biomart[,1]),2]
}else{ylabels<-rownames(vsd_sub)}

hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)

# plots heatmap
png(file="results/heatmap_distances_topsubset.png",height=2000,width=1000)
heatmap.2( vsd_sub,Rowv=F,Colv=F,scale="none",dendrogram="none",trace="none",key=F,col=hmcol,
        rowsep=heatmap_number,sepwidth=c(1,2),margin=c(15,20),lhei=c(1,10),cexCol=2,cexRow=2,
        labRow=c(ylabels[1:heatmap_number],"","",ylabels[(heatmap_number+1):(heatmap_number*2)])
        )
title(main=paste("Top",heatmap_number*2,"Differentially\nExpressed Genes"),line=-2,cex.main=3)
mtext("Samples",side=1,line=3,cex=3)
mtext(paste("Top",heatmap_number,"\nup-\nregulated\nGenes"),side=2,line=-10,las=1,at=0.75,cex=3)
mtext(paste("Top",heatmap_number,"\ndown-\nregulated\nGenes"),
        side=2,line=-10,las=1,at=0.25,cex=3)
dev.off()


## Expression level plot

# splits data into up and downreg
upreg_sig_res <- sig_res[which(sig_res$log2FoldChange > 0),]
upreg_sig <- counts(dds2,normalized=T)[rownames(upreg_sig_res),]
downreg_sig_res <- sig_res[which(sig_res$log2FoldChange < 0),]
downreg_sig <- counts(dds2,normalized=T)[rownames(downreg_sig_res),]
# uses normalized data, can also use variance stabilized:
#upreg_sig <- vsd[rownames(upreg_sig_res),]
#downreg_sig <- vsd[rownames(downreg_sig_res),]

# averages samples using first column of colData
upreg_means <- c()
downreg_means <- c()
for(i in unique(colData[,1])){
    
    samples <- rownames(colData)[which(colData[,1]==i)]
    meanname <- paste(unique(unlist(lapply( strsplit(samples,"_"), `[[`, 1))),collapse="+")
    
    upreg_sample_means <- apply(upreg_sig[,samples],1,mean)
    upreg_means <- cbind(upreg_means,upreg_sample_means)
    colnames(upreg_means)[ncol(upreg_means)] <- meanname
    
    downreg_sample_means <- apply(downreg_sig[,samples],1,mean)
    downreg_means <- cbind(downreg_means,downreg_sample_means)
    colnames(downreg_means)[ncol(downreg_means)] <- meanname
}

# standardizes data for directly comparable plots
upreg_std <- t(apply(upreg_means,1,function(p) (p-mean(p))/sd(p)))
downreg_std <- t(apply(downreg_means,1,function(q) (q-mean(q))/sd(q)))

# sorts data by decreasing p-value (so heat colour increases as p-value decreases)
upreg <- upreg_std[ rownames(upreg_sig_res)[order(upreg_sig_res$padj,decreasing=T)] ,]
downreg <- downreg_std[ rownames(downreg_sig_res)[order(downreg_sig_res$padj,decreasing=T)] ,]

# plots data
ylims <- c( floor(min(c(upreg,downreg))) , ceiling(max(c(upreg,downreg))) )
upreg_cols <- colorRampPalette( brewer.pal(7,"YlOrRd") )( nrow(upreg) )
downreg_cols <- colorRampPalette( brewer.pal(7,"YlGnBu") )( nrow(downreg) )

png(file="results/DEgenes_expression.png",width=700)

par(mfrow=c(1,2)) 
plot( upreg[1,],type="l",ylim=ylims, col=upreg_cols[1],
        main="Upregulated genes\ncoloured by p-value",
        xaxt="n", xlab="Developmental Timepoint",
        ylab="Standardized Mean Normalized Expression Level"
    )
axis(1,at=axTicks(1),labels=colnames(upreg))
for(j in 2:nrow(upreg)){ lines( upreg[j,],type="l", col=upreg_cols[j] )}

plot( downreg[1,],type="l",ylim=ylims, col=downreg_cols[1],
        main="Downregulated genes\ncoloured by p-value",
        xaxt="n", xlab="Developmental Timepoint",
        ylab="Standardized Mean Normalized Expression Level"
    )
axis(1,at=axTicks(1),labels=colnames(upreg))
for(k in 2:nrow(downreg)){ lines( downreg[k,],type="l", col=downreg_cols[k] )}

dev.off()




