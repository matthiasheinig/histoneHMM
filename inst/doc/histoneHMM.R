### R code from vignette source 'histoneHMM.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: foo
###################################################
options(keep.source = TRUE, width = 60)
foo <- packageDescription("histoneHMM")


###################################################
### code chunk number 3: histoneHMM.Rnw:36-37
###################################################
library(histoneHMM)


###################################################
### code chunk number 4: histoneHMM.Rnw:41-42
###################################################
data(rat.H3K27me3)


###################################################
### code chunk number 5: histoneHMM.Rnw:46-48
###################################################
posterior = run.univariate.hmm("BN_H3K27me3.txt", data=BN, n.expr.bins=5, 
    em=TRUE, chrom="chr19", maxq=1-1e-3, redo=F)


###################################################
### code chunk number 6: histoneHMM.Rnw:53-55
###################################################
n = load("BN_H3K27me3-zinba-params-em.RData")
zinba.mix = get(n)


###################################################
### code chunk number 7: histoneHMM.Rnw:59-62
###################################################
max = 100
x = BN[BN[,"chrom"] == "chr19","signal"]
x[x > max] = NA


###################################################
### code chunk number 8: histoneHMM.Rnw:66-68
###################################################
plotDensity(model=zinba.mix, x=x, xlim=c(0, max), add=FALSE, 
            col="black", alpha=1)


###################################################
### code chunk number 9: histoneHMM.Rnw:72-87
###################################################
posterior = run.univariate.hmm("SHR_H3K27me3.txt", data=SHR, n.expr.bins=5, 
    em=TRUE, chrom="chr19", maxq=1-1e-3, redo=F)

## load the fit from the file
n = load("BN_H3K27me3-zinba-params-em.RData")
zinba.mix = get(n)

## for the visualization we ignore the long tail of the distribution
max = 100
x = SHR[SHR[,"chrom"] == "chr19","signal"]
x[x > max] = NA

## visualize the fit
plotDensity(model=zinba.mix, x=x, xlim=c(0, max), add=FALSE, 
            col="black", alpha=1)


###################################################
### code chunk number 10: histoneHMM.Rnw:92-103
###################################################
bivariate.posterior = run.bivariate.hmm("BN_H3K27me3.txt", "SHR_H3K27me3.txt", 
    outdir="BN_vs_SHR_H3K27me3/", data1=BN, data2=SHR, 
    sample1="BN", sample2="SHR", n.expr.bins=5, maxq=1-1e-3, 
    em=TRUE, chrom="chr19")

## call regions
BN.regions = callRegions(bivariate.posterior, 0.5, "BN", NULL)
# GR2gff(BN.regions, "BN_vs_SHR_H3K27me3/BN_specific.gff")

SHR.regions = callRegions(bivariate.posterior, 0.5, "SHR", NULL)
# GR2gff(SHR.regions, "BN_vs_SHR_H3K27me3/SHR_specific.gff")


###################################################
### code chunk number 11: histoneHMM.Rnw:111-116
###################################################
system("wget http://histonehmm.molgen.mpg.de/data/expression.txt")
system("wget http://histonehmm.molgen.mpg.de/data/BN.bam")
system("wget http://histonehmm.molgen.mpg.de/data/BN.bam.bai")
system("wget http://histonehmm.molgen.mpg.de/data/chroms.txt")
system("wget http://histonehmm.molgen.mpg.de/data/ensembl59-genes.gff")


###################################################
### code chunk number 12: histoneHMM.Rnw:120-122
###################################################
expr = read.csv("expression.txt", sep="\t")
mean.expr = apply(expr[,1:5], 1, mean)


###################################################
### code chunk number 13: histoneHMM.Rnw:127-131
###################################################
chroms = read.table("chroms.txt", stringsAsFactors=FALSE)
sl = chroms[,2]
names(sl) = chroms[,1]
genome = GRanges(seqlengths=sl)


###################################################
### code chunk number 14: histoneHMM.Rnw:136-140
###################################################
genes = gff2GR("ensembl59-genes.gff", "ID")
data = preprocess.for.hmm("BN.bam", genes, bin.size=1000, genome, 
    expr=mean.expr, chr=c("chr19", "chr20"))
write.table(data, "BN.txt", sep="\t", quote=F, row.names=F)


