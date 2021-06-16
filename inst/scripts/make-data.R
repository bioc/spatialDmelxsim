# run with R 4.0.4 and R 3.12
# June 12 2021

### get metadata from GEO

library(GEOquery)
geo <- getGEO("GSE102233")
titles <- c(pData(geo[[2]])$title,pData(geo[[3]])$title)
rearrange <- c(1:27, 53:105, 28:52, 106:132) # for later
slice <- as.numeric(sub(".*sl(.*?)$","\\1",titles))
table(slice)
strain <- substr(titles, 1, 7)
table(strain)
rep <- as.numeric(sub(".*rep(.*?)_.*","\\1",titles))
rep <- ifelse(strain == "simXmel", rep + 3, rep)
table(rep, strain)
table(rep, slice)
length(titles)

# normalized slice, bc rep2 has 26 slices and rep4 has 25
normSlice <- slice
normSlice[rep == 2] <- slice[rep == 2] * 27/26
normSlice[rep == 4] <- slice[rep == 4] * 27/25
#plot(slice, normSlice, col=rep); abline(h=27, v=27, col="grey")

### download and read the counts

## for (i in 2:3) {
##   for (j in seq_along(pData(geo[[i]])$title)) {
##     download.file(pData(geo[[i]])$supplementary_file_1[j],
##                   destfile=paste0("ase_files/",i,"_",j,".txt.gz"), method="wget")
##   }
## }

n <- length(titles)
g <- 13620 # number of genes per ase file
a1 <- matrix(nrow=g, ncol=n) # effect / alt allele (D sim)
a2 <- matrix(nrow=g, ncol=n) # non-effect / ref allele (D mel)
i <- 1
j <- 1
x <- read.delim(paste0("ase_files/",i+1,"_",j,".txt.gz"))
rownames(a1) <- rownames(a2) <- x$FEATURE
for (i in 1:2) {
  for (j in 1:(c(52,80)[i])) {
    cat(j)
    x <- read.delim(paste0("ase_files/",i+1,"_",j,".txt.gz"))
    jj <- j + (i-1)*52
    a1[,jj] <- x$ALT_COUNTS
    a2[,jj] <- x$REFERENCE_COUNTS
  }
}

# fix CG4500 / hll
rownames(a1)[rownames(a1) == "FBgn0028519"] <- "FBgn0286723"
rownames(a2)[rownames(a2) == "FBgn0028519"] <- "FBgn0286723"

### build the rowRanges

suppressPackageStartupMessages(library(TxDb.Dmelanogaster.UCSC.dm6.ensGene))
suppressPackageStartupMessages(library(org.Dm.eg.db))
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
suppressMessages(g <- genes(txdb))
g$symbol <- mapIds(org.Dm.eg.db, names(g), "SYMBOL", "FLYBASE")

# issue with CG4500 / hll
sf6 <- c("CG43110", "veil", "hb", "CG4500", "CG4594", "Cyp4p2", 
         "dream", "CG13384", "CG34266", "l(1)sc", "Adgf-A", "Cht3", "scw",
         "CG3502", "pxb", "CG15628", "CG8960", "CG43085", "path", "CG8147",
         "lea", "Bsg25D", "CG14915", "CG10035", "bmm", "Elba2", "CG14767",
         "prd", "CG15480", "bnb", "comm", "Pino", "fz2", "NimC4", "CG15479",
         "slp1", "mira", "CG14937", "CG14317", "unc-5", "ps", "CG12581", "rst",
         "Surf1", "mas", "Esp", "cnc", "uif", "CG12730", "Ance", "Bsg25A",
         "pcs", "CG13454", "CG43394", "CG9775", "sro", "CG2930",
         "l(2)08717", "scb", "srw", "CG30015", "sprt", "edl", "Mipp1",
         "CG9863", "Pepck")

# 8 are missing
table(sf6 %in% g$symbol)
miss <- sf6[!sf6 %in% g$symbol]

# use one of the files to link missing symbols to Ensembl gene IDs
## url <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2731nnn/GSM2731036/suppl/GSM2731036%5FmelXsim%5Fcyc14C%5Frep1%5Fsl01%2Efpkm%2Etxt%2Egz"
## download.file(url, "fpkm.txt.gz")
fpkm1 <- read.delim("fpkm.txt.gz")
symbol <- fpkm1$gene_short_name
names(symbol) <- fpkm1$gene_id
any(is.na(symbol))
any(symbol == "")

# for later, define ranges that didn't move from dm5.57 -> dm6
# the fpkm data is in dm5.57
library(IRanges)
spl <- strsplit(fpkm1$locus, ":|-")
starts <- as.numeric(sapply(spl, `[`, 2)) + 1 # 0-based
ends <- as.numeric(sapply(spl, `[`, 3))
fpkm_rng <- IRanges(starts, ends)
names(fpkm_rng) <- fpkm1$gene_id

# fix CG4500 / hll
which(names(fpkm_rng) == "FBgn0028519")
fpkm_rng[""] "FBgn0286723"

# 8 are here, so this work to bridge the gap
table(miss %in% symbol)
table(names(symbol[symbol %in% miss]) %in% names(g))

# manually fix CG4500 / hll
which(symbol == "CG4500")
names(symbol)[1840] <- "FBgn0286723"

# fixed
table(names(symbol[symbol %in% miss]) %in% names(g))

g$paper_symbol <- g$symbol
miss_symbol <- symbol[symbol %in% miss]
g[names(miss_symbol)]$paper_symbol <- miss_symbol
g[names(miss_symbol)]

table(sf6 %in% g$paper_symbol)
g$svASE <- g$paper_symbol %in% sf6
g <- g[!is.na(g$symbol)]

### save SummarizedExperiment

common <- intersect(rownames(a1), names(g))
a1 <- a1[common,]
a2 <- a2[common,]
genes <- g[common]
suppressPackageStartupMessages(library(SummarizedExperiment))
coldata <- data.frame(strain=factor(strain), slice, normSlice, rep)
se <- SummarizedExperiment(assays=list(a1=a1, a2=a2), rowRanges=genes, colData=coldata)
colnames(se) <- titles
metadata(se) <- list(
  author="Combs PA, Fraser HB",
  title="Spatially varying cis-regulatory divergence in Drosophila embryos elucidates cis-regulatory logic",
  journal="PLOS Genetics",
  year="2018",
  additional="14(11):e1007631",
  doi="10.1371/journal.pgen.1007631",
  alleles=c(a1="alt / D simulans",
            a2="ref / D melanogaster")
)
mcols(se)
table(mcols(se)$svASE) # 65 out of 66
sum(is.na(mcols(se)$symbol)) # 0
se <- se[,rearrange]
plot(se$slice)
plot(se$rep)
spatialDmelxsim <- se
save(spatialDmelxsim, file="../../spatialDmelxsim_original.rda")

#############

# below requires 'fpkm_rng' and 'symbol'
suppressPackageStartupMessages(library(SummarizedExperiment))
load("../../spatialDmelxsim_original.rda")

# the FPKM matrix
#url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102233/suppl/GSE102233%5Fexpr%5Fsummary%2Etsv%2Egz"
#download.file(url, "expr_summary.tsv.gz")
fpkm <- read.delim("expr_summary.tsv.gz", row.names=1)
idx <- grepl("melXsim",colnames(fpkm)) | grepl("simXmel",colnames(fpkm))
table(idx)
fpkm <- fpkm[,idx]
fpkm <- dplyr::na_if(fpkm, "---")
fpkm <- as.matrix(fpkm)
mode(fpkm) <- "numeric"
colnames(fpkm) <- sub("_FPKM","",colnames(fpkm))
all_na <- apply(fpkm, 2, function(x) all(is.na(x)))
which(all_na)
se <- spatialDmelxsim
# note the 14 slices which had missing FPKM data
se$fpkmMiss <- FALSE
se$fpkmMiss[all_na] <- TRUE

# impute the FPKM values using the average of left- and
# right-most non-missing slices (none were boundary slices)
left <- c(3,3,7,19,22,41,73,86,92,92,92,96,114,129)
right <- c(6,6,9,21,24,43,75,88,96,96,96,98,116,131)
for (j in seq_along(which(all_na))) {
  fpkm[,which(all_na)[j]] <- .5 * fpkm[,left[j]] + .5 * fpkm[,right[j]]
}

# predict total count (the a1 + a2 count, not the gene total count)
assay(se, "total") <- assay(se, "a1") + assay(se, "a2")
common <- intersect(rownames(se), names(fpkm_rng))
fpkm_rng <- fpkm_rng[common]
mcols(se)$matchDm557 <- FALSE
mcols(se)[common,]$matchDm557 <- ranges(rowRanges(se)[common]) == fpkm_rng
table(mcols(se)$matchDm557) # 2959

# change fpkm rownames to Ensembl
geneid <- names(symbol)
names(geneid) <- symbol
rownames(fpkm) <- geneid[rownames(fpkm)]
table(rownames(fpkm) %in% rownames(se)[mcols(se)$matchDm557]) # 2959
fpkm_sub <- fpkm[rownames(se)[mcols(se)$matchDm557],]
nrow(fpkm_sub)
all(colnames(fpkm_sub) == colnames(se))

# gene length
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
suppressMessages(g <- genes(txdb))
txp <- transcripts(txdb)
ebt <- exonsBy(txdb, by="tx")
names(ebt) <- txp$tx_name
txp_len <- sum(width(ebt))
txdf <- select(txdb,
               rownames(se),
               "TXNAME", "GENEID")
txdf$len <- txp_len[txdf$TXNAME]
gene_len <- tapply(txdf$len, txdf$GENEID, mean)
gene_len_nms <- names(gene_len)
gene_len <- as.vector(gene_len)
names(gene_len) <- gene_len_nms
table(rownames(fpkm_sub) %in% names(gene_len))

idx <- mcols(se)$matchDm557
idx2 <- rownames(fpkm_sub)
plot(fpkm_sub[,132]*gene_len[idx2]+1000, assay(se[idx,], "total")[,132]+1, log="xy")
cor_mat <- cor(log10(fpkm_sub*gene_len[idx2]+1000), log10(assay(se[idx,], "total")+1), use="pairwise")
image(cor_mat) # looks ok
median(cor_mat) # .82
coefs <- matrix(nrow=2, ncol=ncol(fpkm_sub))
adjr2 <- numeric(ncol(fpkm_sub))
for (j in seq_len(ncol(fpkm_sub))) {
  tot <- assay(se[idx,], "total")[,j]
  nuc <- fpkm_sub[,j]*gene_len[idx2]
  finite <- tot > 0 & nuc > 0
  fit <- lm(log(tot[finite]) ~ log(nuc[finite]))
  coefs[,j] <- coef(fit)
  adjr2[j] <- summary(fit)$adj.r.squared
}
summary(adjr2)
which.min(adjr2) # 42
summary(t(coefs))
save(coef, adjr2, file="coefs.rda")

fpkm_len <- fpkm[rownames(fpkm) %in% names(gene_len),]
pred_tot <- matrix(nrow=nrow(fpkm_len), ncol=ncol(se))
rownames(pred_tot) <- rownames(fpkm_len)
colnames(pred_tot) <- colnames(fpkm_len)
for (j in seq_len(ncol(fpkm))) {
  finite <- fpkm_len[,j] > 0
  pred_tot[finite,j] <- exp(coefs[1,j] + coefs[2,j] * log(fpkm_len[finite,j] * gene_len[rownames(fpkm_len)][finite]))
}
save(pred_tot, file="pred_tot.rda")

# see how it does on training:
train <- assay(se[mcols(se)$matchDm557,], "total")
j <- 132
plot(train[,j], pred_tot[rownames(train),j], log="xy")

# the ASE matrix
#url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102233/suppl/GSE102233%5Fase%5Fsummary%2Etsv%2Egz"
#download.file(url, "ase_summary.tsv.gz")
ase <- read.delim("ase_summary.tsv.gz", row.names=1)
idx <- grepl("melXsim",colnames(ase)) | grepl("simXmel",colnames(ase))
table(idx)
ase <- ase[,idx]
ase <- dplyr::na_if(ase, "---")
ase <- as.matrix(ase)
mode(ase) <- "numeric"
ase <- (ase + 1)/2
colnames(ase) <- sub("_ase_value","",colnames(ase))
all_na <- apply(ase, 2, function(x) all(is.na(x)))
table(all_na)
which(all_na) # 42 is here also, so we won't impute for that one

# subset to the 10347 genes we don't have counts for (no match with dm5.57)
pred_sub <- pred_tot[intersect(rownames(pred_tot), rownames(se[!mcols(se)$matchDm557])),]
rownames(ase) <- geneid[rownames(ase)]
# subset to the 5672 genes we have predicted total counts for
table(rownames(pred_sub) %in% rownames(ase))
pred_sub <- pred_sub[intersect(rownames(pred_sub), rownames(ase)),]
ase <- ase[rownames(pred_sub),]

pred_a1 <- round(pred_sub * ase)
pred_a2 <- round(pred_sub * (1 - ase))
save(pred_a1, pred_a2, file="pred_allelic_counts.rda")

assay(se, "a1")[rownames(pred_a1),] <- pred_a1
assay(se, "a2")[rownames(pred_a2),] <- pred_a2
mcols(se)$predicted <- FALSE
mcols(se)$predicted[rownames(se) %in% rownames(pred_a1)] <- TRUE
table(mcols(se)$predicted) # 5672
assays(se) <- assays(se)[1:2] # remove total assay

# pull out the FPKM matrix matching SE rows
table(rownames(se) %in% rownames(fpkm))
fpkm_miss <- rownames(se)[!rownames(se) %in% rownames(fpkm)]
fpkm_match <- rbind(fpkm,
                    matrix(0, nrow=length(fpkm_miss), ncol=ncol(fpkm),
                           dimnames=list(fpkm_miss, colnames(fpkm))))
fpkm_match <- fpkm_match[rownames(se),]
all(colnames(fpkm_match) == colnames(se))
all(rownames(fpkm_match) == rownames(se))
tpm <- apply(fpkm_match, 2, function(x) 1e6*x/sum(x))
assay(se, "tpm") <- tpm # save the TPM data

spatialDmelxsim <- se
save(spatialDmelxsim, file="../../spatialDmelxsim.rda")
