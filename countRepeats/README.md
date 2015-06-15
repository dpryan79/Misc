`countRepeats` uses libGTF and HTSlib to process groups of alignments to the genome produced by STAR. STAR can be run as follows to produce up to 1000 alignments for each entry:

    STAR ...threading and file options... --outSAMstrandField intronMotif \
        --outFilterMultimapNmax 1000 --outFilterMismatchNoverLmax 0.1 \
        --outSAMattributes Standard --outSAMunmapped Within

The resulting alignments are processed as follows:

1. Alignments are grouped according to their read name, so that all alignments produced by a read are processed at once.
2. If any alignment from a group overlaps a gene (technically an `exon`) then the entire group is excluded from further analysis.
3. The remaining groups are then overlapped with the UCSC repeatMasker track, after its chromosome names have been converted to Ensembl.
    1. Any group with no repeat overlaps is ignored
    2. Any group that overlaps more than one repeat is ignored
    3. If a group overlaps only a single repeat type (repName, repClass, or repFamily), then the appropriate count is incremented.

The resulting unique counts are written to standard output. They consist of three columns: repeat name, count, and repeat type (repName, repClass, or repFamily).

These can be imported into DESeq2 with the following command:

```r
fromRepeats<- function(sampleTable, type = "repFamily", directory = "", design,
    ignoreRank = FALSE, ...)
{
    if (missing(design)) {
        stop("design is missing")
    }
    l <- lapply(as.character(sampleTable[, 2]), function(fn) read.table(file.path(directory, fn)))
    if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1))))
        stop("Gene IDs (first column) differ between files.")
    tbl <- sapply(l, function(a) a$V2)
    tbl <- tbl[which(read.table(file.path(directory, sampleTable[1,2]))[,3] == type),]
    colnames(tbl) <- sampleTable[, 1]
    rownames(tbl) <- l[[1]]$V1[which(read.table(file.path(directory,
        sampleTable[1,2]))[,3] == type)]
    rownames(sampleTable) <- sampleTable[, 1]
    oldSpecialNames <- c("no_feature", "ambiguous", "too_low_aQual",
        "not_aligned", "alignment_not_unique")
    specialRows <- (substr(rownames(tbl), 1, 1) == "_") | rownames(tbl) %in%
        oldSpecialNames
    tbl <- tbl[!specialRows, , drop = FALSE]
    dds <- DESeqDataSetFromMatrix(countData = tbl, colData = sampleTable[,
        -(1:2), drop = FALSE], design = design, ignoreRank, ...)
    return(dds)
}

ddsName <- fromRepeats(sampleTable, type="repName", directory="~/", design=~Group)
```

`repName` can alternatively be `repClass` or `repFamily`.


