# README

This file contains information regarding documentation of certain scripts used in this thesis. Shell scripts have self-contained documentation, and you may refer to the file headers for such information on usage.

## dge-analysis.r Documentation

### txt2counts()

Reads a featurecounts txt file into an R data frame
Requires:

* featurecounts.txt: Name/path to a featurecounts count file
    
Return value: Flat matrix of counts

### dge_analysis()
Given either A) a `featurecounts` flat text file or B) a `data frame` like object containing raw read counts, and a design matrix `coldata`, `dge_analysis` will automate the process of calling DESeq and perform all pairwise combinations of `contrasts` (by default, the levels of elements in the design matrix column corresponding to `contrast.var` - see below) for DGE hypothesis testing.

The function (by default) will assume that all columns of the design matrix are to be used in the call to `DESeqDataSetFromMatrix`, with the order of the columns being the order of main effects. E.g., suppose that a user passes the following design matrix:

```
coldata
  Culture Type
1     786  MET
2     786  MET
3     786  PRT
4     786  PRT
5      OS  MET
6      OS  MET
7      OS  PRT
8      OS  PRT
```
The user may optionally pass a vector for all elements to be used in the formula for DESeq - for the formula `~ Culture + Type + Culture:Type`, one would pass `formula.vec = c("Culture", "Type", "Culture:Type")`. If no formula vector is passed to `dge_analysis()`, the function assumes that the columns of the design matrix indicate the design. Therefore, for the above `coldata`, the formula used in the function will be `Culture + Type`, searching for DGE as a main effect of Type.

Currently, the function only supports DGE as a main effect of a single variable (`contrast.var`), which can be passed to the function by the user, but by default it is simply the last column in `coldata`. Future support for performing pairwise DGE for multiple variables (e.g., discovering differentially expressed genes for both `786 vs. OS` and `MET vs. PET` above) will hopefully be available in the future.

Additionally, this function assumes that `contrasts` are the levels of the design matrix with column name given by `contrast.var`. The levels are assumed to be in order of numerator -> ... -> denominator. In other words, given a leveling scheme of `levels = c("l1", "l2", "l3")`, `dge_analysis` will perform DGE for `l1 vs. l2`, `l1 vs. l3`, and `l2 vs. l3`. This is contradictory to how DESeq2 suggests to level, since `dge_analysis()` assumes the reference level is the "highest" level. More work to make this in-line with the DESeq2 vignette is forthcoming.