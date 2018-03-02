# metaxcanr
An R alternative to MetaXcan 

## more details on method
for more details on metaxcan method refer to Metaxcan/Predixcan authors' repo (https://github.com/hakyimlab/MetaXcan)

## installation
```
library(devtools)
install_github("drveera/metaxcanr")
```

## usage 
```
library(metaxcanr)

results <- metaxcan(gwas.file = <gwas.filename>, db.file = <dbfile.name>, snpcov.file = <snpcov.filename>)
```

## examples

### download sample files

```
wget https://github.com/drveera/metaxcanexample/raw/master/sample.db
wget https://raw.githubusercontent.com/drveera/metaxcanexample/master/sample.gwas.summary
wget https://github.com/drveera/metaxcanexample/raw/master/sample.snpcov
```

### open R and run the following commands

```
library(metaxcanr)
results <- metaxcan(gwas.file="sample.gwas.summary", db.file = "sample.db", snpcov.file="sample.snpcov")

```

## Format of GWAS data

The gwas data should have columns: SNP, A1, A2, BETA, SE 

If OR, change that to BETA by taking log(OR)

**Its very important the above columns should be there in the gwas file**
 
## For Citation

1. Barbeira, A. N. et al. Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics. bioRxiv 045260 (2017). doi:10.1101/045260

(link: https://www.biorxiv.org/content/early/2017/10/03/045260)
