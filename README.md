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
