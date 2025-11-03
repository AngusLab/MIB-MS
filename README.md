# MIB-MS

MIB-MS kinome analysis package and protocols

This package uses either MaxQuant LFQ or DIA-NN protein level data and filters for: 1) kinases (either human or mouse),2) proteins with >1 unique peptide, and 3) kinases quantified in all replicates of at least one sample.

Users must also create a metadata file to match the protein column names with the correct treatment replicate. An example can be found in the data folder.

To force the correct comparison, users may put "Z." in front of the control sample treatment name, thus creating the Treatment vs. Control comparison.

For DIA-NN, the report.parquet file is also required.

Below is how to set up the env and run the script:


```r
library(mibms)
library(arrow)

data("mouse.kinome") #Load pre-installed kinome data


df<-read.csv("proteinGroups.csv", check.names = F)
metadata<-read.csv("Sample.csv")
peptide<- read_parquet("report.parquet")
kinome.analysis<- diann.cleanup(df = df, sample = metadata, mouse.kinome)
stats<- statistical.testing(DF = kinome.analysis)
```


To install in R, please run: devtools::install_github("AngusLab/MIB-MS")

