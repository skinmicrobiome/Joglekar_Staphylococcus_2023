# Phyloseq Object of ASVs
The PhyloSeq object used for analysis is exported as: ps_HV.PJ.recoded.rds. There are two important points:
1. The Subject IDs correspond to the BioSampleID column of Joglekar_amplicon.xlsx, **not** the HV# used in the manuscript
2. BioSampleID 1209 was not used in the analysis due to anomalously high _S. aureus_ burden 

To use this data object you will need a copy of R and the DADA2 and PhyloSeq packages. You can load the object and export the sample metadata using something like:

````
working_path <- "/Users/username/Desktop"
psobject <- "ps_HV.PJ.recoded.rds"
ps<-readRDS(file=file.path(working_path,psobject))
samp<-sample_data(ps)
head(data.frame(samp))
````
