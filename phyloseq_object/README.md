# Phyloseq Object 1 of ASVs
The PhyloSeq object with ASVs from all bacterial phyla used for analysis is exported as: ps_HV_site_time_final.ps.rds. This object was used for Figure S2. There are two important points:
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

# Phyloseq Object 2 of ASVs
The PhyloSeq object with ASVs from only the Staphylococcus genus used for analysis is exported as: ps_staph_ASV_filter_final.rds. This object was used for Figure 1 and Figure S3. The object retains only those Staphylococcus ASVs that fulfil the following two criteria:
1. ASV was present at >= 1% abundance in atleast one sample
2. ASV was present in > 3 individuals, irrespective of the body site at which it was detected 
