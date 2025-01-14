# molBV
This is a brief instruction showing several options for installing and running the molBV algorithm (which calculated Nugent-like scores ranging from 1-10 using 16S rRNA data sequenced from cervicovaginal samples). See original publication for full details, validation and connection between molBV and cervicovaignal cytokines: https://www.nature.com/articles/s41467-021-27628-3

# molBV Package
## Installation
Download the source file `molBV_1.0.tar.gz` and run the following R commands:
```R
install.packages("molBV_1.0.tar.gz", repos = NULL, type = "source")
browseVignettes("molBV")
??molBV
```
This will install the package and allow you to browse a sample tutorial using the above command. 
At the bottom of this page there is additional instructions for running the molBV algorithm manually (i.e. see "molBV (manual instructions)", which is more clunky, but you get to see all of the individual commands). 


Please cite original paper if using molBV:

Usyk, M., Schlecht, N. F., Pickering, S., Williams, L., Sollecito, C. C., & Gradissimo, A. molBV reveals immune landscape of bacterial vaginosis and predicts human papillomavirus infection natural history. Nat Commun. 2022; 13: 233.
https://www.nature.com/articles/s41467-021-27628-3 

# Example Usage 
###### created by Usyk, Mykhaylo: mykhaylo.usyk@einsteinmed.edu and Ling, Wodan: wling@fredhutch.org
## Load the Package and Required Libraries

```R
library(molBV)
library(phyloseq)
library(ggplot2)
library(ggpubr)
```
## Load Example Data

```R
data(exampleMicrobiomeData)
dim(exampleMicrobiomeData)
head(exampleMicrobiomeData)
range(exampleMicrobiomeData)
```

## Calculate molBV Scores

```R
map_with_molBVScores <- molBV(exampleMicrobiomeData, outputTable = NULL)
str(map_with_molBVScores)
range(map_with_molBVScores$molBV)

```
## Use BIOM File for Input
```R
exampleBiomFile <- system.file("extdata", "converted_otu.biom", package = "molBV")
map_with_molBVScores <- molBV(exampleBiomFile, outputTable = "map_with_molBV.txt")

```
## Merge with Clinical Metadata

```R
exampleMetaData <- system.file("extdata", "Mt_Sinai_sample_map.txt", package = "molBV")
clinical_map <- read.table(exampleMetaData, header = TRUE, sep = "\t", as.is = TRUE)
str(clinical_map)

full_set <- merge(x = map_with_molBVScores, y = clinical_map, by.x = "SID", by.y = "Cerv_SampleID")
```
## Visualize molBV Scores compared to Clinical BV diagnosis 
```R
basic_boxes <- function(df_in=full_set, x_variable="Amsel_Diagnosis"){
	p <- ggplot(df_in, aes_string(x = x_variable, y = "molBV", fill = x_variable)) +
		geom_boxplot()+
		theme_minimal() + stat_compare_means() +
		ggtitle(label = paste("molBV across", x_variable))
	print(p)
}

basic_boxes()  # molBV across Amsel_Diagnosis
basic_boxes(x_variable = "Nugent_diagnosis")  # molBV across Nugent_diagnosis

```

# molBV (manual instructions)
This is a brief instruction for the "rlm_molBV.R" script. The purpose of the script is to perform the calculations described in the manuscript tittled "molBV reveals immune landscape of bacterial vaginosis and predicts human papillomavirus infection natural history".

To run the script you should have R version 4.2.0 installed. Additionally if you will be importing data using the biom format, you will need to install the phyloseq package. The script also uses the ggplot2 and ggpubr packages to visualize the association of molBV with clinical variables and you will need to install these packages if you would like to generate the boxplots. Typical installation time for the packages should take under 10 minutes on MacBook Pro laptop with an i7 processor and 16 Gb of RAM.  


Please cite original paper if using molBV:
Usyk, M., Schlecht, N. F., Pickering, S., Williams, L., Sollecito, C. C., & Gradissimo, A. molBV reveals immune landscape of bacterial vaginosis and predicts human papillomavirus infection natural history. Nat Commun. 2022; 13: 233.
https://www.nature.com/articles/s41467-021-27628-3 

