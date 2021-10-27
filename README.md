# molBV
This is a brief instruction for the "rlm_molBV.R" script. The purpose of the script is to perform the calculations described in the manuscript tittled "The Immune Landscape of Molecular Bacterial Vaginosis and HPV Natural History".

To run the script you should have R version 4.2.0 installed. Additionally if you will be importing data using the biom format, you will need to install the phyloseq package. The script also uses the ggplot2 and ggpubr packages to visualize the association of molBV with clinical variables and you will need to install these packages if you would like to generate the boxplots. Typical installation time for the packages should take under 10 minutes on MacBook Pro laptop with an i7 processor and 16 Gb of RAM.  

To run the script, open it within Rstudio (or which ever IDE you prefer) and set the working directory where you have the contents of the script directory saved. Then you can run through the calculation using the original Mt. Sinai BV dataset (files: "converted_otu.biom" and "Mt_Sinai_sample_map.txt"). 

If you wish to use the script with your own data, please import using phyloseq as is done in the main script. If your data is already stored in a matrix with sample IDs as columns and genera as rows, you may wish to directly import into the script on line 17. Please make sure that the data frame looks like the testing dataset before proceeding. Also make sure that the genera labels match the green-genes convention, which is necessary when calculating the score for each reference frame. 

Further details about specific steps of the script can be found as comments throughout the script. 
