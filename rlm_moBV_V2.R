# set the working directory to where the contents of the script are
setwd("~/Desktop/nature_communications_script")

# load the required packages 
# Only the phyloseq is necessary for the core calculation if using the biom format for input
library(ggplot2)
library(phyloseq)
library(ggpubr)

# load the data
clst.expt <- import_biom(BIOMfilename = "test_data/converted_otu.biom", parseFunction = parse_taxonomy_greengenes)
clst.expt_gen <- tax_glom(physeq = clst.expt, taxrank = "Genus")
taxa_names(clst.expt_gen) <- make.unique(tax_table(clst.expt_gen)[,"Genus"])

# Extract the OTU table
# If you have you'r data as a standard matrix with genus names as rows and samples as columns, you cna just import here
mbiome_data_in <- as.data.frame(as.matrix(otu_table(clst.expt_gen)))

# This is a function to make it easier to get specific taxa out
bug_returner <- function(bug = "Gardnerella", samp_order, df_in = mbiome_data_in){
  if(sum(!samp_order%in%colnames(df_in))>1){
    df_in[,samp_order[!samp_order%in%colnames(df_in)]] <- 0
  }
  get_bug <- grepl(bug, rownames(df_in))
  return(colSums(df_in[get_bug, samp_order]))
}


# Importing the beta coefficients for each BV association ratio
molecuL <- read.csv(file = "dependencies/robust_regression_coefficient_map.csv", as.is = T) 
rownames(molecuL) <- molecuL$bug2

# Making a df to store the scores
map_in <- as.data.frame(sample_names(clst.expt_gen))
colnames(map_in)[1] <- "SID" 

# Setting up a vector to loop over
top_bogsels <- c("Lactobacillus", molecuL$bug2)

# now to get out beta coefficients
# load all of them pugs into the map
for(i in 1:length(top_bogsels)){
  map_in[[top_bogsels[i]]] <- bug_returner(bug = top_bogsels[i], samp_order = map_in[[1]])
}

# Adding the entry for molBV
map_in$molBV <- NA

for(i in 1:dim(map_in)[1]){
  calci_vec <- sapply(molecuL$bug2, FUN = function(x) molecuL[x,"intercept"] + 
                        log((map_in[i,"Lactobacillus"]+1)/(map_in[i,x]+1))*molecuL[x,"beta"])
  map_in$molBV[i] <- mean(calci_vec)
}

# truncate the results at molBV 10 to match Nugent range
# map_in$molBV <- ifelse(test = map_in$molBV>10, yes = 10, no = map_in$molBV)
# scaling based on Reviewer #2
map_in$molBV <-  10*(map_in$molBV - min(map_in$molBV, na.rm = T))/(max(map_in$molBV, na.rm = T)-min(map_in$molBV, na.rm = T))

# Write out a sample map with the molBV socre 
write.table(x = map_in, file = "expected_output/map_with_molBV.txt", quote = F, sep = "\t", row.names = F)

##############################################################################################################
############ This part is just to visualize the distribution of molBV across Clinical BV diagnoses ###########
##############################################################################################################

# import clinical map 
clinical_map <- read.csv(file = "test_data/Mt_Sinai_sample_map.txt", sep = "\t", as.is = T)

# merge the two dfs
full_set <- merge(x = map_in, y = clinical_map, by.x = "SID", by.y = "Cerv_SampleID")

# make a nice plotter function to see the distribution of molBV across amsel and nugent diagnoses
basic_boxes <- function(df_in=full_set, x_variable="Amsel_Diagnosis"){
  p <- ggplot(df_in, aes_string(x = x_variable, y = "molBV", fill = x_variable)) +
    geom_boxplot()+
    theme_minimal() + stat_compare_means() + 
    ggtitle(label = paste("molBV across", x_variable))
  print(p)
}

# plot for Amsel
basic_boxes()

# Plot for Nugent Diagnoses
basic_boxes(x_variable = "Nugent_diagnosis")




