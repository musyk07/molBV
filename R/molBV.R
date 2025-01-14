#' Perform the calculation of the molBV scores
#'
#' @import phyloseq
#' @import ggplot2
#' @import ggpubr
#'
#'
#' @param inputData Either a biom file or a matrix-like R object with sample IDs as columns and genera as rows.
#' @param outputTable The path and name of the output table to be generated. Default is NULL.
#'
#' @details
#' \itemize{
#'   \item If \code{inputData} is biom file, package 'phyloseq' is needed to process it. Some warnings related to the greengenes database can be ignored.
#'   \item If \code{inputData} is a matrix/data frame, please make sure that the genera labels (row names) match the green-genes convention, which is necessary when calculating the score for each reference frame.
#'   \item If \code{outputTable} is specified, a table contains the sample map with the molBV socres will be generated.
#' }
#'
#' @return A data frame contains the sample map with the molBV socres.
#' \itemize{
#'   \item SID - Sample IDs.
#'   \item Column 2 to 12 - Counts of the genus that used to calculated the scores.
#'   \item molBV - Calculated molBV scores.
#' }
#'
#' @references Usyk et al. The Immune Landscape of Molecular Bacterial Vaginosis and HPV Natural History. https://doi.org/10.21203/rs.3.rs-613095/v1
#'
#' @examples
#' map_with_molBVScores <- molBV(exampleMicrobiomeData, outputTable = NULL)
#' print(map_with_molBVScores)
#' range(map_with_molBVScores$molBV)
#'
#' exampleBiomFile <- system.file("extdata", "converted_otu.biom", package = "molBV")
#' map_with_molBVScores <- molBV(exampleBiomFile, outputTable = "map_with_molBV.txt")
#' #A txt file "map_with_molBV.txt" will be generated in current working directory.
#'
#' @export

molBV <- function(inputData, outputTable = NULL){
	if(is.character(inputData)){
		# load the biom file:
		clst.expt <- import_biom(BIOMfilename = inputData, parseFunction = parse_taxonomy_greengenes)
		#warnings ignored

		# aggregate to genus level:
		clst.expt_gen <- tax_glom(physeq = clst.expt, taxrank = "Genus")

		# rename duplicated genus names:
		taxa_names(clst.expt_gen) <- make.unique(tax_table(clst.expt_gen)[, "Genus"])

		mbiome_data_in <- as.data.frame(as.matrix(otu_table(clst.expt_gen)), check.names = FALSE)

	}else if(is.matrix(inputData) || is.data.frame(inputData)){
		if(any(is.na(inputData) | inputData < 0)){
			inputData[is.na(inputData) | inputData < 0] <- 0
			warning("Negative or N/A values detected in the inputData! set to 0.")
		}
		mbiome_data_in <- as.data.frame(inputData, check.names = FALSE)

	}else stop("'inputData' should be either a path to a biom file or a matrix-like R object!")

	# importing the beta coefficients for each BV association ratio from dataset 'molecuL'
	if(!exists("molecuL"))data(molecuL)

	# making a data frame to store the scores
	map_in <- data.frame(SID = colnames(mbiome_data_in), stringsAsFactors = FALSE)

	# setting up a vector to loop over
	top_bogsels <- c("Lactobacillus", molecuL$bug2)

	# this is a function to make it easier to get specific taxa out
	bug_returner <- function(bug = "Gardnerella", samp_order, df_in = mbiome_data_in){
		if(sum(!samp_order %in% colnames(df_in)) > 1){
			df_in[, samp_order[!samp_order %in% colnames(df_in)]] <- 0
		}
		get_bug <- grepl(bug, rownames(df_in))
		return(colSums(df_in[get_bug, samp_order]))
	}

	# apply the above function to load all of the bugs into the map
	for(i in 1:length(top_bogsels))map_in[[top_bogsels[i]]] <- bug_returner(bug = top_bogsels[i], samp_order = map_in$SID)

	# adding the entry for molBV and apply the formula to calculated the scores
	map_in$molBV <- NA
	for(i in 1:dim(map_in)[1]){
		calci_vec <- sapply(molecuL$bug2, FUN = function(x)molecuL[x, "intercept"] + log((map_in[i, "Lactobacillus"] + 1) / (map_in[i, x] + 1)) * molecuL[x, "beta"])
		map_in$molBV[i] <- mean(calci_vec)
	}

	# truncate the results of molBV at 10 to match Nugent range
	map_in$molBV[map_in$molBV > 10] <- 10

	# write out the sample map with the molBV score
	if(!is.null(outputTable))write.table(map_in, file = outputTable, quote = FALSE, sep = "\t", row.names = FALSE)

	# the sample map with molBV scores is returned
	return(map_in)
}
