##########################################
###FILTERING OUT LOW ABUNDANCE FEATURES###
##########################################

#This function filters out low abundance features found within an entire
#dataset. It accepts as input a matrix or dataframe containing features
#normalized to relative abundances and returns a table without the filtered
#features. For instance, running this function under the default arguments will
#remove all features that are not at least 0.01% abundant across 10% of all
#samples. Samples must be arranged by row. 

#If sample identifiers are held in a single column rather than row names,
#provide the column index in ID_index. This is not 100% necessary unless you
#want to renormalize the dataset after filtering.

#Set renormalize = TRUE if you want to renormalize features after filtering.



abundance_filter <- function(dataset, 
                             dataset_threshold = 0.10, 
                             sample_threshold = 0.0001, 
                             ID_index = NULL,
                             renormalize = FALSE) {
        if (dataset_threshold == 0){
                print("No filtering will be applied due to dataset_threshold being set to zero.")
                return(dataset)
        }
        if (!is.null(ID_index)) {
                ID <- dataset[, ID_index]
                dataset <- dataset[, -ID_index]
        }
        if (sum(is.na(dataset)) > 0) {
                print('Changing all NA values to zero')
                dataset[is.na(dataset)] <- 0
        }
        retained <- c()
        threshold <- ceiling(dataset_threshold * nrow(dataset))
        for (i in 1:ncol(dataset) ) {
                counter <- sum(dataset[, i] > sample_threshold, na.rm = TRUE)
                if (counter > threshold) {  
                        retained <- c(retained , i)
                }
        }
        if (renormalize == TRUE) {
                dataset_filt <- dataset[, retained]
                dataset_filt <- relative_abundance(dataset_filt)
                }
        else {
                dataset_filt <- dataset[, retained]
        }
        if (!is.null(ID_index)) {
                dataset_filt <- cbind(ID, dataset_filt)
        }
        print(paste(ncol(dataset_filt), 'out of', ncol(dataset), 'columns were retained.'))
        return(dataset_filt)
}