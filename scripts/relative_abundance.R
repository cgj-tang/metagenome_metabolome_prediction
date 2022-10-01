###################################################
###CONVERTING GENE COUNTS TO RELATIVE ABUNDANCES###
###################################################

#This function converts raw counts of microbial functions into relative
#abundances. It accepts as input a matrix or dataframe containing unnormalized
#gene abundances and returns a table of the same size containing abundances
#normalized by sample. All NA values are converted into zeroes in the process.
#Samples must be arranged by row. 

#If sample identifiers are held in a single column rather than row names, 
#provide the column index in ID_index, else the function will throw back an 
#error.


relative_abundance <- function(dataset, 
                               ID_index = NULL) {
        if (!is.null(ID_index)) {
                        ID <- dataset[, ID_index]
                        dataset <- dataset[, -ID_index]
        }
        if (sum(is.na(dataset)) > 0) {
                print('Changing all NA values to zero')
                dataset[is.na(dataset)] <- 0
        }
        dataset_norm <- sweep(dataset, MARGIN = 1, STATS = rowSums(dataset), FUN = '/')
        if (!is.null(ID_index)) {
                dataset_norm <- cbind(ID, dataset_norm)
        }
        return(dataset_norm)
}