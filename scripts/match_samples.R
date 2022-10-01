#Samples must be arranged by row
#Accepts as input two lists containing individual dataframes with sample x feature abundance data
#Sample names must be in first column of each dataframe
#Matches dataframes in lists

match_samples <- function(list1, list2) {
        if(!all(names(list1) == names(list2))) {
                break('Error: Names of in both lists must match exactly.')
        }
        list_names <- names(list1)
        sample_list <- list()
        for (i in 1:length(list_names)) {
                list1_samples <- list1[[list_names[i]]][[1]]
                list2_samples <- list2[[list_names[i]]][[1]]
                common_samples <- intersect(list1_samples, list2_samples)
                sample_list[[list_names[i]]] <- common_samples
        }
        return(sample_list)
}