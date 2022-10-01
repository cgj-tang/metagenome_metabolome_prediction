#####################################
###LOG AND AUTO-SCALING ABUNDANCES###
#####################################

log_autoscale <- function(dataset,
                          impute = FALSE,
                          glog = FALSE,
                          ID_index = NULL) {
        if (!is.null(ID_index)) {
                ID <- dataset[, ID_index]
                dataset <- dataset[, -ID_index]
        }
        if (sum(is.na(dataset)) > 0) {
                print('Changing all NA values to zero')
                dataset[is.na(dataset)] <- 0
        }
        if (impute == TRUE) {
                for (i in 1:ncol(dataset)) {
                        min_nonzero <- min(dataset[, i][dataset[, i] > 0])
                        dataset[, i][dataset[, i] == 0] <- min_nonzero / 4
                }
        }
        if (glog == TRUE) {
                min_dataset <- min(dataset)
                dataset_logged <- log2((dataset + sqrt(dataset^2 + min_dataset^2)) / 2)
        }
        else {
                dataset_logged <- log10(dataset)
        }
        dataset_scaled <- as.data.frame(scale(dataset_logged))
        if (!is.null(ID_index)) {
                dataset_scaled <- cbind(ID, dataset_scaled)
        }
        return(dataset_scaled)
}