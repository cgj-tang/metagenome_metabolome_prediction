#####################################
###DIFFERENTIAL ABUNDANCE FUNCTION###
#####################################

#This function identifies differentially abundant metabolites through t-test
#with FDR correction
#TO DO: finish writing up comments

differential_abundance <- function(dataset, 
                                   ID_index = NULL, 
                                   dep_index = NULL,
                                   correction_method = 'fdr') {
        if (is.null(dep_index)) {
                stop('Error: Please specify the column index of your dependent variable')
        }
        pv <- 
                dataset[-c(ID_index, dep_index)] %>%
                lapply(function(x) t.test(x ~ dataset[[dep_index]])) %>%
                lapply('[', 'p.value') %>%
                unlist()
        diff_feat <- colnames(dataset[-c(ID_index, dep_index)])[which(p.adjust(pv, method = correction_method) < 0.05)]
        return(diff_feat)
}