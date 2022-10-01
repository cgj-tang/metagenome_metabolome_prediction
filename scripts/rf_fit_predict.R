##########################################################
###RANDOM FOREST METABOLOMICS MACHINE LEARNING PIPELINE###
##########################################################

#This function encodes the random forest regression pipeline used to train
#models to predict metabolite abundances using microbial sequencing features
#(e.g., taxonomic or functional abundances). It accepts accepts two lists, one
#for the metabolite data another for the features to be used as predictors. 
#Each list must contain two elements each, corresponding to two separate
#metabolite-sequencing feature datasets. For each metabolite, a model will be
#trained using one dataset, and used to predict metabolite abundances in the
#other. This function will output metabolite abundance predictions using the
#sequencing features derived from each sample.

#Most of the time will probably be spent fitting models optimizing the 
#hyperparameters. Note to self: make hyperparameter tuning optional and add in
#argument to set the number of cross-validation folds.

#Elements across both lists must share the exact same names in the same order.
#It is intended that list elements be named after the datasets they were derived
#from.

#The first column in each element must contain sample IDs. Sample IDs in one
#element must match the sample IDs of the corresponding element of the same
#name in the other list. For instance then, an element named 'Lloyd' in the
#metabolite list must contain sample IDs that match those found in the element
#bearing the same name in the features list.

#All elements in both lists must be arranged with samples in rows and features
#in columns.

#Features in the metabolite list (i.e., metabolites) must be shared across all
#elements in that list. 

#You may set a custom seed for the model if you like.

#Required libraries
require(tidymodels) #Model fitting
require(tidyverse)

rf_fit_predict <- function(metab_list, 
                           feature_list, 
                           custom_seed = 123,
                           tuning = TRUE) {
        
        #Catching any potential errors
        if (length(metab_list) != 2 | length(feature_list) != 2) {
                stop('Error: Only two datasets are accepted.')
        }
        if (!all((names(metab_list)) == names(feature_list))) {
                stop('Error: Elements in both lists must match by name.')
        }
        if (!all((lapply(metab_list, colnames)[[1]] == lapply(metab_list, colnames)[[2]]))) {
                stop('Error: Metabolites across both datasets must match.')
        }
        if(!all((metab_list[[1]][[1]] == feature_list[[1]][[1]]) | !all(metab_list[[2]][[1]] == feature_list[[2]][[1]]))) {
                stop('Error: Sample names must match.')
        }
        
        #Initializing various objects
        pred_list <- list() #Empty list to hold predictions
        dataset_names <- names(metab_list) #Vector of dataset names to iterate over
        metabolite_names <- colnames(metab_list[[dataset_names[1]]])[-1] #Gathering metabolite names
        
        for (i in 1:length(dataset_names)) { #Iterating over datasets
                
                k <- ifelse(i == 1, 2, 1)
                
                #Empty dataframe containing sample IDs to capture predictions
                assign(paste(dataset_names[k], '_pred', sep = ''),
                       metab_list[[dataset_names[k]]][1])
                
                for (j in 1:length(metabolite_names)) { #Iterating over metabolites
                        
                        #System messages
                        start_time <- Sys.time()
                        print(paste('Fitting model for ', metabolite_names[j], ' (', j, ')', ' using dataset: ', dataset_names[i], sep = ''))
                        
                        #Assembling training data
                        train_dt <- 
                                cbind(metab_list[[dataset_names[i]]][c(1, j + 1)], 
                                      feature_list[[dataset_names[i]]][-1]) %>%
                                rename(Metabolite = 2) %>%
                                mutate(Case = 'train')
                        
                        #Assembling test data
                        test_dt <- 
                                cbind(metab_list[[dataset_names[k]]][c(1, j + 1)], 
                                      feature_list[[dataset_names[k]]][-1]) %>%
                                rename(Metabolite = 2) %>%
                                mutate(Case = 'test')
                        
                        #Combining training and testing data together
                        dt <- 
                                train_dt %>%
                                bind_rows(test_dt) %>%
                                replace(is.na(.), 
                                        0)
                        
                        #Splitting data into training and test sets
                        data_split <- 
                                make_splits(dt %>% filter(Case == 'train'),
                                            dt %>% filter(Case == 'test'))
                        train_data <- training(data_split)
                        test_data <- testing(data_split)
                        rf_folds <- vfold_cv(train_data, v = 10) #Cross-validation folds
                        
                        #Initializing random forest model and recipe
                        rf_rec <- 
                                recipe(Metabolite ~ ., data = train_data) %>%
                                update_role(ID, Case, new_role = 'Sample_IDs') %>% #Masking ID and Case columns from the model
                                step_zv(all_predictors()) #Removing zero variance features
                        
                        if (tuning == FALSE) {
                                rf_mod <- 
                                        rand_forest(trees = 1001) %>% #Odd number to resolve tie votes
                                        set_engine('ranger', 
                                                   seed = custom_seed) %>%
                                        set_mode('regression')
                        }
                        else {
                                rf_mod <- 
                                        rand_forest(mtry = tune(),
                                                    min_n = tune(),
                                                    trees = 1001) %>% #Odd number to resolve tie votes
                                        set_engine('ranger', 
                                                   seed = custom_seed) %>%
                                        set_mode('regression')
                        }
                        
                        #Assembling workflow
                        rf_wflow <- 
                                workflow() %>%
                                add_model(rf_mod) %>%
                                add_recipe(rf_rec)
                        
                        #Skip right to fitting and predicting if no tuning is to occur
                        if (tuning == FALSE) {
                                #Predicting metabolite abundances in test set
                                final_fit <- 
                                        rf_wflow %>%
                                        fit(data = train_data)
                                
                                #Collecting predictions
                                final_pred <-
                                        predict(final_fit, test_data) %>%
                                        as.data.frame()
                                assign(paste(dataset_names[k], '_pred', sep = ''),
                                       get(paste(dataset_names[k], '_pred', sep = '')) %>%
                                               mutate(!!metabolite_names[j] := final_pred[, 1]))
                                
                                #System messages
                                print(paste('Generated predicted metabolite abundance for ', metabolite_names[j], ' (', j, ')', ' in dataset: ', dataset_names[k], sep = ''))
                                end_time <- Sys.time()
                                print(end_time - start_time)
                                next
                        }
                        
                        #Tuning hyperparameters
                        rf_res <- 
                                rf_wflow %>%
                                tune_grid(resamples = rf_folds,
                                          grid = 10,
                                          control = control_grid(save_pred = TRUE))
                        rf_best <- 
                                rf_res %>%
                                select_best(metric = 'rmse') #Select optimal hyperparameters
                        
                        #Finalizing workflow and predicting metabolite abundances in test set
                        final_wflow <- 
                                rf_wflow %>%
                                finalize_workflow(rf_best)
                        final_res <- 
                                final_wflow %>%
                                last_fit(data_split)
                        
                        #Collecting predictions
                        #Note to self: Find a more memory-efficient way to do this
                        final_pred <- 
                                final_res %>%
                                collect_predictions()
                        assign(paste(dataset_names[k], '_pred', sep = ''),
                               get(paste(dataset_names[k], '_pred', sep = '')) %>%
                                       mutate(!!metabolite_names[j] := final_pred$.pred))
                        
                        #System messages
                        print(paste('Generated predicted metabolite abundance for ', metabolite_names[j], ' (', j, ')', ' in dataset: ', dataset_names[k], sep = ''))
                        end_time <- Sys.time()
                        print(end_time - start_time)
                }
                
                #Populating empty list
                #Note to self: Find a more memory-efficient way to do this
                pred_list[[dataset_names[k]]] <- get(paste(dataset_names[k], '_pred', sep = ''))
        }
        return(pred_list)
}