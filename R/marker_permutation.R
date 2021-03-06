#' marker_permutation
#'
#' @description Creates random combinations of phenotypes by shuffling markers and
#' calculates the enrichment and depletion p values
#' @param sce_object SingleCellExperiment object in the form of output from format_image_to_sce
#' @param tumour_marker String with the name of the marker used for tumour cells
#' @param num_iter Integer specifying the number of iterations for bootstrapping
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom stats complete.cases
#' @importFrom utils combn
#' @export

marker_permutation <- function(sce_object, num_iter) {

    formatted_data <- data.frame(colData(sce_object))
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    expression_matrix <- assay(sce_object)

    markers <- rownames(expression_matrix)
    cell_ids <- colnames(expression_matrix)

    rownames(expression_matrix) <- NULL
    colnames(expression_matrix) <- NULL
    expression_matrix_t <- t(expression_matrix)
    expression_df <- data.frame(expression_matrix_t)
    colnames(expression_df) <- markers

    formatted_data <- cbind(formatted_data, expression_df)
    formatted_data <- formatted_data[complete.cases(formatted_data),]

    markers <- markers[markers != "DAPI"]

    #generate all combinations of markers into a vector
    marker_combinations <- vector()
    for(i in 1:length(markers)) {
        comb_matrix <- combn(markers, i)
        for(j in 1:ncol(comb_matrix)) {
            comb <- paste0(comb_matrix[,j], collapse = '', sep=',')
            comb <- gsub(",$", "", comb)
            marker_combinations <- c(marker_combinations, comb)
        }
    }

    #create the results df to store the output of every bootstrap iteration
    results <- data.frame(matrix(0, nrow = length(marker_combinations), ncol = num_iter))
    rownames(results) <- marker_combinations
    colnames(results) <- 1:num_iter

    #count the markers and put counts into a df
    marker_count <- vector()
    for (marker in markers) {
        count <- nrow(formatted_data[grepl(marker, formatted_data$Phenotype), ])
        marker_count <- c(marker_count, count)
    }
    count_df <- data.frame(matrix(marker_count, ncol=length(marker_count), nrow = 1))
    colnames(count_df) <- markers

    total_cells <- nrow(formatted_data)

    #ITERATION LOOP here...##############
    for (iter_num in 1:num_iter) {
        #bootstrap_df is used in the randomization of markers to generate random phenotypes
        bootstrap_df <- data.frame(matrix(0, nrow = nrow(formatted_data), ncol = length(markers)))
        rownames(bootstrap_df) <- formatted_data$Cell.ID
        colnames(bootstrap_df) <- markers

        #randomly assign marker expressions
        for (marker in markers) {
            marker_count <- count_df[,marker]
            rand <- sample(total_cells, size = marker_count)
            bootstrap_df[rand,marker] <- 1
        }

        #start a new column for phenotype
        bootstrap_df$Phenotype <- ""
        for (i in 1:length(markers)) {
            marker <- paste(markers[i], ",", sep="")
            #select the marker column that was assigned to be true (1) for the marker and add the marker name as phenotype
            bootstrap_df[bootstrap_df[, i] == 1, ]$Phenotype <- paste(bootstrap_df[bootstrap_df[, i] == 1, ]$Phenotype, marker, sep="")
        }

        #get rid of comma at the end
        bootstrap_df$Phenotype <- gsub(",$", "", bootstrap_df$Phenotype)
        #add "OTHER" as phenotype for those without a phenotype, since they're all DAPI positive
        bootstrap_df[bootstrap_df$Phenotype == "", ]$Phenotype <- "OTHER"

        #get all unique phenotypes generated
        phenotypes_generated <- unique(bootstrap_df$Phenotype)
        #count the phenotypes and read it into result df
        for (phenotype in phenotypes_generated) {
            count <- nrow(bootstrap_df[bootstrap_df$Phenotype == phenotype, ])
            results[phenotype, iter_num] <- count
        }

    }

    #start a summary dataframe
    summary_df <- data.frame(matrix(nrow = length(marker_combinations), ncol=5))
    rownames(summary_df) <- marker_combinations
    colnames(summary_df) <- c("Percentage_of_occurrence", "Observed_cell_number", "Average_bootstrap_cell_number",
                              "Enrichment.p", "Depletion.p")


    #calculate the percentage, enrichment, depletion, observed_cell_number, average_bootstrap_cell_number
    for (combination in marker_combinations) {

        #percentage
        combination_scores <- results[combination, ]
        num_non_zero <- length(combination_scores[combination_scores != 0])
        percentage_presence <- num_non_zero/num_iter * 100
        summary_df[combination, "Percentage_of_occurrence"] <- percentage_presence

        #observed_cell_number
        num_observed <- nrow(formatted_data[formatted_data$Phenotype == combination, ])
        summary_df[combination, "Observed_cell_number"] <- num_observed

        #average_bootstrap_cell_number
        average_bootstrap_cell_number <- mean(unlist(combination_scores))
        summary_df[combination, "Average_bootstrap_cell_number"] <- average_bootstrap_cell_number

        #enrichment
        num_greater <- sum(num_observed > combination_scores)
        if (length(num_greater) == 0) {
            enrichment_score <- 1
        } else {
            enrichment_score <- 1-(num_greater/num_iter)
            if(enrichment_score == 0){
              enrichment_score <- 1/num_iter
            }
        }
        summary_df[combination, "Enrichment.p"] <- enrichment_score

        #depletion
        num_lesser <- sum(num_observed < combination_scores)
        if (length(num_lesser) == 0) {
            depletion_score <- 1
        } else {
            depletion_score <- 1-(num_lesser/num_iter)
            if(depletion_score == 0){
              depletion_score <- 1/num_iter
            }
        }
        summary_df[combination, "Depletion.p"] <- depletion_score
    }

    return(summary_df)

}
