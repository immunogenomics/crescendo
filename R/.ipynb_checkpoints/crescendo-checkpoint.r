#' Implementation of the Crescendo algorithm.
#'
#' Runs the Crescendo algorithm by performing the downsampling (if necessary), estiimation, marginalization, and matching steps.

#' @inheritParams get_downsampled_idx
#' @inheritParams estimate_betas
#' @inheritParams marginalize_effects
#' @inheritParams match_cdfs
#' @param impute_reads If not NULL, supplies a vector of UMI to perform imputation at the supplied read depths.
#' @param return_obj Return only corrected counts, or results from each step.
#' @import tidyverse
#' @import magrittr
#' @import parallel
#' @import data.table
#' @import Matrix
#' @import glmnet
#' @importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importFrom stats coef median model.matrix sd
#' @return A dataframe containing corrected gene count expressions, or a list containing results from each Crescendo step.
#' @export
crescendo <- function(Ycounts, meta_data, R, batch_var = 'batch', genes_use = NULL, 
                      prop = 1, min_cells = 50, lambda = NULL, alpha = 0, seed = 2, 
                      mc.cores = NULL, impute_reads = NULL, return_obj = FALSE, verbose = FALSE){
    if(!identical(rownames(meta_data), colnames(Ycounts))){
        stop('Metadata does not have rownames equal to the colnames of the expression matrix. Please ensure these are identical.')
    }
    if(!(batch_var %in% colnames(meta_data))){
        stop('Provided batch variable does not exist. Please input name of batch column.')
    }
    if(is.null(mc.cores)){
        mc.cores <- max(1, detectCores() - 1)
    }
    meta_data$batch <- meta_data[[batch_var]]

    # Merge redundant clusters
    R <- t(merge_redundant_clusters(R, 0.99))
    colnames(R) <- paste0('R', 1:ncol(R))

    # Use selected genes if applicable
    if(!is.null(genes_use)){
        Ycounts <- Ycounts[genes_use, , drop = FALSE]
    } else{
        genes_use <- rownames(Ycounts)
    }
    Ycounts <- t(Ycounts)

    offset <- log(meta_data$nUMI)
    # If custom read depths are provided for "impute_reads", use these as the final reads for imputation
    # instead of a constant read depth equal to the median read depth
    if(is.null(impute_reads)){
        final_reads <- log(rep(median(meta_data$nUMI), length(meta_data$nUMI)))
    } else{
        final_reads <- impute_reads
    }

    if(verbose){
        message('Creating design matrix')
    }
    design_full <- create_design_matrix(meta_data = meta_data, R = R)
    # Designate terms to keep (e.g. cluster terms)
    terms_keep <- grep('^R\\d+$|disease|(Intercept)', colnames(design_full), value = TRUE)
    # Designate terms to remove (e.g. interaction terms)
    terms_remove_crossed <- terms_remove_crossed <- grep(':batch', colnames(design_full), value = TRUE)
    
    # Downsampling
    if(prop != 1){
        idx <- get_downsampled_idx(meta_data = meta_data, R = R, mc.cores = 20, verbose = verbose)
    } else{
        if(verbose){
            message('No downsampling specified')
        }
        idx <- seq(nrow(meta_data))
    }
    
    # Estimation
    if(verbose){
        message('Estimating')
    }
    betas <- estimate_betas(
        Ycounts = Ycounts[idx, ],
        R = R[idx, ],
        design_full = design_full[idx, ],
        lambda = lambda,
        alpha = alpha,
        genes_use = genes_use,
        offset = offset[idx],
        mc.cores = mc.cores
    )
    
    # Marginalization
    if(verbose){
        message('Marginalizing')
    }
    expected_counts_marginalized <- marginalize_effects(
        R = R, 
        design_full = design_full, 
        betas = betas, 
        terms_keep = terms_keep, 
        terms_remove_crossed = terms_remove_crossed, 
        offset = final_reads
    )
    
    # Matching
    if(verbose){
        message('Matching')
    }
    corrected_matching <- match_cdfs(
        expected_counts_marginalized = expected_counts_marginalized, 
        Ycounts = Ycounts, 
        design_full = design_full, 
        betas = betas, 
        offset = final_reads
    )
    
    if (return_obj) {
        list(
            betas = betas, # from estimation
            expected_counts_marginalized = expected_counts_marginalized, # from marginalization
            corrected_matching = corrected_matching # from matching
        )
    } else {
        return(corrected_matching)            
    }
}