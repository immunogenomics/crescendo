#' Obtain corrected counts.
#'
#' Obtain corrected counts by performing inverse cumulative distribution matching.
#'
#' @param expected_counts_marginalized A dataframe containing rates from the marginalized model for each cell.
#' @param Ycounts A matrix containing un-normalized gene count information (cells x genes).
#' @param design_full A design matrix that contains batch and cluster membership information.
#' @param betas A dataframe containing parameter estimates for each gene.
#' @param offset A vector containing the nUMI to sample up to.
#' @return A dataframe containing corrected gene expression counts.
#' @export
match_cdfs <- function(
    expected_counts_marginalized, 
    Ycounts, 
    design_full, 
    betas, 
    offset
){
    expected_counts_full <- exp(design_full %*% betas + offset) %>% data.matrix
    if(!('matrix' %in% class(Ycounts))){
        Ycounts <- Ycounts %>% data.matrix
    }
    
    corrected_matching = as(parallelicdfPoisson(
        Ycounts,
        expected_counts_full, 
        expected_counts_marginalized
    ), 'dgCMatrix')
    
    colnames(corrected_matching) <- colnames(Ycounts)
    rownames(corrected_matching) <- rownames(Ycounts) 
    return(corrected_matching)
}