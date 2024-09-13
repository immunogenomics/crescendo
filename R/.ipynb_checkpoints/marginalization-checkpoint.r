#' Remove effects of confounding factors from a Poisson generalized linear model.
#'
#' Remove effects of confounding factors from a Poisson generalized linear model (obtained from the estimation step) to obtain a 
#' marginalized model.
#'
#' @param R A matrix containing cluster assignments for each cell (cluster x cell). 
#' @param design_full A design matrix that contains batch and cluster membership information.
#' @param betas A dataframe containing parameter estimates for each gene.
#' @param terms_keep A vector containing the terms that we want to keep in the marginalized model.
#' @param terms_remove_crossed A vector containing the terms we want to marginalize out (e.g. crossed terms between a confounding
#' factor and clusters).
#' @param offset A vector containing the nUMI to sample up to.
#' @return A dataframe containing rates from the marginalized model for each cell.
#' @importFrom purrr map map_chr 
#' @importFrom stringr str_split
#' @export
marginalize_effects <- function(
    R, 
    design_full, 
    betas, 
    terms_keep, 
    terms_remove_crossed, 
    offset
){
    expected_logcounts_marginalized <- solve_mu1(
        design_full[, terms_keep, drop = FALSE],
        betas[terms_keep, , drop = FALSE]
    )
    colnames(expected_logcounts_marginalized) <- colnames(betas)    
    if(ncol(R) == 1){
        mu_list <- split(terms_remove_crossed, map_chr(str_split(terms_remove_crossed, ':'), 1)) %>% 
                map(function(terms) {
                    colMeans(betas[terms, , drop = FALSE])
                })
        mu <- Reduce(Matrix::rbind2, mu_list) %>% t %>% data.matrix
        rownames(mu) <- colnames(R)

        sigmas_list <- split(terms_remove_crossed, map_chr(str_split(terms_remove_crossed, ':'), 1)) %>% 
            map(function(terms) {
                apply(betas[terms, , drop = FALSE], 2, sd)
            })
        sigmas <- Reduce(Matrix::rbind2, sigmas_list) %>% t %>% data.matrix
        rownames(sigmas) <- colnames(R)
    } else{
        mu_list <- split(terms_remove_crossed, map_chr(str_split(terms_remove_crossed, ':'), 1)) %>% 
            map(function(terms) {
                colMeans(betas[terms, , drop = FALSE])
            })
        mu <- Reduce(Matrix::rbind2, mu_list) 
        rownames(mu) <- names(mu_list)
        mu <- mu[colnames(R), , drop = FALSE]

        sigmas_list <- split(terms_remove_crossed, map_chr(str_split(terms_remove_crossed, ':'), 1)) %>% 
            map(function(terms) {
                apply(betas[terms, , drop = FALSE], 2, sd)
            })
        sigmas <- Reduce(Matrix::rbind2, sigmas_list) 
        rownames(sigmas) <- names(sigmas_list)
        sigmas <- sigmas[colnames(R), , drop = FALSE]
    }
    
    solve_mu2(expected_logcounts_marginalized, R, sigmas, mu) ## in memory 

    expected_counts_marginalized <- exp(expected_logcounts_marginalized + offset)
    rownames(expected_counts_marginalized) <- rownames(design_full)
    return(expected_counts_marginalized)
}