#' Estimate gene structue with a Poisson generalized linear model.
#'
#' Fit a gene's expression distribution with a Poisson generalized linear model. Default is ridge penalty.
#'
#' @param Ycounts A matrix containing un-normalized gene count information (cells x genes).
#' @param R A matrix containing cluster assignments for each cell (cluster x cell). 
#' @param design_full A design matrix that contains batch and cluster membership information.
#' @param lambda Regularization parameter.
#' @param alpha Elastic net regularization parameter. 1 indicates lasso penalty, 0 indicates ridge penalty.
#' @param genes_use Genes to estimate. If none supplied, will automatically estimate all genes in Ycounts.
#' @param offset A vector containing the total nUMI for each cell.
#' @param mc.cores Number of cores to use for parallel functions.
#' @return A dataframe containing parameter estimates for each gene.
#' @export
estimate_betas <- function(
    Ycounts, ## cells x genes
    R, 
    design_full, 
    lambda = NULL,
    alpha = 0,
    genes_use = NULL,
    offset = NULL,
    mc.cores = NULL
){
    if(is.null(mc.cores)){
        mc.cores <- max(1, detectCores() - 1)
    }
    if(is.null(genes_use)){
        genes_use <- colnames(Ycounts)
    }
    
    betas_list <- mclapply(genes_use, function(i){
        y <- Ycounts[, i]
        x <- design_full[, -1]
        fit <- glmnet(
            x = x,
            y = y,
            family = 'poisson',
            lambda = lambda,
            offset = offset,
            alpha = alpha
        )
        return(coef(fit, s = 1)[, 1])
    }, mc.preschedule = TRUE, mc.cores = mc.cores)
    
    names(betas_list) <- genes_use
    betas <- betas_list %>% data.frame %>% data.matrix
    return(betas)
}