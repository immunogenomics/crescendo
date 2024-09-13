#' Sample a discrete harmony cluster based off of soft-cluster probabilities.
#'
#' Harmony soft-clusters give probabilities for each cell. In order to downsample data, we sample a discrete cluster for each cell
#' based on its harmony cluster probabilities.
#'
#' @param R A matrix of harmony clusters (cluster x cell) that are received as output from the Harmony algorithm.
#' @param seed Seed for RNG.
#' @param mc.cores Number of cores to use for parallel functions.
#' @return A vector containing a discrete cluster assignment for each cell design matrix that contains batch and cluster membership 
#' information.
#' @importFrom tidyselect all_of
#' @export
sample_harmony_clus <- function(R, seed = NULL, mc.cores = NULL){
    if(is.null(mc.cores)){
        mc.cores <- max(1, detectCores() - 1)
    }
    if(!is.null(seed)){
        RNGkind("L'Ecuyer-CMRG")
        set.seed(seed)
    }
    sampled_clus <- factor(mclapply(seq(nrow(R)), function(x){
        sample(x = colnames(R), size = 1, prob = R[x, ])
    }, mc.cores = mc.cores, mc.preschedule = TRUE), levels = colnames(R))
    return(sampled_clus)
}

#' Downsample data in a batch and cell-type aware manner.
#'
#' To downsample data, we ensure that each cell-type in each batch is represented by a minimum number of cells (or if lower than
#' the minimum, all of the cells).
#'
#' @param meta_data Metadata table that contains cell information (cell x info).
#' @param batch_var Name of column containing batch information (e.g. batch, dataset, technology).
#' @param prop Fraction of cells to downsample to (e.g. 0.5 means downsample to arouond 50\% of original cells, 0.25 means 25\%, etc.).
#' @param min_cells Number of minimum cells to keep from each cell type in each batch.
#' @param return_idx Return just the indices of the downsampled metadata, or return the metadata table.
#' @param clus_col Name of column containing cluster information.
#' @param seed Seed for RNG.
#' @param mc.cores Number of cores to use for parallel functions.
#' @return A metadata table containing the downsampled cells.
#' @importFrom dplyr slice_sample
#' @export
downsample_batch_clus <- function(
    meta_data, 
    batch_var = 'batch', 
    prop = 0.1,
    min_cells = 50,
    return_idx = TRUE,
    clus_col = 'sampled_clus', 
    seed = NULL,
    mc.cores = NULL
){
    col_get <- c(batch_var, clus_col)
    temp_meta <- meta_data %>% dplyr::select(all_of(col_get))
    temp_meta$idx <- seq(nrow(temp_meta))
    temp_meta <- split(temp_meta, temp_meta[[batch_var]])
    temp_meta <- lapply(temp_meta, function(x){
        split(x, x[[clus_col]])
    })
    
    if(!is.null(seed)){
        RNGkind("L'Ecuyer-CMRG")
        set.seed(seed)
    }
    temp_meta <- mclapply(temp_meta, function(x){
        lapply(x, function(y){
            if(nrow(y) > min_cells){
                samp_n <- max(min_cells, ceiling(nrow(y) * prop))
                new_frame <- y %>% dplyr::slice_sample(n = samp_n)
                return(new_frame)
            } else{
                return(y)
            }
        })
    }, mc.cores = mc.cores)
    
    temp_meta <- unlist(temp_meta, recursive = FALSE)
    temp_meta <- rbindlist(temp_meta)
    if(return_idx){
        return(temp_meta$idx)
    } else{
        return(temp_meta)
    }
}

#' Get the indices from downsampled data
#'
#' To downsample data, we ensure that each cell-type in each batch is represented by a minimum number of cells (or if lower than
#' the minimum, all of the cells). This function can return just the indices, or the whole medatadata table.
#'
#' @inheritParams sample_harmony_clus 
#' @inheritParams downsample_batch_clus
#' @param verbose Print out progress messages.
#' @return Indices of downsample data (or the whole metadata table if return_idx = FALSE)
#' @export
get_downsampled_idx <- function(
    meta_data,
    R, 
    batch_var = 'batch', 
    prop = 0.1,
    min_cells = 50, 
    return_idx = TRUE,
    clus_col = 'sampled_clus',
    seed = 2, 
    mc.cores = NULL,
    verbose = FALSE
){
    RNGkind("L'Ecuyer-CMRG")
    if(!is.null(seed)){
        set.seed(seed)
    }
    if(verbose){
        message("Getting discrete Harmony cluster")
    }
    meta_data$sampled_clus <- sample_harmony_clus(R = R, mc.cores = mc.cores, seed = NULL)
    if(verbose){
        message("Downsampling")
    }
    idx <- downsample_batch_clus(
        meta_data = meta_data, batch_var = batch_var, clus_col = "sampled_clus", 
        min_cells = min_cells, prop = prop, return_idx = return_idx, mc.cores = mc.cores, seed = NULL
    )
    return(idx)
}