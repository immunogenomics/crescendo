#' Create a design matrix.
#'
#' This creates a desgin matrix from a meta-data table and cluster identities.
#'
#' @param meta_data Metadata table that contains cell information (cell x info).
#' @param R A matrix containing cluster assignments for each cell (cluster x cell). 
#' @return A design matrix that contains batch and cluster membership information.
#' @export
create_design_matrix <- function(meta_data, R){
    design_full <- cbind(
        stats::model.matrix(~1, meta_data),
        stats::model.matrix(~0 + R),
        stats::model.matrix(~0 + R:batch, meta_data),
        stats::model.matrix(~0 + batch, meta_data)
    )
    
    colnames(design_full) <- gsub('^R?(R.*)', '\\1', colnames(design_full))
    rownames(design_full) <- rownames(meta_data)
    design_full <- as(design_full, 'dgCMatrix')
    return(design_full)
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
