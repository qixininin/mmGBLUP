#' calculateKdXd function
#'
#' @param geno_data genotype data
#' @param pheno_data phenotype data
#' @param D dominance relationship matrix
#' @param site_qtl_dom a vector, names of dom QTL sites
#' @param m_qtl_dom the number of dom QTLs
#'
#' @return list(Kd = Kd,Xd = Xd)
#' @export
#' @import dplyr
#'
#' @examples \dontrun{KdXd = calculateKdXd(mmgeno_data, mmpheno_data, D, site_qtl, m_qtl)}
calculateKdXd <- function(geno_data, pheno_data, D, site_qtl_dom, m_qtl_dom)
{
  # extract qtl information from what we got from GWAS
  if(m_qtl_dom>0)
  {
    Kd = D.mat(as.matrix(geno_data[,!colnames(geno_data) %in% site_qtl_dom]))
    Xd = 1-abs(geno_data[as.character(pheno_data$GID), site_qtl_dom]) # transform -1,0,1 to 0,1,0 code
    Xd = as.matrix(Xd)
    colnames(Xd) = paste0(site_qtl_dom,"_D")
  } else {
    Kd = D
    Xd = c()
  }

  colnames(Kd) = rownames(Kd) = rownames(geno_data)

  return(list(Kd = Kd,Xd = Xd))
}
