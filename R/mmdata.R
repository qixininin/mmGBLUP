#' mmdata function
#' Format data for mmGEBLUP analysis
#'
#' @param geno_data input genotype data
#' @param pheno_data input phenotype data
#' @param qtl_data input additice qtl data
#' @param qtl_dom_data input dominance qtl data
#'
#' @return a list with QCed genotype and phenotype data, and other data required for GS analysis
#'         $ mmgeno_data an n*m matrix
#'         $ mmpheno_data an n_obs*3 data frame
#'         $ Ka additive relationship matrix for mmGBLUP
#'         $ Xa additive fixed effect coefficient for mmGBLUP
#'         $ Kd dominance relationship matrix for mmGBLUP
#'         $ Xd dominance fixed effect coefficient for mmGBLUP
#'         $ A additive relationship matrix for GBLUP
#'         $ D additive relatinoship matrix for GBLUP
#'         $ mmsummary data summary
#' @export
#' @import dplyr
#'
#' @examples \dontrun{qcdata = dataqc(geno_data, pheno_data, qtl_data, qtl_env_data)}
mmdata <- function(geno_data, pheno_data, qtl_data, qtl_env_data)
{
  # Transform genotype
  rownames(geno_data) = geno_data[,2]
  mmgeno_data = t(geno_data[,-c(1:3)])

  # Summary
  traitName = colnames(pheno_data)[2]
  lineName = intersect(unique(pheno_data$GID),rownames(mmgeno_data))
  lineNum = length(lineName)

  # Prepare genotype
  mmgeno_data = as.matrix(mmgeno_data[lineName, ])

  # Prepare phenotype
  mmpheno_data = pheno_data %>%
    dplyr::filter(GID %in% lineName) %>%
    droplevels()
  colnames(mmpheno_data)[colnames(mmpheno_data) == traitName] <- deparse(substitute(trait))

  # Prepare qtl
  if(!is.null(qtl_data)){
    site_qtl = qtl_data$QTL
    m_qtl = length(site_qtl)
  } else {
    m_qtl = 0
  }

  # Prepare qtl_dom
  if(!is.null(qtl_dom_data)){
    site_qtl_dom = qtl_dom_data$QTL
    m_qtl_dom = length(site_qtl_dom)
  } else {
    m_qtl_dom = 0
  }

  # Additive relationship matrix
  A = A.mat(mmgeno_data)
  colnames(A) = rownames(A) = rownames(mmgeno_data)

  # Dominance relationship matrix
  D = D.mat(mmgeno_data)
  colnames(D) = rownames(D) = rownames(mmgeno_data)

  # Calculate reduced additive relationship matrix
  KaXa = calculateKaXa(mmgeno_data, mmpheno_data, A, site_qtl, m_qtl)

  # Calculate reduced dominance relationship matrix
  KdXd = calculateKdXd(mmgeno_data, mmpheno_data, D, site_qtl_dom, m_qtl_dom)

  return(list(mmgeno_data = mmgeno_data,
              mmpheno_data = mmpheno_data,
              Ka = KaXa$Ka,
              Xa = KaXa$Xa,
              Kd = KdXd$Kd,
              Xd = KdXd$Xd,
              A = A,
              D = D,
              mmsummary = list(traitName = traitName,
                               lineName  = lineName)))
}
