#' qtxnetwork.output.trans
#'
#' @param pheno_data a data frame of phenotypic data, one of the outputs from simulation
#' @param pre_file path and name for .pre file
#'
#' @return a list
#'         $ qtl_data a data frame with additive qtl information
#'         $ qtl_dom_data a data frame with dominance qtl information
#' @export
#'
#' @examples \dontrun{qtl = qtxnetwork.output.trans(pheno_data, pre_file)
#'                    qtl$qtl_data
#'                    qtl$qtl_dom_data}
qtxnetwork.output.trans <- function(pheno_data, pre_file)
{
  traitName = colnames(pheno_data)[-c(1)]
  traitNum = length(traitName)

  start_title = "_1D_effect"
  end_title = "_1D_heritability"
  df.qtl = data.frame()
  for(c in 1:traitNum)
  {
    trait = traitName[c]
    file_lines = readLines(pre_file)
    # locate target lines
    start_line <- 0
    end_line <- 0
    for (i in 1:length(file_lines)) {
      line <- file_lines[i]
      if (line == start_title) {
        start_line <- i
        next
      }
      if (line == end_title) {
        end_line <- i
        break
      }
    }
    # store in data.frame
    if(start_line) {
      tmp <- as.data.frame(do.call(rbind, strsplit(file_lines[(start_line + 2):(end_line - 2)], "\\s+")))
      df.qtl = rbind(df.qtl, data.frame(trait, tmp))
    }
  }
  colnames(df.qtl) = c("TRAIT", "QTL", "SNPID", "A", "SE", "P-Value","D", "DSE", "DP-Value")

  dt = df.qtl %>% dplyr::select(c("TRAIT","SNPID","A", "D"))
  dt_reshape=reshape(dt,
                     idvar=c("TRAIT", "SNPID"),
                     varying=c("A","D"),
                     v.names="Effect",
                     timevar="Type",
                     times=c("A","D"),
                     direction="long")
  # Additive qtl
  qtl_data = dt_reshape %>%
    dplyr::filter(Type == "A") %>%
    dplyr::filter(Effect != "---") %>%
    dplyr::select(c("TRAIT", "SNPID")) %>%
    dplyr::rename(QTL = SNPID)
  # AE qtl
  qtl_dom_data = dt_reshape %>%
    dplyr::filter(Type == "D") %>%
    dplyr::filter(Effect != "---") %>%
    dplyr::select(c("TRAIT", "SNPID")) %>%
    dplyr::rename(QTL = SNPID)

  if(nrow(qtl_data)>0) rownames(qtl_data) = 1:nrow(qtl_data)
  if(nrow(qtl_dom_data)>0) rownames(qtl_dom_data) = 1:nrow(qtl_dom_data)

  return(list(qtl_data = qtl_data, qtl_dom_data = qtl_dom_data))
}
