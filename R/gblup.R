#' gblup function
#'
#' @param data a data frame
#' @param A additive genetic relationship matrix
#' @param D dominance genetic relationship matrix
#'
#' @return list(mod, BV)
#' @export
#' @import sommer
#' @import dplyr
#'
#' @examples \dontrun{rst = gblup(model = "AD", data = dt, A = mmdata$A, D = mmdata$D)}
gblup <- function(model, data, A, D)
{
  ## Evaluate input
  if (missing(A)) {
    stop("Error: no input for additive kinship matrix")
  }

  if (model == "A") {
    mod = mmer(reformulate("1", "trait"),
               random = ~vsr(GID, Gu = A),
               rcov = ~units,
               data = data,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    BV = subset.data.frame(data, select = "GID")
    BV$mu = mod$Beta$Estimate[1]  # mu
    BV$A = mod$U$`u:GID`$trait[BV$GID]  # A
    BV$pre = BV$mu + BV$A
  } else if (model == "AD") {
    if (missing(D)) {
      stop("Error: no input for dominance kinship matrix")
    }

    datafake = data %>% dplyr::mutate(GID1=GID)

    mod = mmer(reformulate("1", "trait"),
               random = ~vsr(GID, Gu = A) + vsr(GID1, Gu = D),
               rcov = ~units,
               data = datafake,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    BV = subset.data.frame(data, select = "GID")
    BV$mu = mod$Beta$Estimate[1]  # mu
    BV$A = mod$U$`u:GID`$trait[BV$GID]  # A
    BV$D = mod$U$`u:GID1`$trait[BV$GID]  # D
    BV$pre = BV$mu + BV$A + BV$D
  } else {
    stop("Error: Model type not recognized. Please provide either 'A' or 'AD'.")
  }

  return(list(mod, BV))
}
