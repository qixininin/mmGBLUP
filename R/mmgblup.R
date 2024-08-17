#' mmgblup function
#'
#' @param data a data frame
#' @param Ka additive genetic relationship matrix
#' @param Kd dominance genetic relationship matrix
#'
#' @return list(mod, BV)
#' @export
#' @import sommer
#' @import dplyr
#'
#' @examples \dontrun{rst = mmgblup(data = cbind(dt, mmdata$Xa, mmdata$Xd), Ka = mmdata$Ka, Kd = mmdata$Ka)}
mmgblup <- function(data, Ka, Kd)
{
  ## Evaluate input
  if(missing(Ka) | missing(Kd)){
    stop("Error: no input for additive kinship matrix or dominance covariance matrix")
  }

  ## Receive fixed effect column names
  ## If no fixed effect columns in data, then only intercept is fixed effect
  fixColName = names(data)[!(names(data) %in% c("GID","trait"))]
  if(length(fixColName)==0){
    fixColName = "1"
  }

  datafake = data %>% dplyr::mutate(GID1=GID)

  mod = mmer(reformulate(fixColName, "trait"),
             random = ~vsr(GID, Gu=Ka) + vsr(GID1, Gu=Kd),
             rcov = ~units,
             data = datafake,
             verbose = FALSE, date.warning = FALSE)

  BV = subset.data.frame(data, select = "GID")
  fixColName_remain = as.vector(mod$Beta$Effect)
  fixColName_remain_Aidx = grep(pattern = "_A$", x = fixColName_remain)
  fixColName_remain_Didx = grep(pattern = "_D$", x = fixColName_remain)

  BV$mu  = mod$Beta$Estimate[1]                                                                     # mu
  BV$A_l = as.matrix(data[,fixColName_remain[fixColName_remain_Aidx]]) %*%
    as.vector(mod$Beta$Estimate)[fixColName_remain_Aidx]                                            # A-major
  BV$D_l = as.matrix(data[,fixColName_remain[fixColName_remain_Didx]]) %*%
    as.vector(mod$Beta$Estimate)[fixColName_remain_Didx]                                            # D-major
  BV$A_s = mod$U$`u:GID`$trait[BV$GID]                                                              # A-minor
  BV$D_s = mod$U$`u:GID1`$trait[BV$GID]                                                             # D-minor
  BV$pre = BV$mu + BV$A_l + BV$D_l + BV$A_s + BV$D_s

  return(list(mod, BV))

}
