#' snp.effect function
#'
#' @param model a character that indicate the genetic effects included ("A" or "AD")
#' @param snpNum The number of SNPs
#' @param major_a_idx The index for major additive SNP
#' @param major_d_idx The index for major dominance SNP
#' @param variance_a_major The variance for major additive SNP effect
#' @param variance_a_minor The variance for minor additive SNP effect
#' @param variance_d_major The variance for major dominance SNP effect
#' @param variance_d_minor The variance for minor dominance SNP effect
#'
#' @return effect an M*p matrix p is the number of genetic effects included (such as "A" for 1, "AD" for 2)
#'
#' @export
#' @importFrom stats rnorm
#'
#' @examples snp.effect(model = "AD", snpNum = 2000,
#'                      major_a_idx = c(500, 750, 1000, 1250, 1500), variance_a_major = 0.02, variance_a_minor = 0.002,
#'                      major_d_idx = c(500, 750, 1000, 1250, 1500), variance_d_major = 0.02, variance_d_minor = 0.002)
snp.effect <- function(model, snpNum,
                       major_a_idx, variance_a_major, variance_a_minor,
                       major_d_idx, variance_d_major, variance_d_minor) {

  if(model == "A") {
    if(max(major_a_idx) > snpNum) {
      stop("Error: snp.effect(). The input major_a_idx is out of the snpNum range.")
    }

    if(missing(variance_a_major) | missing(variance_a_minor)) {
      stop("Error: snp.effect(). Please provide variance values for both major and minor additive effects.")
    }

    # determine minor index for additive effects
    minor_a_idx = setdiff(1:snpNum, major_a_idx)

    # store main effect for all SNPs
    effects = data.frame(A = rep(0, snpNum))

    # Simulate additive genetic effect
    effects[major_a_idx, "A"] <- rnorm(length(major_a_idx), mean = 0, sd = sqrt(variance_a_major))
    effects[minor_a_idx, "A"] <- rnorm(length(minor_a_idx), mean = 0, sd = sqrt(variance_a_minor))
  } else if(model == "AD") {
    if(max(major_a_idx) > snpNum | max(major_d_idx) > snpNum) {
      stop("Error: snp.effect(). The input major_a_idx or major_d_idx is out of the snpNum range.")
    }

    if(missing(variance_d_major) | missing(variance_d_minor)) {
      stop("Error: snp.effect(). Please provide variance values for both major and minor dominance effects.")
    }

    # determine minor index for additive and dominance effects
    minor_a_idx = setdiff(1:snpNum, major_a_idx)
    minor_d_idx = setdiff(1:snpNum, major_d_idx)

    # store main and interaction effect for all SNPs
    effects = data.frame(A = rep(0, snpNum), D = rep(0, snpNum))

    # Simulate genetic effects
    effects[major_a_idx, "A"] <- rnorm(length(major_a_idx), mean = 0, sd = sqrt(variance_a_major))
    effects[minor_a_idx, "A"] <- rnorm(length(minor_a_idx), mean = 0, sd = sqrt(variance_a_minor))

    effects[major_d_idx, "D"] <- rnorm(length(major_d_idx), mean = 0, sd = sqrt(variance_d_major))
    effects[minor_d_idx, "D"] <- rnorm(length(minor_d_idx), mean = 0, sd = sqrt(variance_d_minor))
  } else {
    stop("Error: snp.effect(). Unsupported model type. Please use 'A' for additive effects only or 'AD' for both additive and dominance effects.")
  }

  # Return effects matrix
  return(effects)
}
