#' pheno.generate function
#'
#' @param model a character that indicate the genetic effects included ("A" or "AD")
#' @param geno_data geno_data given from geno.generate() function
#' @param effects an M*p matrix p is the number of genetic effects included (such as "A" for 1, "AD" for 2)
#' @param sigma.error the variance for residual effects.
#'
#' @return pheno_data a data frame, with two columns $GID, and $SimTrait
#' @export
#'
#' @examples \dontrun{pheno.generate(genotypes = t(geno_data[-c(1:3)]), effects = b,
#'                    indNum = indNum, sigma.error = sigma_error)}
pheno.generate <- function(model, geno_data, effects, sigma.error){

  Ga = t(geno_data[-c(1:3)])
  indNum = nrow(Ga)

  if(model == "AD") {
    Gd = 1 - abs(Ga)
    g = tcrossprod(Ga, t(effects[,"A"])) + tcrossprod(Gd, t(effects[,"D"]))
  } else if(model == "A") {
    g = tcrossprod(Ga, t(effects[,"A"]))
  } else {
    stop("Error: Model type not recognized. Please provide either 'A' or 'AD'.")
  }

  pheno_data = data.frame(GID = as.factor(rownames(Ga)))

  g = as.vector(g)
  error = rnorm(indNum, mean = 0, sd = sqrt(sigma.error))
  pheno_data = cbind(pheno_data, g + error)

  colnames(pheno_data) = c("GID","SimTrait")

  pheno_data$GID = as.factor(pheno_data$GID)

  return(pheno_data)
}
