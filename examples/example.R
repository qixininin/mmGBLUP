library(mmGBLUP)
set.seed(421)

# Simulation ----
## (A) additive  effects ----
mod = "A"
prefix = paste0("./inst/simulation-", mod)
indNum = 1000
snpNum = 2000
rho_a = 0.8
h2_a = 0.6
sigma_a = h2_a
sigma_error = 1-h2_a

major_a_idx = c(500, 750, 1000, 1250, 1500)
snpNum_a_major = length(major_a_idx)
snpNum_a_minor = snpNum-snpNum_a_major

sigma_a_major = rho_a * sigma_a/ snpNum_a_major
sigma_a_minor = (1-rho_a) * sigma_a / snpNum_a_minor

# Effect generation
effects = snp.effect(model = mod, snpNum = snpNum,
                     major_a_idx = major_a_idx, variance_a_major = sigma_a_major, variance_a_minor = sigma_a_minor,)

# Genotype generation
geno_data = geno.generate(indNum = indNum, snpNum = snpNum, maf.min = 0.05, maf.max = 0.5,
                          chr.snpNum = c(500, 500, 500, 500))

# Phenotype generation
pheno_data = pheno.generate(model = mod, geno_data = geno_data, effects = effects, sigma.error = sigma_error)

# Save Rdata
save(geno_data, pheno_data, file = paste0(prefix, "-genphe.Rdata"))

## (AD) additive and dominance effects ----
mod = "AD"
prefix = paste0("./inst/simulation-", mod)
indNum = 3000
snpNum = 2000
rho_a = 0.5
rho_d = 0.5
h2_a = 0.3
h2_d = 0.2
sigma_a = h2_a
sigma_d = h2_d
sigma_error = 1-h2_a-h2_d

major_a_idx = c(500, 750, 1000, 1250, 1500)
major_d_idx = c(200, 750, 1250, 1800)
snpNum_a_major = length(major_a_idx)
snpNum_a_minor = snpNum-snpNum_a_major
snpNum_d_major = length(major_d_idx)
snpNum_d_minor = snpNum-snpNum_d_major

sigma_a_major = rho_a * sigma_a/ snpNum_a_major
sigma_a_minor = (1-rho_a) * sigma_a / snpNum_a_minor
sigma_d_major = rho_d * sigma_d/ snpNum_d_major
sigma_d_minor = (1-rho_d) * sigma_d / snpNum_d_minor

# Effect generation
effects = snp.effect(model = mod, snpNum = snpNum,
                     major_a_idx = major_a_idx, variance_a_major = sigma_a_major, variance_a_minor = sigma_a_minor,
                     major_d_idx = major_d_idx, variance_d_major = sigma_d_major, variance_d_minor = sigma_d_minor)

# Genotype generation
geno_data = geno.generate(indNum = indNum, snpNum = snpNum, maf.min = 0.05, maf.max = 0.5,
                          chr.snpNum = c(500, 500, 500, 500))


# Phenotype generation
pheno_data = pheno.generate(model = mod, geno_data = geno_data, effects = effects, sigma.error = sigma_error)

# Save Rdata
save(geno_data, pheno_data, file = paste0(prefix, "-genphe.Rdata"))
# save(effects, major_a_idx, major_d_idx, file = paste0(prefix, "-eff.Rdata"))

# QTS ----
prefix = "./inst/simulation-AD"
load(paste0(prefix, "-genphe.Rdata"))

# Prepare QTXNetwork
qtxnetwork.input.trans(geno_data, pheno_data,
                       geno_output_prefix = prefix,
                       pheno_output_prefix = prefix)

# Perform QTXNetwork
qtxnetwork_path = "/public3/zqx/QTXNetwork_4.0/build/QTXNetwork"
gen_file = paste0(prefix, ".gen")
phe_file = paste0(prefix, "_SimTrait.phe")
pre_file = paste0(prefix, "_SimTrait.pre")

qtxnetwork.perform(qtxnetwork_path, gen_file, phe_file, pre_file)

# Extract QTXNetwork
qtl = qtxnetwork.output.trans(pheno_data, pre_file)
qtl_data = qtl$qtl_data
qtl_dom_data = qtl$qtl_dom_data

save(qtl_data, qtl_dom_data, file = paste0(prefix, "-qtl.Rdata"))

# GS ----
## Load data
prefix = "./inst/simulation-AD"
load(paste0(prefix, "-genphe.Rdata"))
load(paste0(prefix, "-qtl.Rdata"))

# geno_data[1:5,1:10]
# pheno_data[1:5,]
# qtl_data
# qtl_dom_data

## Data QC and Summary
mmdata = mmdata(geno_data, pheno_data, qtl_data, qtl_dom_data)
mmgeno_data = mmdata$mmgeno_data
mmpheno_data = mmdata$mmpheno_data

# mmgeno_data[1:5,1:10]
# mmpheno_data[1:5,]

## Set cross-validation sets
cvNum = 4
naNum = ifelse(length(mmdata$mmsummary$lineName)%%cvNum==0, 0, cvNum-length(mmdata$mmsummary$lineName)%%cvNum)
cvSet = matrix(c(sample(mmdata$mmsummary$lineName), rep(NA, naNum)), nrow = cvNum)

# + GBLUP ------------------------------------------------------
a <- 0
gblup_list <- list()
gblup_ad_list <- list()
for(i in 1:cvNum) # loop for cross validation fold
{
  a <- a + 1

  # set phenotype in validation set and in validation env to be NA
  cv = as.vector(na.omit(unique(cvSet[i,])))
  dt = mmpheno_data %>% dplyr::mutate(trait = ifelse(GID %in% cv, NA, trait))
  dt = as.data.frame(dt)

  # GBLUP model
  rst = gblup(model = "A", data = dt, A = mmdata$A)
  BV = rst[[2]]
  cor <- BV %>% dplyr::mutate(obs = mmpheno_data$trait) %>%
    dplyr::filter(GID %in% cv) %>%
    dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>%
    as.numeric()
  gblup_list[[a]] = data.frame(TRAIT = mmdata$mmsummary$traitName, CV = i, COR = cor, R2 = cor^2)

  # GBLUP-AD model
  rst = gblup(model = "AD", data = dt, A = mmdata$A, D = mmdata$D)
  BV = rst[[2]]
  cor <- BV %>% dplyr::mutate(obs = mmpheno_data$trait) %>%
    dplyr::filter(GID %in% cv) %>%
    dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>%
    as.numeric()
  gblup_ad_list[[a]] = data.frame(TRAIT = mmdata$mmsummary$traitName, CV = i, COR = cor, R2 = cor^2)

  print(a)
}

rstGBLUP <- dplyr::bind_rows(gblup_list)
rstGBLUPAD <- dplyr::bind_rows(gblup_ad_list)

# + mmGBLUP ------------------------------------------------------
a <- 0
mmgblup_ad_list <- list()
for(i in 1:cvNum) # loop for cross validation fold
{
  a <- a + 1

  # set phenotype in validation set and in validation env to be NA
  cv = as.vector(na.omit(unique(cvSet[i,])))
  dt = mmpheno_data %>% dplyr::mutate(trait = ifelse(GID %in% cv, NA, trait))
  dt = as.data.frame(dt)

  # mmGBLUP model
  if(!is.null(mmdata$Xa)) {dt = cbind(dt, mmdata$Xa)}
  if(!is.null(mmdata$Xd)) {dt = cbind(dt, mmdata$Xd)}
  rst = mmgblup(data = dt, Ka = mmdata$Ka, Kd = mmdata$Kd)
  BV = rst[[2]]

  # Calculate correlation
  cor <- BV %>% dplyr::mutate(obs = mmpheno_data$trait) %>%
    dplyr::filter(GID %in% cv) %>%
    dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>%
    as.numeric()

  mmgblup_ad_list[[a]] = data.frame(TRAIT = mmdata$mmsummary$traitName, CV = i, COR = cor, R2 = cor^2)
  print(a)
}

rstmmGBLUPAD <- dplyr::bind_rows(mmgblup_ad_list)

## Save correlation results
rst = data.frame(rbind(rstGBLUP, rstGBLUPAD, rstmmGBLUPAD),
                 MODEL = rep(c("GBLUP","GBLUP-AD","mmGBLUP-AD"), each = cvNum))

save(rst, file = paste0(prefix, "-rst.Rdata"))
