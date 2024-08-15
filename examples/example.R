set.seed(421)

# Simulation ----
## (AD) additive and dominance effects ----
indNum = 100
snpNum = 200
rho_a = 0.5
rho_d = 0.25
h2_a = 0.3
h2_d = 0.15
sigma_a = h2_a
sigma_d = h2_d
sigma_error = 1-h2_a-h2_d

major_a_idx = c(50, 75, 100, 125, 150)
major_d_idx = c(50, 75, 100, 125, 150)
snpNum_a_major = length(major_a_idx)
snpNum_a_minor = snpNum-snpNum_a_major
snpNum_d_major = length(major_d_idx)
snpNum_d_minor = snpNum-snpNum_d_major

sigma_a_major = rho_a * sigma_a/ snpNum_a_major
sigma_a_minor = (1-rho_a) * sigma_a / snpNum_a_minor
sigma_d_major = rho_d * sigma_d/ snpNum_d_major
sigma_d_minor = (1-rho_d) * sigma_d / snpNum_d_minor

# Effect generation
effects = snp.effect(model = "AD", snpNum = snpNum,
                     major_a_idx = major_a_idx, variance_a_major = sigma_a_major, variance_a_minor = sigma_a_minor,
                     major_d_idx = major_d_idx, variance_d_major = sigma_d_major, variance_d_minor = sigma_d_minor)

# Genotype generation
geno_data = geno.generate(indNum = indNum, snpNum = snpNum, maf.min = 0.05, maf.max = 0.5, chr.snpNum = c(50, 50, 50, 50))

# Phenotype generation
pheno_data = pheno.generate(model = "AD", geno_data = geno_data, effects = effects, sigma.error = sigma_error)

