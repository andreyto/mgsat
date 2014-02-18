### species matrix with 20 species abundances (mean = 50, sd = 10)
### one time variable, with 3 timepoints, which should be tested
### and a factor denoting sites that were repeatedly sampled (site)

## Load packages
require(vegan)

### Data:
sp <- matrix(rnorm(3 * 6 * 20, 50, 10), nrow = 3 * 6, ncol = 20,
             dimnames = list(1:18, paste("Sp", 1:20, sep = "")))

time <- as.ordered(rep(1:3, 6))
site <- gl(6, 3)
cbind(site, time, sp)

### add time effect at timepoint 3,
### this will effect will be tested by adonis():
sp_1 <- sp
sp_1[time==3,] <-  sp[time==3,] + rnorm(20, 10, 1)
cbind(site, time, sp_1)

### choose which species set to test:
test_sp <- sp_1

### computing the true R2-value

### (btw, using dist() defaults to euclidean distance):
print(fit <- adonis(dist(test_sp) ~ time, permutations=1))

### number of perms
B <- 1999

### setting up frame which will be populated by
### random r2 values:
pop <- rep(NA, B + 1)

### the first entry will be the true r2:
pop[1] <- fit$aov.tab[1, 5]

### set up a "permControl" object:
### we turn off mirroring as time should only flow in one direction
ctrl <- permControl(strata = site, within = Within(type = "series", mirror = FALSE))

### Number of observations:
nobs <- nrow(test_sp)

### check permutation (...rows represent the sample id):
### ..they are ok!
### within in each repeated sample (= sites) timepoints are shuffled,
### with keeping the sequence intact (e.g., for  site 1: 1,2,3 - 2,3,1 - 3,2,1)
shuffle(nobs, control = ctrl)

### loop:
### in adonis(...) you need to put permutations = 1, otherwise 
### adonis will not run
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand <- adonis(dist(test_sp) ~ time[idx], permutations = 1)
  pop[i] <- fit.rand$aov.tab[1, 5]
}

### get the p-value:
print(pval <- sum(pop >= pop[1]) / (B + 1))
### [1] 0.0035

### the sign. p-value supports the H1 (->there is a time effect).
### ..and the fact that samples are not iid is allowed by
### the customized perms - so this p-value is trustworthy as opposed
### to tests not acknowledging dependency of data points..

### test sp set without an effect:
### replace test_sp with sp set without effect:
test_sp <- sp

### now re-run the script and see the result:
### it is insign. - as expected:

### setting up frame which will be populated by
### random r2 values:
pop <- rep(NA, B + 1)

### computing the true R2-value:
print(fit <- adonis(dist(test_sp) ~ time, permutations = 1))

### the first entry will be the true r2:
pop[1] <- fit$aov.tab[1, 5]

### run the loop:
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand <- adonis(dist(test_sp) ~ time[idx], permutations = 1)
  pop[i] <- fit.rand$aov.tab[1, 5]
}
print(pval <- sum(pop >= pop[1]) / (B + 1))
### [1] 0.701

## make a histogram to see random R2-values and the true one:
hist(pop, xlab = "Population R2")
abline(v = pop[1], col = 2, lty = 3)
text(0.08, 300, paste("true R2,\np = ", pval, sep = ""))