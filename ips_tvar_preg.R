

###########################################################################
#
# Project: Estimation of time-varying incremental effect in EAGeR
#
# Authors: Kwangho Kim, Jacqueline Rudolph
#
# Last Update: 11 Feb 2022
#
###########################################################################

# Necessary packages
packages <- c("tidyverse", "reshape2", "latex2exp", "SuperLearner", "KernelKnn", "ranger")
for (package in packages) {
  library(package, character.only=T)
}

# Functions to run analysis
source("./inc_dropout_utils_preg.R")

# What exposure?
exposure <- "aspirin" # or "placebo"


# Read in data ------------------------------------------------------------

eager <- read.table("../data/eager_tvar.txt", header = TRUE) 


# Elongate data -----------------------------------------------------------

#Make it so everyone has maximum # of weeks of follow-up (here, 26 weeks)
  #If someone had event or was censored before 26 weeks, duplicate their last record
  #This reorders the data as well, by time then ID
eager_elg <- proprocessing(data=eager, save.to.rds=TRUE) 

# If data has already been set up
# eager_elg <- readRDS(file = "../data/eager_elg.rds")

# Create Dropout (Rt) variable (indicating subject is still in the study):
eager_elg <- eager_elg %>% 
  mutate(Rt = 1 - drop) %>% 
  rename("time" = "week") 
  

# Algorithm set-up --------------------------------------------------------

# Lower bound for omega: 
    # Essentially being used to truncate the weights to avoid extreme values
delta.omega = 0.01 

# Subject info & outcome of interest
dat.names = c("id", "time", "treatment", "compliance", "conception", "last", "Rt")

# Baseline covariates
cov.bs.names = c("id", "age", "BMI", "smoke", "eligibility")

# Time-dependent covariates
cov.td.names = c("id", "time", "nausea", "bleed")

# Prepare dataframe
dat <- select(eager_elg, all_of(dat.names)); dim(dat)
  x.cov.bs <- select(filter(eager_elg, time==1), all_of(cov.bs.names)); dim(x.cov.bs)
  x.cov.td <- select(eager_elg, all_of(cov.td.names)); dim(x.cov.td)

# Sort in time & id
dat <- dat[order(dat$time, dat$id), ]
  x.cov.td <- x.cov.td[order(x.cov.td$time, x.cov.td$id),]
  x.cov.bs <- x.cov.bs[order(x.cov.bs$id),]

# Treatment variable
if (exposure=="aspirin") {
  dat$A <- dat$compliance * dat$treatment # Comply with aspirin
} else if (exposure=="placebo") {
  dat$A <- dat$compliance * (1 - dat$treatment) # Comply with placebo
}

# Set of unique subject id's
uniq_ids <- unique(eager_elg$id)

# Number of subjects
n <- length(uniq_ids)

# Maximum timepoint (T)
max.time <- max(dat[dat$last==1, "time"])

# Reset maximum timepoint to be evaluated if necessary
ntimes <- 26 # <= max.time

# Timepoint vector (T >>> 1)
times.vec <- ntimes:1

# Outcome of interest (incidence of hCG confirmed pregnancy)
outcome.name <- "conception"

# Use Random Forest or SuperLearner ensemble for nuisance component modeling
unisance.est.md <- "SuperLearner" # "SuperLearner", "RF"

# Final dataframe to be analayzed after fixing our outcome and maximum timepoint
datf.names <- c("id", "time", outcome.name, "A", "Rt")
dat.f <-  dat[ , datf.names]
colnames(dat.f)[3] <- "Y"
dat.f <- dat.f[dat.f$time <= ntimes, ]; dim(dat.f)
  x.cov.td <- x.cov.td[x.cov.td$time <= ntimes, ]; dim(x.cov.td)

# Delta vector: this is the vector of increments to use as interventions
delta.seq <- c(seq(0.30, 1-0.05, by=0.05), seq(1, 3, by=0.1))

# Sample splitting
nsplits <- 2 # > 1


# Implement algorithm -----------------------------------------------------

# Mean value estimation
set.seed(123)
estimation.result <- 
  estimation.sample.splitting(dat.f, x.cov.bs, x.cov.td, uniq_ids, n, ntimes, delta.seq, nsplits)
kvals <- estimation.result$kvals # Mean estimate on group k units
ifvals <- estimation.result$ifvals # Influence function ests
est.eff <- colMeans(kvals)

# Compute asymptotic variance
interval.list <- variance.bootstrap(ifvals)
eff.ll <- interval.list$eff.ll
eff.ul <- interval.list$eff.ul
eff.ll2 <- interval.list$eff.ll2
eff.ul2 <- interval.list$eff.ul2

# 95% CI for the risk difference between two deltas
pt.ci.diff(delta.tgt=0.3, delta.ref = 1.0)


# Summarize results -------------------------------------------------------

# Output results
res <- data.frame(delta = delta.seq, 
                  est = est.eff, 
                  lower = eff.ll, 
                  upper = eff.ul)
write_csv(res, file=paste0("../results/ips_tvar_", exposure, ".csv"))
