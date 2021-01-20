

####################################################################################################################################
#
# Project: time-varying incremental propensity scores (IPS)
#
# Purpose: Apply time-varying IPS to the EAGeR data (exposure = aspirin, outcome= hCG pregnancy)
#
# Authors: Kwangho Kim, Ashley Naimi, Jacqueline Rudolph
#
# Last Update: 13 Jan 2020
#
###################################################################################################################################


packages <- c("tidyverse", "ranger", "reshape2", "latex2exp")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}


# -----------------------------------------------------------------------------------
# Import data
# -----------------------------------------------------------------------------------

source("./inc_dropout_utils_hcg.R")
eager <- read.table("../data/eager_tvar.txt", header = TRUE) 

# -----------------------------------------------------------------------------------
# setup: manipulate data structure so that it becomes balanced across id
# output = aspirin dataset (elongated)
# -----------------------------------------------------------------------------------

#Make it so everyone has maximum # of weeks of follow-up (here, 26 weeks)
  #If someone had event or was censored before 26 weeks, duplicate their last record
  #This reorders the data as well, by time then ID
#eager_elg <- proprocessing(data=eager, save.to.rds=TRUE) 
eager_elg <- readRDS(file = "../data/eager_elg.rds")

# Create Dropout (Rt) variable (indicating subject is still in the study):
    # JR: Rt is the opposite of the variable "drop" I created when setting up the data
    # Rename week variable to time
eager_elg <- eager_elg %>% 
  mutate(Rt = 1 - drop) %>% 
  rename(c("time" = "week",
           "delta" = "conception")) 
  

# -----------------------------------------------------------------------------------
#  Dataframe preparation & parameter setup 
# -----------------------------------------------------------------------------------

# lower bound for omega: 
    #Essentially being used to truncate the weights to avoid extreme values
delta.omega = 0.01 

# subject info & outcome of interest (delta = hCG pregnancy)
dat.names = c("id", "time", "treatment", "compliance", "delta", "last", "Rt")

# baseline covariates
cov.bs.names = c("id", "age", "BMI", "smoke", "eligibility")

# time-dependent covariates
cov.td.names = c("id", "time", "nausea", "bleed")

# prepare dataframe
dat <- select(eager_elg, all_of(dat.names)); dim(dat)
  x.cov.bs <- select(filter(eager_elg, time==1), all_of(cov.bs.names)); dim(x.cov.bs)
  x.cov.td <- select(eager_elg, all_of(cov.td.names)); dim(x.cov.td)

# sort in time & id
dat <- dat[order(dat$time, dat$id), ]
  x.cov.td <- x.cov.td[order(x.cov.td$time, x.cov.td$id),]
  x.cov.bs <- x.cov.bs[order(x.cov.bs$id),]

# ---  create treatment variable ----

#dat$A <- dat$compliance #Main exposure definition
dat$A <- dat$compliance * dat$treatment

# set of unique subject id's
uniq_ids <- unique(eager_elg$id)

# number of subjects
n <- length(uniq_ids)

# maximum timepoint (T)
max.time <- max(dat[dat$last==1, "time"]) # T

# reset maximum timepoint to be evaluated if necessary
ntimes <- 26 # <= max.time

# timepoint vector
times.vec <- ntimes:1

# outcome of interest (hCG pregnancy)
outcome.name <- "delta"

# for now, use Random Forest for computational efficiency
unisance.est.md <- "RF" 

# final dataframe to be analayzed after fixing our outcome and maximum timepoint
datf.names <- c("id", "time", outcome.name, "A", "Rt")
dat.f <-  dat[ , datf.names]
colnames(dat.f)[3] <- "Y"
dat.f <- dat.f[dat.f$time <= ntimes, ]; dim(dat.f)
  x.cov.td <- x.cov.td[x.cov.td$time <= ntimes, ]; dim(x.cov.td)

# delta vector: this is the vector of increments to use as interventions
delta.seq <- c(seq(0.30, 1-0.05, by=0.05), seq(1, 3, by=0.1))

# sample splitting
nsplits = 2 # > 1


# -----------------------------------------------------------------------------------
#  Implementation of Algorithm 1
# -----------------------------------------------------------------------------------

# mean value estimation
set.seed(123) # There is sampling, but I didn't notice a seed set anywhere?
estimation.result <- 
  estimation.sample.splitting(dat.f, x.cov.bs, x.cov.td, n, ntimes, delta.seq, nsplits)
kvals <- estimation.result$kvals #mean estimate on group k units
ifvals <- estimation.result$ifvals #influence function ests
est.eff <- colMeans(kvals)

# compute asymptotic variance
interval.list <- variance.bootstrap(ifvals)
eff.ll <- interval.list$eff.ll
eff.ul <- interval.list$eff.ul
eff.ll2 <- interval.list$eff.ll2
eff.ul2 <- interval.list$eff.ul2

# results
res <- data.frame(delta = delta.seq, 
                  est = est.eff, 
                  lower = eff.ll, 
                  upper = eff.ul)
#write_csv(res, file="../results/ips_tvar.csv")

# Draw figure
jpeg(file="../figures/ips_tvar.jpeg", height=5, width=6, unit="in", res=300)
plot.delta.curve(delta.seq, est.eff, eff.ll, eff.ul, eff.ll2, eff.ul2, outcome.name)
dev.off()
