####################################################################################
#   A set of utility functions for
#    "Incremental Intervention Effects in Studies 
#     with Many Timepoints, Repeated Outcomes, and Dropout." 
####################################################################################

# -----------------------------------------------------------------------------------
# setup and preprocessing
# -----------------------------------------------------------------------------------

proprocessing <- function(aspirin, save.to.rds=TRUE) {
  # preprocess the given long-form aspirin dataset
  # so that it can be balanced in the number of timepoints
  # across different units
  
  # ARGUMENTS:
  #   aspirin: [dataframe] long-form aspirin dataset
  
  # RETURN(S):
  #   aspirin_elg: [dataframe] new aspirin dataset that has been
  #                balanced in the number of timepoints across 
  #                different units
  
  
  aspirin <- create.last.variable(aspirin)
  
  # create separete columns for pregnancy outcome
  aspirin$live_birth <- ifelse(aspirin$status_new == "live birth" & aspirin$last == 1,1,0)
  aspirin$fetal_loss <- ifelse(aspirin$status_new == "pregnancy loss" & aspirin$last == 1,1,0)
  aspirin$efuwp <- ifelse(aspirin$status_new == "efuwp" & aspirin$last == 1,1,0)
  aspirin$withdraw <- ifelse(aspirin$status_new == "withdrawal" & aspirin$last == 1,1,0)
  
  max.time <- max(aspirin[aspirin$last==1,"j"]) # T
  
  aspirin_elg <- elongate.final.val(aspirin, max.time, save.to.rds=save.to.rds)
  
  return(aspirin_elg)
}


create.last.variable <- function(aspirin) {
  # Create new 'last' indicator variable to denote that either 
  # last_id==1 (outcome is realized) or the study ends (t=T)
  # We do this to include subjects whose outcome is not observed but status still exists
  
  # ARGUMENTS:
  #   aspirin: aspirin dataset in the long form
  
  # RETURN(S):
  #   aspirin dataset with the added 'last' indicator
  
  aspirin_ <- aspirin
  uniq_ids <- unique(aspirin$id)
  N <- dim(aspirin)[1]
  
  aspirin_$last <- 0
  for (i in 1:(N-1)) {
    if (aspirin_[(i+1),"first_id"] == 1) {
      aspirin_[i, "last"] <- 1
    }
    if (i == N-1) {
      aspirin_[i+1, "last"] <- 1
    }
  }
  
  if (sum(aspirin_$last) != length(uniq_ids)) {
    stop("sum(aspirin$last) != length(uniq_ids)")
  } # [1] 1224
  
  return (aspirin_)
}


elongate.final.val <- function(aspirin, max.time, save.to.rds=TRUE) {
  
  # elongate data for each subject up to the maximun timepoint (T) 
  # by padding its (dummy) final value after the time point at which no more 
  # data is collected on the subject (van der Laan & Robins, 2003 p.314)
  
  # ARGUMENTS:
  #   aspirin: [dataframe] aspirin dataset in the long form
  #            Must have the 'last' indicator
  #   max.time: [int] T
  
  # RETURNS:
  #   aspirin_elg: [dataframe] new aspirin dataset that has been
  #                balanced in the number of timepoints across 
  #                different units
  
  uniq_ids <- unique(aspirin$id)
  n <- length(uniq_ids)
  aspirin_elg <- aspirin
  
  print ("Data elongation Started")
  for (i in uniq_ids) {
    if (i %% 100 == 0) {
      print (paste("task ", toString(i)," complete", sep=""))
    }
    # final value for subject i
    last_dat <- aspirin[aspirin$id == i & aspirin$last == 1, ]
    if (last_dat$j == max.time) { 
      # skip if j = max_T
      next
    } else {
      # add dummy values
      last_dat$last <- 0
      final.val.dummy.df <- last_dat[rep(1, each=(max.time - last_dat$j)), ]
      final.val.dummy.df$j <- (last_dat$j+1):max.time
      aspirin_elg <- rbind(aspirin_elg, final.val.dummy.df)
    }
  }
  print ("Finished")
  
  # inconsistency check
  if (dim(aspirin_elg)[1] != (n * max.time)) {
    stop ("data size inconsistent")
  }
  
  # sorting by (time, id)
  aspirin_elg <- aspirin_elg[order(aspirin_elg$j, aspirin_elg$id),] 
  if (save.to.rds) {
    saveRDS(aspirin_elg, file = "aspirin_new_elg.rds")
    print (paste("saved to: ", "aspirin_new_elg.rds"))
  }
  
  return(aspirin_elg)
}


# -----------------------------------------------------------------------------------
# Implementation of Algorithm 1
# -----------------------------------------------------------------------------------

####################################################################################

# Compute the proposed estimator based on Algorithm 1
estimation.sample.splitting <- 
  function(dat.f, x.cov.bs, x.cov.td, uniq_ids, n, ntimes, delta.seq, nsplits){
    
    # ARGUMENTS:
    #   dat.f [dataframe] dataframe that contains subject info & outcome of interest
    #   x.cov.bs [dataframe] dataframe that contains baseline covariates (id, cov.name.1, ...)
    #   x.cov.td [dataframe] dataframe in the long form 
    #            that contains time-varying covariates (id, time, cov.name.1, ...)
    #   n [int] numher of the subject
    #   ntimes [int] the maximum number of timepoints (T)
    #   delta.seq [vector] a sequence of delta values
    #   nsplits [int] the number of splits (K)
    
    # RETURNS:
    #   ifvals: influence function estimates
    #   kvals: mean estimate on group k units
    #   caches
    
    # setup storage
    n.delta <- length(delta.seq)
    ifvals <- matrix(nrow=n, ncol=n.delta)
    kvals <- matrix(nrow=nsplits,ncol=n.delta)
    est.eff <- rep(NA,n.delta)
    
    W.long.delta = list()
    cumW.long = list()
    V.long.delta = list()
    M.long.delta = list()
    cumW.long.delta = list()
    
    s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])
    # slong <- rep(s, ntimes) # for (time, id) order
    
    Rt.prev <- c(rep(1,n),dat.f$Rt[1:(length(dat.f$Rt)-n)])
    
    ## Big loop started
    for (split in 1:nsplits){ 
      print(paste("split",split)); flush.console()
      id.split <- x.cov.bs[s == split,"id"]
      # all(id.split, uniq_ids[s == split])
      
      ## Step 1 & 2: estimate pi and omega functions
      print("  fitting pi and omega functions.."); flush.console()
      trt.mat <- data.frame(id = matrix(uniq_ids, nrow=n, ncol=1))
      omega.mat <- data.frame(id = matrix(uniq_ids, nrow=n, ncol=1))
      
      for (t in 1:ntimes) {
        id.t <- obsvbleIdSet(t-1, dat.f)
        x.trt.t <- buildCovMat(t, id.t, dat.f, x.cov.bs, x.cov.td, maxT = ntimes)
        R.next.t <- dat.f[dat.f$time==t & dat.f$id %in% id.t, c("id", "Rt")]
        x.trt.R.t <- merge(x.trt.t, R.next.t, by="id")
        
        trt.t.formula <- as.formula(paste(paste("A_",t,sep=""),"~.",sep = ""))
        omega.t.formula <- as.formula(paste("Rt","~.",sep = ""))
        if (nsplits == 1) {
          warning("nsplits should be greater than 1: naive Z-estimator used", call. = FALSE)
          trt.mod <- ranger(trt.t.formula, dat=x.trt.t[,-1])
          omega.mod <- ranger(omega.t.formula, dat=x.trt.R.t[,-1])
        } else {
          if (unisance.est.md == "RF"){
            trt.mod <- ranger(trt.t.formula, dat=x.trt.t[!(x.trt.t$id %in% id.split), -1] , write.forest = TRUE)
            omega.mod <- ranger(omega.t.formula, dat=x.trt.R.t[!(x.trt.R.t$id %in% id.split), -1] , write.forest = TRUE)
          } else if (unisance.est.md == "SuperLearner") {
            trt.mod <- SuperLearner(Y = x.trt.t[!(x.trt.t$id %in% id.split), c(paste("A_",t,sep=""))],
                                    X = x.trt.t[!(x.trt.t$id %in% id.split), !names(x.trt.t) %in% c(paste("A_",t,sep=""), "id")], 
                                    family = gaussian(), verbose = FALSE, cvControl = list(V=2), 
                                    SL.library = c("SL.ksvm", "SL.ranger", "SL.kernelKnn"))
            omega.mod <- SuperLearner(Y = x.trt.R.t[!(x.trt.R.t$id %in% id.split), c("Rt")],
                                     X = x.trt.R.t[!(x.trt.R.t$id %in% id.split), !names(x.trt.R.t) %in% c("Rt", "id")], 
                                     family = gaussian(), verbose = FALSE, cvControl = list(V=2), 
                                     SL.library = c("SL.ranger", "SL.kernelKnn"))
            
          } else {
            stop("other methods to be implemented")
          }
        }
        
        ps.t <- if (unisance.est.md == "RF"){ 
          predict(trt.mod, data=x.trt.t[,!names(x.trt.t) %in% c(paste("A_",t,sep=""), "id")])$predictions
        } else if (unisance.est.md == "SuperLearner") {
          as.vector(predict(trt.mod, x.trt.t[,!names(x.trt.t) %in% c(paste("A_",t,sep=""), "id")], onlySL = TRUE)$pred)
        } else {
          stop("other methods to be implemented")
        }
        
        omega.t <- if (unisance.est.md == "RF"){ 
          predict(omega.mod, data=x.trt.R.t[,!names(x.trt.R.t) %in% c("Rt", "id")])$predictions
        } else if (unisance.est.md == "SuperLearner") {
          as.vector(predict(omega.mod, x.trt.R.t[,!names(x.trt.R.t) %in% c("Rt", "id")], onlySL = TRUE)$pred)
        } else {
          stop("other methods to be implemented")
        }
        
        # set the lower bound for omega.t
        omega.t[which(omega.t<delta.omega)] <- delta.omega
        
        id.ps.t <- cbind(id = id.t, ps.t); colnames(id.ps.t)[2] <- t; #colnames(id.ps.t)[2] <- paste("ps_",t,sep="")
        id.omega.t <- cbind(id = id.t, omega.t); colnames(id.omega.t)[2] <- t; #colnames(id.omega.t)[2] <- paste("omega_",t,sep="")
        trt.mat = merge(trt.mat, id.ps.t, by="id", all.x = TRUE)
        omega.mat = merge(omega.mat, id.omega.t, by="id", all.x = TRUE)
      }
      
      trt.long <- melt(trt.mat, id.vars = "id")
      omega.long <- melt(omega.mat, id.vars = "id")
      
      ## Step 3 - 5
      for (j in 1:n.delta){
        print(paste("  delta",j)); flush.console()
        delta <- delta.seq[j]
        
        # Step 3:
        W.long <- omega.long
        W.long$value <- 
          (delta*dat.f$A + 1-dat.f$A)/(delta*trt.long$value + 1-trt.long$value) * 1/omega.long$value
        # (delta*dat.f$A + 1-dat.f$A)/(delta*trt.long$value + 1-trt.long$value) * Rt.prev/omega.long$value
        
        cumW.mat <- as.data.frame(aggregate(W.long$value, by=list(W.long$id), cumprod)[[2]])
        cumW.mat$id <- uniq_ids; colnames(cumW.mat)[1:ntimes] <- 1:ntimes;
        cumW.long <- melt(cumW.mat, id.vars = "id")
        
        # Step 4: fit outcome models
        M.mat <- data.frame(id = matrix(uniq_ids, nrow=n, ncol=1))
        M.mod <- vector("list", ntimes) 
        id.tp1 <- obsvbleIdSet(ntimes, dat.f)
        id.M.tp1 <- dat.f[dat.f$time==ntimes & dat.f$id %in% id.tp1, c("id", "Y")]
        print("    fitting psuedo regressions.."); flush.console()
        
        for (t in times.vec){
          colnames(id.M.tp1)[2] <- "Mt"
          x.trt.t <- buildCovMat(t, id.tp1, dat.f, x.cov.bs, x.cov.td, maxT = ntimes)
          x.M.t <- merge(x.trt.t, id.M.tp1, by="id")
          if (nsplits == 1) {
            M.mod[[t]] <- ranger(Mt ~ ., dat=x.M.t[,-1], write.forest = TRUE)
          } else {
            if (unisance.est.md == "RF"){
              M.mod[[t]] <- ranger(Mt ~ ., dat=x.M.t[!(x.M.t$id %in% id.split), -1], write.forest = TRUE)
            } else if (unisance.est.md == "SuperLearner") {
              M.mod[[t]] <- SuperLearner(Y = x.M.t[!(x.M.t$id %in% id.split), c("Mt")],
                                          X = x.M.t[!(x.M.t$id %in% id.split), !names(x.M.t) %in% c("Mt", "id")], 
                                          family = gaussian(), verbose = FALSE, cvControl = list(V=2), 
                                          SL.library = c("SL.ksvm", "SL.ranger", "SL.kernelKnn"))
            } else {
              stop("other methods to be implemented")
            }
          }
          
          id.t <- obsvbleIdSet(t-1, dat.f)
          x.trt.t.pred <- buildCovMat(t, id.t, dat.f, x.cov.bs, x.cov.td, maxT = ntimes)
          a.idx = paste(c("A", t), collapse = "_")
          newx0 <- newx1 <- x.trt.t.pred[,-1]
          newx1[,a.idx] <- 1 # (H_t, 1)
          m1 <- if (unisance.est.md == "RF") {
            predict(M.mod[[t]], data=newx1)$predictions
          } else if (unisance.est.md == "SuperLearner") {
            as.vector(predict(M.mod[[t]], newx1, onlySL = TRUE)$pred)
          } else {
            stop("other methods to be implemented")
          }
          newx0[,a.idx] <- 0  # (H_t, 0)
          m0 <- if (unisance.est.md == "RF") {
            predict(M.mod[[t]], data=newx0)$predictions
          } else if (unisance.est.md == "SuperLearner") {
            as.vector(predict(M.mod[[t]], newx0, onlySL = TRUE)$pred)
          } else {
            stop("other methods to be implemented")
          }
          
          pi.t <- trt.mat[trt.mat$id %in% id.t,t+1]
          # recursive regression formula in E.1
          M.tp1 <- (delta*m1*pi.t + m0*(1-pi.t))/(delta*pi.t + 1-pi.t) 
          id.M.tp1 <- cbind(id = id.t, M.tp1); colnames(id.M.tp1)[2] <- t; 
          
          # step 4-b
          omega.t <- omega.mat[omega.mat$id %in% id.t,t+1]
          A.t <- dat.f[dat.f$id %in% id.t & dat.f$time == t, c("A")]
          R.tp1 <- dat.f[dat.f$id %in% id.t & dat.f$time==t, c("Rt")]
          M.t <- (m1 - m0)*delta*(A.t - pi.t)*omega.t/(delta*pi.t + 1-pi.t) +
            delta*m1*(pi.t*omega.t - A.t*R.tp1) + m0*((1-pi.t)*omega.t - (1-A.t)*R.tp1)
          id.M.t <- cbind(id = id.t, M.t); colnames(id.M.t)[2] <- t; 
          
          M.mat = merge(M.mat, id.M.t, by="id", all.x = TRUE)
          
          id.tp1 <- id.t
        }
        
        M.long <- melt(M.mat, id.vars = "id")
        
        # step 5:
        V.long <- W.long
        V.long$value <- 1/(delta*dat.f$A + 1-dat.f$A)
        
        # cache
        W.long.delta[[j]] <- W.long
        cumW.long.delta[[j]] <- cumW.long
        V.long.delta[[j]] <- V.long
        M.long.delta[[j]] <- M.long
        
        # Step 6:
        phi.val <- ((cumW.long$value*dat.f$Y)[dat.f$time==ntimes] +
                      aggregate(cumW.long$value*V.long$value*M.long$value, by=list(dat.f$id), sum)[[2]])[s==split]
        
        ifvals[s==split,j] <- phi.val
        print(mean(phi.val, na.rm = T))
        kvals[split,j] <- mean(phi.val, na.rm = T)
      } 
    }
    
    return(list(ifvals=ifvals, kvals=kvals, s=s,
                W.long.delta=W.long.delta, cumW.long.delta=cumW.long.delta,
                V.long.delta=V.long.delta, M.long.delta=M.long.delta))
  }

####################################################################################

####################################################################################

# Function to create (\bar{X}_t, \bar{A}_t) = (\bar{H}_t, A_t)
buildCovMat <- 
  function(t, id.t, dat, x.cov.bs, x.cov.td, maxT, td.cov.names = NULL, append = FALSE, no.last.trt = FALSE)
  {
    # DESCRYPTION:
    #    Returns cbind(\bar{H}_t, A_t) with "id" column
    #    which is the historical values of covariates and exposures
    #    of given id set up to time t, only with the specified id set (id.t)
    #    if no.last.trt = TRUE, then we will just get \bar{H}_t
    
    # ARGUMENTS:
    #    id.t =  index set of "id"'s; data for each id should be available at time t
    
    # RETURNS:
    #   HA.t.bar [dataframe] cbind(\bar{H}_t, A_t) with id column
    
    if (t<1) return(NULL)
    
    # time-independent covariates
    if (!is.null(x.cov.bs)) {
      x.bs.t <- x.cov.bs[x.cov.bs$id %in% id.t, ]
      x.bs.t <- x.bs.t[order(x.bs.t$id),]
      # x.bs.t <- x.bs.t[,-1];
    } else {
      x.bs.t <- data.frame(matrix(NA, nrow=length(id.t), ncol=0))
    }
    
    # time-dependent covariates
    if (dim(x.cov.td)[2] >= 3) {
      if (is.null(td.cov.names)) {
        td.cov.names <- colnames(x.cov.td)[3:length(colnames(x.cov.td))]
      }
    } else {
      warning('No time-dependent variable used')
      td.cov.names <- NULL
    }
    x.td.t.pool <- x.cov.td[x.cov.td$id %in% id.t,]
    x.td.t.pool <- x.td.t.pool[order(x.td.t.pool$time, x.td.t.pool$id),]
    
    # treatment vector
    A.t.pool <- dat[dat$id %in% id.t, c("time", "id", "A")]
    A.t.pool <- A.t.pool[order(A.t.pool$time, A.t.pool$id),]
    
    # define a local function that returns td cov values at time t_ (<t)
    td.time.t <- function(t_, x.td.t.pool, td.cov.names, keep.id = F) { 
      x.td.t = x.td.t.pool[x.td.t.pool$time == t_, c("id", td.cov.names), drop = F]
      td.names = NULL
      for (nm in td.cov.names) td.names = c(td.names, paste(c(nm, t_), collapse = "_"))
      colnames(x.td.t)[-1] <- td.names
      if (keep.id == T) {
        return(x.td.t)
      } else {
        return(x.td.t[,-1, drop = F])  
      }
    }
    
    # build (\bar{H}_t, A_t)
    x.td.t.bar <- td.time.t(1, x.td.t.pool, td.cov.names, T)
    A.t.bar <- td.time.t(1, A.t.pool, "A", T)
    if (t>1) {
      for (k in 2:t) {
        k.temp.x <- td.time.t(k, x.td.t.pool, td.cov.names, T)
        x.td.t.bar <- merge(x.td.t.bar, k.temp.x, by="id")
        k.temp.a <- td.time.t(k, A.t.pool, "A", T)
        A.t.bar <- merge(A.t.bar, k.temp.a, by="id")
      }
    }
    
    HA.td.t.bar <- merge(x.td.t.bar, A.t.bar, by="id")
    
    if (append == T & t < maxT) { # force to repeat the final value until maxT 
      for (k in (t+1):maxT) {
        k.temp <- td.time.t(t, x.td.t.pool, td.cov.names)
        x.td.t.bar <- cbind(x.td.t.bar, k.temp)
      }
    }
    
    # drop the last column if needed
    if (no.last.trt) HA.td.t.bar <- HA.td.t.bar[,-ncol(HA.td.t.bar), drop = F]
    
    HA.t.bar <- if(is.null(x.cov.bs)) HA.td.t.bar else merge(x.bs.t, HA.td.t.bar, by="id")
    return(HA.t.bar)
    
  }

####################################################################################
####################################################################################

# Function to identify 'id's still observable at time t (id_t: R_t == 1)
obsvbleIdSet <- 
  function(t, dat)
  {
    # DESCRYPTION:
    #    Returns the vector of id's that have not withdrawn at t
    #    and will be staying in the next timepoint 
    #    (Typically we assume R_1 (or R_0) = 1)
    if (t==1 || t==0) {
      # we can always assume R_1 = 1
      id.t <- dat[dat$time == 1, "id"]  
    } else{
      id.t <- dat[dat$time == t & dat$Rt == 1, "id"]  
    }
    
    return(id.t)
    
  }

####################################################################################

# compute 95% asymptotic variance via bootstrapping
variance.bootstrap <- function(ifvals) {
  # compute both pointwise confidence interval (CI) and confidence band (CB)
  # where alpha = 0.05
  
  # ARGUMENTS:
  #   ifvals: value of estimated influence function 
  
  # RETURN(S):
  # eff.ll: lower bound of 95% pointwise CI
  # eff.ul: upper bound of 95% pointwise CI
  # eff.ll2: lower bound of 95% CB
  # eff.ul2: upper bound of 95% CB
  
  # 95% pointwise asymptotic variance
  sigma <- sqrt(apply(ifvals,2,var,na.rm=TRUE))
  ifvals.0 <- na.omit(ifvals)
  n0 <- dim(ifvals.0)[1]
  eff.ll <- est.eff-1.96*sigma/sqrt(n0); eff.ul <- est.eff+1.96*sigma/sqrt(n0)
  
  # confidence band based on multiplier bootstrapping
  eff.mat <- matrix(rep(est.eff,n0), nrow=n0, byrow=T)
  sig.mat <- matrix(rep(sigma,n0), nrow=n0, byrow=T)
  ifvals2 <- (ifvals.0 - eff.mat)/sig.mat
  nbs <- 10000 
  mult <- matrix(2*rbinom(n0*nbs,1,.5) - 1, nrow=n0, ncol=nbs)
  
  maxvals <- sapply(1:nbs, function(col){
    max(abs(apply(mult[,col]*ifvals2,2,sum)/sqrt(n0))) 
  } 
  )
  calpha <- quantile(maxvals, 0.95)
  eff.ll2 <- est.eff - calpha*sigma/sqrt(n0); eff.ul2 <- est.eff + calpha*sigma/sqrt(n0)
  
  return(list(eff.ll=eff.ll, eff.ul=eff.ul, eff.ll2=eff.ll2, eff.ul2=eff.ul2))
}

####################################################################################

# draw curve for estimates at varying delta with 95% pointwise CI and confidence band
plot.delta.curve <- 
  function(delta.seq, est.eff, eff.ll, eff.ul, eff.ll2, eff.ul2, outcome.name) {
    
    outcome.name.t <- paste(strsplit(outcome.name, "_")[[1]], collapse=" ")
    outcome.name.t <- simpleCap(outcome.name.t)
    
    max.y <- min(1.5,max(eff.ul2)); min.y <- max(-1,min(eff.ll2)); # min.y <- max(min(eff.ll2),0)
    plot(delta.seq, eff.ll2, type="l", col="firebrick", lty=3, lwd=2, xlab = NA, ylab = NA, ylim=c(min.y, max.y)); par(new=T)
    plot(delta.seq, eff.ul2, type="l", col="firebrick", lty=3, lwd=2, xlab = NA, ylab = NA, ylim=c(min.y, max.y)); par(new=T)
    plot(delta.seq, eff.ll, type="n", lty=3, lwd=2, xlab = NA, ylab = NA, ylim=c(min.y, max.y)); par(new=T)
    plot(delta.seq, eff.ul, type="n", lty=3, lwd=2, xlab = NA, ylab = NA, ylim=c(min.y, max.y)); par(new=T)
    polygon(c(delta.seq, rev(delta.seq)), c(eff.ul2, rev(eff.ll2)),
            col = "gainsboro", border = NA); par(new=T)
    polygon(c(delta.seq, rev(delta.seq)), c(eff.ul, rev(eff.ll)),
            col = "gray70", border = NA); par(new=T)
    plot(delta.seq, est.eff, type="l", lwd=1.5, col="dimgrey", 
         xlab=NA, ylab=NA, ylim=c(min.y, max.y));  par(new=F)
    abline(v = 1, col="darkblue", lwd=1.5, lty=3)
    title(main=paste("Estimated Probability of ",outcome.name.t," (T=", ntimes, ") ", sep=""),
          xlab=expression(paste("odds ratio ",delta,sep="")), ylab=expression(Psi(delta)))
    
  }


simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
