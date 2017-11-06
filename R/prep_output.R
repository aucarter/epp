new.inf <- function(mod, fp) {
  attr(mod, "rvec")[fp$proj.steps %% 1 == 0.5] * (rowSums(mod[,-1,1]) + fp$relinfectART * rowSums(mod[,-1,-1])) / rowSums(mod) * mod[,1,1]
}

suscept.pop <- function(mod) {
  mod[,1,1]
}

plwh <- function(mod) {
  rowSums(mod[,-1,])
}

art <- function(mod) {
  rowSums(mod[,-1,-1])
}

total.pop <- function(mod) {
  rowSums(mod)
}

nat.draws <- function(result) {
  nat.inf <- Reduce('+', lapply(result, function(x){x$new.inf}))
  nat.suscept <- Reduce('+', lapply(result, function(x){x$suscept.pop}))

  nat.incid <- nat.inf/nat.suscept

  nat.plwh <- Reduce('+', lapply(result, function(x){x$plwh}))
  nat.total.pop <- Reduce('+', lapply(result, function(x){x$total.pop}))

  nat.art <- Reduce('+', lapply(result, function(x){x$art}))

  nat.art.cov <- nat.art / nat.plwh
  nat.prev <- nat.plwh/nat.total.pop
  output <- list(prev=nat.prev, incid=nat.incid, art = nat.art.cov, art_num=nat.art, pop=nat.total.pop)
  return(output)
}
# Simulate fit and output incidence and prevalence
simfit.gbd <- function(fit, random_walk = F){
  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))
  # Random-walk projection method
  if(random_walk){
    if(exists("eppmod", where=fit$fp) && fit$fp$eppmod == "rtrend")
      stop("Random-walk projection is only used with r-spline model")

    fit$rvec.spline <- sapply(fit$param, "[[", "rvec")
    firstidx <- (fit$likdat$firstdata.idx-1)/fit$fp$dt+1 # which(fit$fp$proj.steps == fit$fp$tsEpidemicStart) ## assume SD only depends on rvec of years with ANC/survey data
    lastidx <- (fit$likdat$lastdata.idx-1)/fit$fp$dt+1

    firstidx.mean <- which(fit$fp$proj.steps==min(epp.art[(epp.art$m.val+epp.art$f.val)>0,]$year)) # which(fit$fp$proj.steps == fit$fp$tsEpidemicStart) 

    ### Get the year of prev data that starts to decrease in the recent years
    hhs.dt <- fit$likdat$hhslik.dat
    if(nrow(fit$likdat$hhslik.dat)!=0){
      hhs.dt <- hhs.dt[order(hhs.dt$year),]
      hhs.dt$dif <- c(0, diff( hhs.dt$prev))
      proj.year <- hhs.dt[hhs.dt$dif<=0,]$year[nrow(hhs.dt[hhs.dt$dif<=0,])-1]
      proj.year.idx <- ifelse(length(proj.year)>0, which(fit$fp$proj.steps==proj.year), 0)
    } else {proj.year.idx <- 0}
    ## replace rvec with random-walk simulated rvec
    fit$param <- lapply(fit$param, function(par){
      ### compare the prev year with the lowest rvec year  
      rvec.min.idx <- which(par$rvec==min(par$rvec))
      first_projidx <- ifelse(proj.year.idx>rvec.min.idx, proj.year.idx, rvec.min.idx)

      par$rvec <- sim_rvec_rwproj(par$rvec, firstidx, lastidx, first_projidx, firstidx.mean, fit$fp$dt); par}) 
  }
  
  fp.list <- lapply(fit$param, function(par) update(fit$fp, list=par))
  mod.list <- lapply(fp.list, simmod)
  
  fit$rvec <- sapply(mod.list, attr, "rvec")
  fit$prev <- sapply(mod.list, prev)
  fit$incid <- mapply(incid, mod = mod.list, fp = fp.list)
  fit$popsize <- sapply(mod.list, rowSums)
  # fit$pregprev <- mapply(epp::fnPregPrev, mod.list, fp.list)

  ## GBD additions
  fit$new.inf <- mapply(new.inf, mod = mod.list, fp = fp.list)
  fit$suscept.pop <- sapply(mod.list, suscept.pop)
  fit$plwh <- sapply(mod.list, plwh)
  fit$total.pop <- sapply(mod.list, total.pop)
  fit$art <- sapply(mod.list, art)
  return(fit)
}
