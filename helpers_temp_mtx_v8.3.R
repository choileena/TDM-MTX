nwn <- function(x) suppressWarnings(as.numeric(x))

checkInputMatrix <- function(m, def, defDuration = NA) {
  m[,'dose'] <- nwn(m[,'dose'])
  m[!is.na(m[,'dose']) & m[,'dose'] < 0, 'dose'] <- NA
  if(!is.na(defDuration)) {
    if(!'duration' %in% names(m)) {
      m[,'duration'] <- defDuration
    } else {
      m[,'duration'] <- nwn(m[,'duration'])
      m[is.na(m[,'duration']), 'duration'] <- defDuration
      m[m[,'duration'] < 0, 'duration'] <- defDuration
    }
  }
  m[,'time'] <- nwn(m[,'dt']) 
  m <- m[complete.cases(m),]
  if(nrow(m) == 0) {
    m <- def
  }
  m
}

##########################
###### MY MODEL ##########
##########################
model_mtx <- function(usrbsa, usrSCR, usrgender) { #define model parameters 
  THETA1 <- 9.41216979510633  # CL pop
  THETA2 <- 30.666641775106   # V pop
  THETA3 <- 0.866697110621852 # Q2 pop
  THETA4 <- 5.61864647680956 # V2 pop
  THETA5 <- 0.13515317483686 # Q3 pop
  THETA6 <- 10.6049699903601 #V3 pop
  beta_scr_cl = -0.556395940439738 ; beta_bsa_cl = 0.609402364492555 ; beta_gender = 0.013242740163542
  
  sigprop <- 0.2508599^2
  sigadd = 0
  Omega = diag(c(0.242508391327, 0.289352674694, 0.38419083468, #cl v1 q2
                 0.300674583119, 0.466055297477, 0.343596438913)) #v2 q3 v3
  Sigma = diag(c(sigprop,sigadd))
  list(ncpt = 3, Omega = Omega, Sigma = Sigma,
       THETA1 = THETA1, THETA2 = THETA2, THETA3 = THETA3, THETA4 = THETA4, THETA5 = THETA5, THETA6 = THETA6,
       beta_scr_cl = beta_scr_cl, beta_bsa_cl = beta_bsa_cl, beta_gender = beta_gender,
       bsa = usrbsa, SCR_mgdl = usrSCR, pt_gender = usrgender,beta_scr_cl = beta_scr_cl, beta_bsa_cl = beta_bsa_cl, beta_gender = beta_gender )
}

mymodel <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    
    ## --- Interpolate covariates at current time ---
    bsa_t <- cov_bsa[!is.na(cov_bsa)][1]
    #bsa_t     <- approx(covtimes, cov_bsa,    xout = time, rule = 2)$y
    #SCR_mmol <- unique(cov_scr[!is.na(cov_scr)])
    SCR_fun <- approxfun(covtimes, cov_scr, method = "constant", rule = 2)
    SCR_mmol <- SCR_fun(time)    
    pt_gender <- ifelse(cov_gender[1] == "Male", 0 , 1)
    beta_bsa_cl = parms[["beta_bsa_cl"]]
    beta_scr_cl = parms[["beta_scr_cl"]]
    beta_gender = parms[["beta_gender"]]
    

    ## --- Structural parameters ---
    # You can include covariate effects here later (e.g., CL ∝ BSA^0.75)
    CL <- THETA1 * ((bsa_t/1.97)^beta_bsa_cl)*((SCR_mmol/68.08)^(beta_scr_cl))*exp(beta_gender*pt_gender) * exp(ETA1)
    V1 <- THETA2 * (bsa_t/1.97) * exp(ETA2)
    Q2 <- THETA3 * (bsa_t/1.97) * exp(ETA3)
    V2 <- THETA4 * (bsa_t/1.97) * exp(ETA4)
    Q3 <- THETA5 * (bsa_t/1.97) * exp(ETA5)
    V3 <- THETA6 * (bsa_t/1.97) * exp(ETA6)
    
    ## --- Microconstants ---
    k10 <- CL / V1
    k12 <- Q2 / V1
    k21 <- Q2 / V2
    k13 <- Q3 / V1
    k31 <- Q3 / V3
    
    ## --- Infusion rate ---
    rate <- 0
    for (i in seq_along(dose_times)) {
      if (time >= dose_times[i] && time <= (dose_times[i] + dose_durs[i])) {
        rate <- rate + dose_amts[i] / dose_durs[i]
      }
    }
    
    ## --- Differential equations ---
    dA1dt <- rate - (k10 + k12 + k13) * A1 + k21 * A2 + k31 * A3
    dA2dt <- k12 * A1 - k21 * A2
    dA3dt <- k13 * A1 - k31 * A3
    
    list(c(dA1dt, dA2dt, dA3dt)) 
  })
}

simulate_mymodel <- function(dat_subj, THETA1, THETA2, THETA3, THETA4, THETA5, THETA6, 
                                        ETA1 = 0, ETA2 = 0, ETA3 = 0, ETA4 = 0, ETA5 = 0, ETA6 = 0, sigma = 0.01,
                                        beta_bsa_cl, beta_scr_cl, beta_gender) {
  # ensure sorted
  dat_subj <- dat_subj %>% arrange(time)
  
  # extract doses
  dose_amts  <- dat_subj$amt[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  dose_times <- dat_subj$time[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  dose_durs  <- dat_subj$dur[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  if (length(dose_durs) == 0) dose_durs <- rep(0, length(dose_amts))
  
  # pack parameters
  parms <- list(
    THETA1 = THETA1, THETA2 = THETA2, THETA3 = THETA3,
    THETA4 = THETA4, THETA5 = THETA5, THETA6 = THETA6,
    ETA1 = ETA1, ETA2 = ETA2, ETA3 = ETA3,
    ETA4 = ETA4, ETA5 = ETA5, ETA6 = ETA6,
    covtimes   = dat_subj$time,
    cov_bsa    = dat_subj$bsa,
    cov_scr    = dat_subj$SCR_mmol,
    cov_gender = dat_subj$pt_gender,
    dose_times = dose_times,
    dose_amts  = dose_amts,
    dose_durs  = dose_durs,
    beta_bsa_cl = beta_bsa_cl,
    beta_scr_cl = beta_scr_cl,
    beta_gender = beta_gender
  )
  
  # initial state
  state <- c(A1 = 0, A2 = 0, A3 = 0)
  times <- sort(unique(c(seq(0, max(dat_subj$time), by = 0.1), dat_subj$time)))
  
  # integrate
  #figure out tolerance
  out <- ode(y = state, times = times, func = mymodel, parms = parms, hmax = Inf)
  out <- as.data.frame(out)
  out$subject_id <- dat_subj$subject_id[1]
  
  SCR_fun <- approxfun(parms$covtimes, parms$cov_scr, method = "constant", rule = 2)
  SCR_mmol <- SCR_fun(times)    

  out$SCR_mmol = SCR_mmol
  
  # compute concentrations
  out <- out %>%
    rowwise() %>%
    mutate(
      bsa = dat_subj$bsa[1], #approx(parms$covtimes, parms$cov_bsa, xout = time, rule = 2)$y,
      pt_gender =  ifelse(dat_subj$pt_gender[1] == "Male", 0 , 1),
      CL = THETA1 * ((bsa/1.97)^beta_bsa_cl)*((SCR_mmol/68.08)^(beta_scr_cl))*exp(beta_gender*pt_gender)*exp(ETA1),
      V1 = THETA2 * (bsa/1.97)* exp(ETA2),
      Q2 = THETA3 * (bsa/1.97)* exp(ETA3),
      V2 = THETA4 * (bsa/1.97)* exp(ETA4),
      Q3 = THETA5 * (bsa/1.97)* exp(ETA5),
      V3 = THETA6 * (bsa/1.97)* exp(ETA6),
      Conc = A1 / V1
    ) %>%
    ungroup()
  
  out <- subset(out, time %in% dat_subj$time)
  out$evid <- dat_subj$evid
  
  return(out)
}

mapbayes_mymodel = function(eta,dat_subj, y, yt, Omega=mp$Omega, Sigma =mp$Sigma,
                            THETA1 = mp$THETA1, THETA2 = mp$THETA2, THETA3 = mp$THETA3, THETA4 = mp$THETA4, THETA5 = mp$THETA5, THETA6 = mp$THETA6,
                            beta_scr_cl = mp$beta_scr_cl, beta_bsa_cl = mp$beta_bsa_cl, beta_gender= mp$beta_gender){
  #suppressMessages(attach(mp))
  sigprop <- as.numeric(Sigma[1,1])
  sigadd <- as.numeric(Sigma[2,2])
  eta =  eta %>%  as.list
  # names(eta) <- names(init)
  eta_m <- eta %>% unlist %>% matrix(nrow=1)
  
  looppat <- simulate_mymodel(dat_subj, THETA1, THETA2, THETA3, THETA4, THETA5, THETA6, 
                                         ETA1 = eta$ETAcl, ETA2 = eta$ETAv1, 
                                         ETA3=eta$ETAq2, ETA4=eta$ETAv2, 
                                         ETA5 = eta$ETAq3, ETA6 = eta$ETAv3, beta_scr_cl = beta_scr_cl, beta_bsa_cl = beta_bsa_cl, beta_gender= beta_gender)
  
  newobs = looppat$Conc[which(looppat$time%in%yt)]
  
  # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339294/
  sig2j <- (newobs*sqrt(sigprop) + sqrt(sigadd))^2    #<------------------Check against pubmed
  sqwres <- log(sig2j) + (1/sig2j)*(y-newobs)^2
  
  nOn <- diag(eta_m %*% solve(Omega) %*% t(eta_m))
  return(sum(sqwres) + nOn)
}

pkprof_est_mymodel <- function(time_and_covs, mp) { 
  ncpt <- mp$ncpt
  
  y <- time_and_covs$dv[which(!is.na(time_and_covs$dv))[1]]
  yt <- round(time_and_covs$time[which(!is.na(time_and_covs$dv))[1]],1)
  
  init = eta <- c(ETAcl=.05, ETAv1=.05, ETAq2=.05, ETAv2=.05, ETAq3=.05, ETAv3=.05)
#  print(time_and_covs$SCR_mmol)
  dat_subj = subset(time_and_covs, !is.na(dv) | !is.na(amt))
  print(dat_subj)
  fit <- newuoa(par=init,
                fn = mapbayes_mymodel,
                dat_subj = dat_subj,
                y=y,
                yt=yt,
                Omega=mp$Omega,
                Sigma =mp$Sigma,
                THETA1 = mp$THETA1, THETA2 = mp$THETA2, THETA3 = mp$THETA3, THETA4 = mp$THETA4, THETA5 = mp$THETA5, THETA6 = mp$THETA6,
                beta_scr_cl = mp$beta_scr_cl, beta_bsa_cl = mp$beta_bsa_cl, beta_gender= mp$beta_gender)
  
  return(fit$par)
}

##########################
######### MTXPK ##########
##########################
model_mtxpk <- function(usrbsa, usrSCR) { #define model parameters 
  THETA1 <- 11  # CL pop
  THETA2 <- 16.5   # V pop
  THETA3 <- 0.602 # Q2 pop
  THETA4 <- 4.55 # V2 pop
  THETA5 <- 0.111 # Q3 pop
  THETA6 <- 13.1 #V3 pop
  beta_scr_cl = -0.247 
  
  sigprop <- 0.44
  sigadd = 0
  Omega = diag(c(0.08, #cl 
                 0.12, 0.13, 0.10)) #v2 q3 v3
  Sigma = diag(c(sigprop,sigadd))
  list(ncpt = 3, Omega = Omega, Sigma = Sigma,
       THETA1 = THETA1, THETA2 = THETA2, THETA3 = THETA3, THETA4 = THETA4, THETA5 = THETA5, THETA6 = THETA6,
       beta_scr_cl = beta_scr_cl, bsa = usrbsa, SCR_mmol = usrSCR )
}

mtxpk <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    
    ## --- Interpolate covariates at current time ---
    bsa_t <- cov_bsa[!is.na(cov_bsa)][1]
    #bsa_t     <- approx(covtimes, cov_bsa,    xout = time, rule = 2)$y
    SCR_fun <- approxfun(covtimes, cov_scr, method = "constant", rule = 2)
    SCR_mmol <- SCR_fun(time)    
    #SCR_mmol  <- approx(covtimes, cov_scr,    xout = time, rule = 2)$y
    
    beta_scr_cl = parms[["beta_scr_cl"]]

    ## --- Structural parameters ---
    CL <- THETA1 * ((bsa_t/1.73))*((SCR_mmol/29)^(beta_scr_cl))* exp(ETA1)
    V1 <- THETA2 * (bsa_t/1.73) 
    Q2 <- THETA3 * (bsa_t/1.73)
    V2 <- THETA4 * (bsa_t/1.73) * exp(ETA4)
    Q3 <- THETA5 * (bsa_t/1.73) * exp(ETA5)
    V3 <- THETA6 * (bsa_t/1.73) * exp(ETA6)
    
    ## --- Microconstants ---
    k10 <- CL / V1
    k12 <- Q2 / V1
    k21 <- Q2 / V2
    k13 <- Q3 / V1
    k31 <- Q3 / V3
    
    ## --- Infusion rate ---
    rate <- 0
    for (i in seq_along(dose_times)) {
      if (time >= dose_times[i] && time <= (dose_times[i] + dose_durs[i])) {
        rate <- rate + dose_amts[i] / dose_durs[i]
      }
    }
    
    ## --- Differential equations ---
    dA1dt <- rate - (k10 + k12 + k13) * A1 + k21 * A2 + k31 * A3
    dA2dt <- k12 * A1 - k21 * A2
    dA3dt <- k13 * A1 - k31 * A3
    
    list(c(dA1dt, dA2dt, dA3dt))
  })
}

simulate_mtxpk <- function(dat_subj, THETA1, THETA2, THETA3, THETA4, THETA5, THETA6, 
                                        ETA1 = 0, ETA4 = 0, ETA5 = 0, ETA6 = 0, sigma = 0.01,
                                        beta_scr_cl) {
  # ensure sorted
  dat_subj <- dat_subj %>% arrange(time)
  
  # extract doses
  dose_amts  <- dat_subj$amt[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  dose_times <- dat_subj$time[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  dose_durs  <- dat_subj$dur[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  if (length(dose_durs) == 0) dose_durs <- rep(0, length(dose_amts))
  
  # pack parameters
  parms <- list(
    THETA1 = THETA1, THETA2 = THETA2, THETA3 = THETA3,
    THETA4 = THETA4, THETA5 = THETA5, THETA6 = THETA6,
    ETA1 = ETA1,
    ETA4 = ETA4, ETA5 = ETA5, ETA6 = ETA6,
    covtimes   = dat_subj$time,
    cov_bsa    = dat_subj$bsa,
    cov_scr    = dat_subj$SCR_mmol,
    dose_times = dose_times,
    dose_amts  = dose_amts,
    dose_durs  = dose_durs,
    beta_scr_cl = beta_scr_cl
  )
  
  # initial state
  state <- c(A1 = 0, A2 = 0, A3 = 0)
  times <- sort(unique(c(seq(0, max(dat_subj$time), by = 0.1), dat_subj$time)))
  
  # integrate
  #figure out tolerance
  out <- ode(y = state, times = times, func = mtxpk, parms = parms)
  out <- as.data.frame(out)
  out$subject_id <- dat_subj$subject_id[1]
  
  SCR_fun <- approxfun(parms$covtimes, parms$cov_scr, method = "constant", rule = 2)
  SCR_mmol <- SCR_fun(times)    
  
  out$SCR_mmol = SCR_mmol
  
  # compute concentrations
  out <- out %>%
    rowwise() %>%
    mutate(
      bsa = dat_subj$bsa[1], #approx(parms$covtimes, parms$cov_bsa, xout = time, rule = 2)$y,
      pt_gender =  ifelse(dat_subj$pt_gender[1] == "Male", 0 , 1),
      CL = THETA1 * ((bsa/1.7))*((SCR_mmol/29)^(beta_scr_cl))*exp(ETA1),
      V1 = THETA2 * (bsa/1.7),
      Q2 = THETA3 * (bsa/1.7),
      V2 = THETA4 * (bsa/1.7)* exp(ETA4),
      Q3 = THETA5 * (bsa/1.7)* exp(ETA5),
      V3 = THETA6 * (bsa/1.7)* exp(ETA6),
      Conc = A1 / V1
    ) %>%
    ungroup()
  
  out <- subset(out, time %in% dat_subj$time)
  out$evid <- dat_subj$evid
  
  return(out)
}

mapbayes_mtxpk = function(eta,dat_subj, y, yt, Omega=mp$Omega,  Sigma =mp$Sigma,
                          THETA1 = mp$THETA1, THETA2 = mp$THETA2, THETA3 = mp$THETA3, THETA4 = mp$THETA4, THETA5 = mp$THETA5, THETA6 = mp$THETA6,
                          beta_scr_cl = mp$beta_scr_cl){
  #suppressMessages(attach(mp))
  sigprop <- as.numeric(Sigma[1,1])
  sigadd <- as.numeric(Sigma[2,2])
  eta =  eta %>%  as.list
  # names(eta) <- names(init)
  eta_m <- eta %>% unlist %>% matrix(nrow=1)
  
  looppat <- simulate_mtxpk(dat_subj, THETA1, THETA2, THETA3, THETA4, THETA5, THETA6, 
                                         ETA1 = eta$ETAcl, ETA4=eta$ETAv2, 
                                         ETA5 = eta$ETAq3, ETA6 = eta$ETAv3, beta_scr_cl = beta_scr_cl)
  
  newobs = looppat$Conc[which(looppat$time%in%yt)]
  
  # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339294/
  sig2j <- (newobs*sqrt(sigprop) + sqrt(sigadd))^2    #<------------------Check against pubmed
  sqwres <- log(sig2j) + (1/sig2j)*(y-newobs)^2
  
  nOn <- diag(eta_m %*% solve(Omega) %*% t(eta_m))
  return(sum(sqwres) + nOn)
}

pkprof_est_mtxpk <- function(time_and_covs, mp = mp_mtxpk) { 
  ncpt <- mp$ncpt
  
  y <- time_and_covs$dv[which(!is.na(time_and_covs$dv))[1]]
  yt <- round(time_and_covs$time[which(!is.na(time_and_covs$dv))[1]],1)
  
  init = eta <- c(ETAcl=.05, ETAv2=.05, ETAq3=.05, ETAv3=.05)

  dat_subj = subset(time_and_covs, !is.na(dv) | !is.na(amt))
  
  fit <- newuoa(par=init,
                fn = mapbayes_mtxpk,
                dat_subj = dat_subj,
                y=y,
                yt=yt,
                Omega=mp$Omega,
                Sigma =mp$Sigma,
                THETA1 = mp$THETA1, THETA2 = mp$THETA2, THETA3 = mp$THETA3, THETA4 = mp$THETA4, THETA5 = mp$THETA5, THETA6 = mp$THETA6,
                beta_scr_cl = mp$beta_scr_cl)
  
  return(fit$par)
}


##########################
######### overall ##########
##########################
combine_schedule <- function(dat, patdat, pats) { #estimated concs #what is pats
    dat$id <- 'Population'
    patdat$id <- 'Individual'
    schedule1 <- rbind(pats, patdat, dat)
    isIP1 <- schedule1$id %in% c("Individual","Population")
    schedule1$size <- as.numeric(isIP1)
    schedule1$Profile <- factor(ifelse(isIP1, schedule1$id, 'Simulated'))
    schedule1
}

#calculates concs
getpatdat <- function(dat_subj, mp, pk_est, mp_mtxpk, pk_est_mtxpk) {
  out_pop <- simulate_mymodel(dat_subj, mp$THETA1, mp$THETA2, mp$THETA3, mp$THETA4, mp$THETA5, mp$THETA6, 
                                         ETA1 = 0, ETA2 = 0, ETA3=0, ETA4=0, ETA5 = 0, ETA6 = 0, 
                                         beta_scr_cl = mp$beta_scr_cl, beta_bsa_cl = mp$beta_bsa_cl, beta_gender= mp$beta_gender) %>% 
    dplyr::select(-c(A1,A2, A3)) 
  
  out_pop_lower <- simulate_mymodel(dat_subj, THETA1 = 8.72097, THETA2 = 28.0139, THETA3 = 0.732417, THETA4 = 5.08519, THETA5 = 0.109888, THETA6 = 6.88222, 
                                         ETA1 = 0, ETA2 = 0, ETA3=0, ETA4=0, ETA5 = 0, ETA6 = 0, 
                                         beta_scr_cl = 0.416081, beta_bsa_cl = 0.392739, beta_gender= -0.083572) %>% 
    dplyr::select(-c(A1,A2, A3)) 
  
  out_pop_upper <- simulate_mymodel(dat_subj, THETA1 = 10.1582, THETA2 = 33.5705, THETA3 = 1.0256, THETA4 = 6.20806, THETA5 = 0.166228, THETA6 = 16.3414, 
                                               ETA1 = 0, ETA2 = 0, ETA3=0, ETA4=0, ETA5 = 0, ETA6 = 0, 
                                               beta_scr_cl = 0.74403, beta_bsa_cl = 0.945592, beta_gender= 0.110057) %>% 
    dplyr::select(-c(A1,A2, A3)) 
  
  IPRED = simulate_mymodel(dat_subj, mp$THETA1, mp$THETA2, mp$THETA3, mp$THETA4, mp$THETA5, mp$THETA6, 
                                      ETA1 = pk_est[1], ETA2 = pk_est[1], 
                                      ETA3=pk_est[1], ETA4=pk_est[1], 
                                      ETA5 =pk_est[1], ETA6 = pk_est[1], 
                                      beta_scr_cl = mp$beta_scr_cl, beta_bsa_cl = mp$beta_bsa_cl, beta_gender= mp$beta_gender) %>% 
    dplyr::select(-c(A1,A2, A3) )
  
  out_pop_mtxpk <- simulate_mtxpk(dat_subj, mp_mtxpk$THETA1, mp_mtxpk$THETA2, mp_mtxpk$THETA3, mp_mtxpk$THETA4, mp_mtxpk$THETA5, mp_mtxpk$THETA6, 
                                         ETA1 = 0, ETA4=0, ETA5 = 0, ETA6 = 0, 
                                         beta_scr_cl = mp$beta_scr_cl) %>% 
    dplyr::select(-c(A1,A2, A3)) 
  
  IPRED_mtxpk = simulate_mtxpk(dat_subj, mp_mtxpk$THETA1, mp_mtxpk$THETA2, mp_mtxpk$THETA3, mp_mtxpk$THETA4, mp_mtxpk$THETA5, mp_mtxpk$THETA6, 
                                      ETA1 = pk_est_mtxpk[1], ETA4=pk_est_mtxpk[1], 
                                      ETA5 =pk_est_mtxpk[1], ETA6 = pk_est_mtxpk[1], 
                                      beta_scr_cl = mp$beta_scr_cl) %>% 
    dplyr::select(-c(A1,A2, A3) )
  
  out_pop$tag <- "Population level: Blackman et al."
  IPRED$tag <- "Individual level: Blackman et al."
  dat_subj$tag <- "Observed Drug Level"
  out_pop_lower$tag <- "Lower 95% population probability: Blackman et al."
  out_pop_upper$tag <- "Upper 95% population probability: Blackman et al."
  out_pop_mtxpk$tag <- "Population level: MTXPK.org"
  IPRED_mtxpk$tag <- "Individual level: MTXPK.org"
  
  
  
  plt <- rbind(out_pop, IPRED, out_pop_lower, out_pop_upper, out_pop_mtxpk, IPRED_mtxpk)
  
  p <- ggplot(data=plt, aes(x=time,y=Conc, color=tag)) + geom_line() + 
    geom_point(data=subset(dat_subj, is.na(amt)), aes(x=time, y=dv)) 
  return(list(plt, p))
  
}

##########################
######### set up model ##########
##########################
setupModel <- function(in_infusmat1 = in_infusmat1, in_concmat = in_concmat, drug, bsa, scr_vals,  scr_times , pt_gender) { #... is everything not named from my params
  
  infN <- 1 #number of infs 1
  conN <- 1
  allN <- max(infN, conN)
  
  # check for consistent size and columns
  # check for required columns
  df <- in_infusmat1 #infusion in list. df will have dose time and dur
  ndf <- names(df)
  if(!('dose' %in% ndf)) {
    stop('infusion data should contain "dose" column')
  }
  if(!('time' %in% ndf)) {
    dtvar <- which(vapply(df, inherits, logical(1), 'POSIXt'))[1]
    if(is.na(dtvar)) stop('infusion data should contain date-time variable or "time" column')
    df[,'time'] <- unclass(df[,dtvar]) / 3600 # convert from seconds to hours
  }
  if(!('duration' %in% ndf)) { #sets a default for duration
    df[,'duration'] <- 2 # default: 2 hours
  }
  df[is.na(df[,'duration']),'duration'] <- 2
  df <- df[,c('dose','time','duration')]
  in_infusmat1 <- df[complete.cases(df),]
  
  df <- in_concmat
  ndf <- names(df)
  if(!('dose' %in% ndf)) {
    stop('concentration data should contain "dose" column')
  }
  if(!('time' %in% ndf)) {
    dtvar <- which(vapply(df, inherits, logical(1), 'POSIXt'))[1]
    if(is.na(dtvar)) stop('concentration data should contain date-time variable or "time" column')
    df[,'time'] <- unclass(df[,dtvar]) / 3600 # convert from seconds to hours
  }
  df <- df[,c('dose','time')]
  in_concmat <- df[complete.cases(df),]
  
  
  time0 <- min(in_infusmat1$time, in_concmat$time)
  in_infusmat1$time = (in_infusmat1$time - time0)/3600
  in_concmat$time = (in_concmat$time - time0)/3600
  
  # data check for SCR values
  if(is.na(scr_times[1])){
    scr_times[1] = 0
  }
  scr_vals = scr_vals[!is.na(scr_times)]

  scr_times = scr_times[!is.na(scr_times)]
  scr_times = round((scr_times - time0)/3600,1)
  last.t <- max(in_infusmat1$time, in_concmat$time, scr_times, na.rm=TRUE)
  usrendt <- max(last.t + 12, 72) #max for input
  # check covariates
  negNA <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[x < 0] <- NA
    x
  }
  usrgender <- usrSCR <- usrbsa <- NA

  usrgender = pt_gender
  usrSCR = scr_vals*88.4
  usrbsa = bsa
  mp <- model_mtx(usrbsa = usrbsa, usrSCR = usrSCR, usrgender = usrgender)
  mp_mtxpk <- model_mtxpk(usrbsa = usrbsa, usrSCR = usrSCR)
  
  # construct dataframe
  timevec = seq(from = 0, to = usrendt, by =0.1)
  time_and_covs = data.frame(time = timevec, amt = NA, dv= NA, rate = NA, 
                             mdv = 1, evid = 0, dur = NA,
                             pt_gender = mp$pt_gender, bsa = mp$bsa, SCR_mmol = NA, SCR_mgdl = NA)

  for (i in 1:length(scr_times)){
    row = which(abs(time_and_covs$time - scr_times[i]) < 0.001)
    time_and_covs$SCR_mmol[row] = usrSCR[i]
    time_and_covs$SCR_mgdl[row] = usrSCR[i]/88.4
  }

  
  inf = in_infusmat1
  for (i in 1:nrow(inf)){
    time = inf$time[i]
    dose = inf$dose[i]
    duration = inf$duration[i]
    rate = dose/duration
    idx = which(abs(time_and_covs$time - round(time,1)) < 0.01)
    time_and_covs$amt[idx] = dose ; time_and_covs$dur[idx] = duration ; time_and_covs$rate[idx] = rate
    time_and_covs$evid[idx] = 1
  }
  conc = in_concmat
  for (i in 1:nrow(conc)){
    row = i
    time = conc$time[i]
    dv = conc$dose[i]
    idx = which(abs(time_and_covs$time - round(time,1)) < 0.01)
    time_and_covs$dv[idx] = dv ; time_and_covs$mdv[idx] = 0
  }
  
  library(tidyverse)
  time_and_covs = time_and_covs %>% 
    fill(dur,SCR_mmol, SCR_mgdl, .direction = "down") 
  
  dat_subj = subset(time_and_covs, !is.na(dv) | !is.na(amt))

  
  def.par <- list(
    y = conc[,'dose'], yt = conc[,'time']
  )
  
  pk_est <- pkprof_est_mymodel(time_and_covs, mp)
  print("pk_est my model ran")
  pk_est_mtxpk <- pkprof_est_mtxpk(time_and_covs, mp = mp_mtxpk)
  print("pk_est mtxpk ran")
  
  preds  = getpatdat(dat_subj = time_and_covs, mp, pk_est, mp_mtxpk, pk_est_mtxpk)

  concdat <- data.frame(y = in_concmat[,'dose'], yt = in_concmat[,'time'], colour = " Observed Drug Level")

  return(list(concdat, preds[[1]], preds[[2]], mp))
}

makePlots <- function(schedule1, concdat) {
  level_order <- c(
    "Observed Drug Level",              # Col 1, Row 1
    " ",                                # Col 2, Row 1 (The Blank)
    "Serum Creatinine (mg/L)",          # Col 1, Row 2
    "Glucarpidase Consensus Guidelines", # Col 2, Row 2
    "Population level: Blackman et al.",
    "Individual level: Blackman et al.",
    "Population level: MTXPK.org",
    "Individual level: MTXPK.org",
    "Upper 95% population probability: Blackman et al.",
    "Lower 95% population probability: Blackman et al."
  )
  
  cols <- c(
    "Observed Drug Level" = "red",
    " " = "rgba(0,0,0,0)", # Using transparent for Plotly compatibility
    "Serum Creatinine (mg/L)" = "purple",
    "Glucarpidase Consensus Guidelines" = "blue2",
    "Population level: Blackman et al." = "deeppink4",
    "Individual level: Blackman et al." = "deeppink1",
    "Population level: MTXPK.org" = "darkorange1",
    "Individual level: MTXPK.org" = "darkgoldenrod2",
    "Upper 95% population probability: Blackman et al." = "grey50",
    "Lower 95% population probability: Blackman et al." = "grey50"
  )
  # level_order <- c(
  #   "Observed Drug Level",
  #   "Serum Creatinine (mg/L)",
  #   "Glucarpidase Consensus Guidelines",
  #   "Population level: Blackman et al.",
  #   "Individual level: Blackman et al.",
  #   "Population level: MTXPK.org",
  #   "Individual level: MTXPK.org",
  #   "Upper 95% population probability: Blackman et al.",
  #   "Lower 95% population probability: Blackman et al."
  # )
  # 
  # cols <- c(
  #   "Observed Drug Level" = "red",
  #   "Serum Creatinine (mg/L)" = "purple",
  #   "Glucarpidase Consensus Guidelines" = "blue2",
  #   "Population level: Blackman et al." = "deeppink4",
  #   "Individual level: Blackman et al." = "deeppink1",
  #   "Population level: MTXPK.org" = "darkorange1",
  #   "Individual level: MTXPK.org" = "darkgoldenrod2",
  #   "Upper 95% population probability: Blackman et al." = "grey50",
  #   "Lower 95% population probability: Blackman et al." = "grey50"
  # )
  # 
  # Prep datasets with factors to enforce order in Plotly
  gluc_guidelines <- data.frame(
    time = c(24, 36, 42, 48), 
    amt = c(50, 30, 10, 5) * 0.454,
    tag = factor("Glucarpidase Consensus Guidelines", levels = level_order)
  ) 
  
  df_ribbon <- schedule1 %>%
    filter(tag %in% c("Lower 95% population probability: Blackman et al.", 
                      "Upper 95% population probability: Blackman et al.")) %>%
    dplyr::select(time, tag, Conc) %>% 
    pivot_wider(names_from = tag, values_from = Conc) %>% 
    rename(lower = 2, upper = 3)
  
  mult <- 30
  scr_plot_data <- schedule1 %>%
    filter(SCR_mmol != lag(SCR_mmol) | is.na(lag(SCR_mmol))) %>% 
    mutate(SCR_mgdl = SCR_mmol / 88.4,
           y_scaled = SCR_mgdl * mult,
           tag = factor("Serum Creatinine (mg/L)", levels = level_order))
  
  # Ensure the main data is factored
  schedule1$tag <- factor(schedule1$tag, levels = level_order)
  
  ggplot() +
    # Ribbon first (no legend name for ribbon usually)
    geom_ribbon(data = df_ribbon, aes(x = time, ymin = lower, ymax = upper), 
                fill = "grey10", alpha = 0.15) +
    
    # Lines - use 'name = tag' for Plotly
    geom_line(data = schedule1 %>% filter(tag %in% level_order[4:10]), 
              aes(x = time, y = Conc, color = tag, group = tag, name = tag)) +
    
    # Observed Points
    geom_point(data = concdat, 
               aes(x = yt, y = y, color = factor("Observed Drug Level", levels = level_order), 
                   name = factor("Observed Drug Level", levels = level_order)), 
               size = 2.5) +
    
    # Glucarpidase Points
    geom_point(data = gluc_guidelines, aes(x = time, y = amt, color = tag, name = tag), 
               shape = 5, size = 2.5, stroke = 0.6) +
    
    # Serum Creatinine Points
    geom_point(data = scr_plot_data, aes(x = time, y = y_scaled, color = tag, name = tag), 
               size = 3, shape = 2, stroke = 0.6) +
    geom_point(aes(x = 0, y = -1e6, color = factor(" ", levels = level_order)), 
               alpha = 0, show.legend = TRUE) +
    scale_colour_manual(name = NULL, values = cols, breaks = level_order, drop = FALSE) +
    # scale_colour_manual(name = NULL, values = cols, breaks = level_order) +
     scale_x_continuous(name = "Time (hours)", breaks = seq(0, 72, by = 12)) +
     scale_y_continuous(name = "Concentration (mg/L)") 
}
#makePlots <- function(schedule1, concdat) {
  
#   # 1. Define the specific order for the legend
#   level_order <- c(
#     "Observed Drug Level",
#     "Serum Creatinine (mg/L)",
#     "Glucarpidase Consensus Guidelines",
#     "Population level: Blackman et al.",
#     "Individual level: Blackman et al.",
#     "Population level: MTXPK.org",
#     "Individual level: MTXPK.org",
#     "Upper 95% population probability: Blackman et al.",
#     "Lower 95% population probability: Blackman et al."
#   )
#   
#   # 2. Define Colors
#   cols <- c(
#     "Observed Drug Level" = "red",
#     "Serum Creatinine (mg/L)" = "purple",
#     "Glucarpidase Consensus Guidelines" = "blue2",
#     "Population level: Blackman et al." = "deeppink4",
#     "Individual level: Blackman et al." = "deeppink1",
#     "Population level: MTXPK.org" = "darkorange1",
#     "Individual level: MTXPK.org" = "darkgoldenrod2",
#     "Upper 95% population probability: Blackman et al." = "grey50",
#     "Lower 95% population probability: Blackman et al." = "grey50"
#   )
#   
#   line_widths <- c(
#     "Population level: Blackman et al." = 1,
#     "Individual level: Blackman et al." = 1,
#     "Population level: MTXPK.org" = 0.7,
#     "Individual level: MTXPK.org" = 0.7,
#     "Observed Drug Level" = 0.7,
#     "Serum Creatinine (mg/L)" = 0.7,
#     "Glucarpidase Consensus Guidelines" = 0.7
#   )
#   
#   # Explicitly setting 1.0 for everything EXCEPT the CI lines
#   # alpha_values <- c(
#   #   "Observed Drug Level" = 1.0,
#   #   "Serum Creatinine (mg/L)" = 1.0,
#   #   "Glucarpidase Consensus Guidelines" = 1.0,
#   #   "Population level: Blackman et al." = 1.0,
#   #   "Individual level: Blackman et al." = 1.0,
#   #   "Population level: MTXPK.org" = 1.0,
#   #   "Individual level: MTXPK.org" = 1.0,
#   #   "Upper 95% population probability: Blackman et al." = 0.3, 
#   #   "Lower 95% population probability: Blackman et al." = 0.3
#   # )
#   
#   # 3. Data Prep: Glucarpidase Reference Points
#   # (Your schedule1 logic for adding the guideline tag)
#   gluc_guidelines <- data.frame(
#     time = c(24, 36, 42, 48), 
#     amt = c(50, 30, 10, 5) * 0.454,
#     tag = "Glucarpidase Consensus Guidelines"
#   ) 
#   
#   # 4. Data Prep: Ribbon for Blackman CIs
#   df_ribbon <- schedule1 %>%
#     filter(tag %in% c("Lower 95% population probability: Blackman et al.", 
#                       "Upper 95% population probability: Blackman et al.")) %>%
#     select(time, tag, Conc) %>% 
#     pivot_wider(names_from = tag, values_from = Conc) %>% 
#     rename(lower = "Lower 95% population probability: Blackman et al.",
#            upper = "Upper 95% population probability: Blackman et al.")
#   
#   # 5. Data Prep: Serum Creatinine (SCr)
#   # Converting mmol to mg/dL and filtering for changes
#   mult <- 30
#   scr_plot_data <- schedule1 %>%
#     filter(SCR_mmol != lag(SCR_mmol) | is.na(lag(SCR_mmol))) %>% 
#     mutate(SCR_mgdl = SCR_mmol / 88.4,
#            y_scaled = SCR_mgdl * mult,
#            tag = "Serum Creatinine (mg/L)")
#   
#   # 6. Final Factor logic for the main data
#   schedule1 <- schedule1 %>%
#     mutate(tag = factor(tag, levels = level_order))
#   
#   # 7. Construct Plot
#   ggplot() +
#     # Ribbon (Background)
#     geom_ribbon(data = df_ribbon, aes(x = time, ymin = lower, ymax = upper), 
#                 fill = "grey10", alpha = 0.15) +
#     
#     # Model Lines (Population/Individual)
#     geom_line(data = schedule1 %>% filter(!tag %in% c("Serum Creatinine (mg/L)", "Observed Drug Level")), 
#               aes(x = time, y = Conc, color = tag, group = tag, linewidth = tag)) +
#     scale_linewidth_manual(values = line_widths, guide = "none") +
#     # Observed Data Points
#     geom_point(data = concdat, aes(x = yt, y = y, color = "Observed Drug Level"), size = 2.5) +
#     
#     # Glucarpidase Points
#     geom_point(data = gluc_guidelines, aes(x = time, y = amt, color = tag), 
#                shape = 5, size = 2.5, stroke = 0.7) +
#     
#     # Serum Creatinine Points
#     geom_point(data = scr_plot_data, aes(x = time, y = y_scaled, color = tag), 
#                size = 3, shape = 2, stroke = 0.7) +
#     
#     # Scales and Axes
#     scale_colour_manual(name = NULL, values = cols, breaks = level_order) +
#     scale_x_continuous(name = "Time (hours)", breaks = seq(0, 72, by = 12)) +
#     scale_y_continuous(
#       name = "Concentration (mg/L)",
#       sec.axis = sec_axis(~ . / mult, name = "Serum Creatinine (mg/dL)")
#     ) +
#     theme_minimal() +
#     theme(
#       legend.position = "bottom",
#       legend.title = element_blank(),
#       axis.title.y.right = element_text(color = "purple"),
#       axis.text.y.right = element_text(color = "purple")
#     ) +
#     # THIS SECTION FIXES THE LEGEND LABELS
#     guides(color = guide_legend(
#       ncol = 2,
#       override.aes = list(
#         linewidth = 1,   # Force a standard line thickness in the legend
#         alpha = 1        # Force everything to be opaque in the legend
#       )
#     )) # Stack legend for readability
# }

# bldplot <- function(dat) {
# 
#   mx <- max(dat[,'time'], na.rm = TRUE)
#   new = data.frame(time = c(1), bsa = NA, SCR_mmol = NA, pt_gender = NA, 
#                    CL = NA, V1 = NA, Q2 = NA, V2 = NA, Q3 = NA, V3 = NA,
#                    Conc = c(-3 ), evid = NA, tag = c("Glucarpidase Consensus Guidelines" ))
#   dat <- dat %>%
#     # rename(Legend = tag) %>%
#     bind_rows(new) %>%
#     #arrange(tag)%>%
#     mutate(
#       tag = factor(
#         tag,
#         levels = c(         "Observed Drug Level",
#                             "Population level: Blackman et al.",
#                              "Individual level: Blackman et al.",
#                              "Population level: MTXPK.org",
#                              "Individual level: MTXPK.org",
#                              "Upper 95% population probability: Blackman et al.",
#                              "Lower 95% population probability: Blackman et al.",
#           "Serum Creatinine (mg/L)",
#           "Glucarpidase Consensus Guidelines"
#         ), ordered = TRUE
#       ),
#       Legend = case_when(tag == "Population level: Blackman et al." ~ 1,
#                             tag == "Individual level: Blackman et al." ~ 1,
#                             tag ==  "Population level: MTXPK.org"~ 1,
#                             tag == "Individual level: MTXPK.org" ~ 1,
#                             tag ==  "Lower 95% population probability: Blackman et al."~ 1,
#                             tag == "Upper 95% population probability: Blackman et al." ~ 1,
#                             tag ==  "Observed Drug Level" ~ 1,
#                             tag == "Serum Creatinine (mg/L)" ~ 1,
#                             tag ==  "Glucarpidase Consensus Guidelines"~ 1)
#     )
#   
#   ggplot(dat) +
#     geom_line(
#       aes(x = time, y = Conc, color = tag, group = tag, alpha = Legend))#+
#     # scale_colour_manual(
#     #   values = cols,
#     #   limits = c(
#     #     "Observed Drug Level",
#     #     "Population level: Blackman et al.",
#     #     "Individual level: Blackman et al.",
#     #     "Population level: MTXPK.org",
#     #     "Individual level: MTXPK.org",
#     #     "Upper 95% population probability: Blackman et al.",
#     #     "Lower 95% population probability: Blackman et al.",
#     #     "Serum Creatinine (mg/L)",
#     #     "Glucarpidase Consensus Guidelines"
#     #   ),
#     #   drop = FALSE
#     # )
#   
# }
# 
# makePlots <- function(schedule1,  concdat) {
#   
#   gluc = data.frame(time = c(24, 36, 42, 48), amt = c(50, 30, 10, 5)*0.454 ) 
#   schedule1 <- schedule1 %>%
#     mutate(
#       tag = factor(
#         tag,
#         levels = c(
#           "Observed Drug Level",
#           "Population level: Blackman et al.",
#           "Individual level: Blackman et al.",
#           "Population level: MTXPK.org",
#           "Individual level: MTXPK.org",
#           "Upper 95% population probability: Blackman et al.",
#           "Lower 95% population probability: Blackman et al.",
#           "Serum Creatinine (mg/L)",
#           "Glucarpidase Consensus Guidelines"
#         )
#       )
#     )
#   
#   cols <- c(
#     "Observed Drug Level" = "red",
#     "Serum Creatinine (mg/L)" = "purple",
#     "Glucarpidase Consensus Guidelines" = "blue2",
#     "Population level: Blackman et al." = "deeppink4",
#     "Individual level: Blackman et al." = "deeppink1",
#     "Population level: MTXPK.org" = "darkorange1",
#     "Individual level: MTXPK.org" = "darkgoldenrod2",
#     "Lower 95% population probability: Blackman et al." = "grey50",
#     "Upper 95% population probability: Blackman et al." = "grey50"
#   )
#   
#   df_ribbon <- schedule1 %>%
#     filter(tag %in% c("Lower 95% population probability: Blackman et al.", "Upper 95% population probability: Blackman et al.")) %>%
#     dplyr::select(time, tag, Conc) %>% 
#     pivot_wider(
#       names_from  = tag,
#       values_from = Conc
#     )  %>% 
#     rename(lower = "Lower 95% population probability: Blackman et al.",
#            upper = "Upper 95% population probability: Blackman et al.")
#   
#   scr = schedule1 %>%
#     filter(SCR_mmol != lag(SCR_mmol) | is.na(lag(SCR_mmol))) %>% 
#     mutate(SCR_mgdl = SCR_mmol/88.4)
#   
#   mult = 30
#     
#   #y_max <- max(c(schedule1$Conc, concdat$y, gluc$amt, scr$SCR_mgdl * mult, df_ribbon$upper), na.rm = TRUE)
#   print(levels(schedule1$tag))
#   p1 <- bldplot(schedule1) +
#     geom_point(data = concdat, aes(x = yt, y = y, color = "Observed Drug Level"), size  = 2 ) +
#   #   geom_point(data = concdat %>% mutate(tag = "Observed Drug Level"),
#   #              aes(x = yt, y = y, color = tag))+
#   #   geom_point(data = gluc %>% mutate(tag = "Glucarpidase Consensus Guidelines"),
#   #              aes(x = time, y = amt, color = tag))+
#   # 
#   # geom_point(data = scr %>% mutate(tag = "Serum Creatinine (mg/L)"),
#   #            aes(x = time, y = SCR_mgdl * mult, color = tag))+
#     # Glucarpidase dosing
#     geom_point(data = gluc, aes(x = time, y = amt, color = "Glucarpidase Consensus Guidelines"),shape = 5, size  = 2.5) +
#     # Population ribbon
#     geom_ribbon(data = df_ribbon, aes(x = time, ymin = lower, ymax = upper), fill = "grey10",alpha = 0.3) +
#     # Serum creatinine → POINTS, SECOND AXIS
#     geom_point(data = scr, aes(x = time, y = SCR_mgdl * mult, color = "Serum Creatinine (mg/L)"),size = 3, shape =2) +
#     
#     scale_colour_manual(
#       values = cols,
#       breaks = levels(schedule1$tag)
#     ) +
#     scale_x_continuous(
#       name = "Time (hours)",
#       breaks = seq(0, 72, by = 12)
#     ) +
#     
#     scale_y_continuous(
#       name = "Concentration (mg/L)",
#       sec.axis = sec_axis(~ . / 30, name = "Serum Creatinine (mg/L)")
#     ) +
#     # theme(
#     #   axis.title.y.right = element_text(angle = 90, margin = margin(l = 15), vjust = 1),
#     #   plot.margin = margin(10, 60, 10, 10) 
#     # ) +
#     labs(color = "")
#   
# 
# }


arrange_data <- function(df) {
  dfn <- tolower(names(df))
  rateCol <- grep('rate', dfn)
  amntCol <- grep('amt', dfn)
  concCol <- grep('conc', dfn)
  duraCol <- grep('dur', dfn)
  timeCol <- grep('time', dfn)
  timeVar <- 'time'
  if(length(rateCol) == 0) stop('rate column should exist')
  if(length(amntCol) == 0) stop('amt column should exist')
  if(length(concCol) == 0) stop('conc column should exist')
  if(length(timeCol) == 0) {
    addlCol <- setdiff(seq(ncol(df)), c(rateCol, amntCol, concCol, duraCol, timeCol))
    dtCol <- character(0)
    # is additional column a date-time variable?
    if(length(addlCol) > 0) {
      isDT <- apply(df[,addlCol,drop=FALSE], 2, function(i) tryCatch(pkdata::guessDateFormat(i), error=function(e) NA_character_))
      dtCol <- dfn[addlCol[!is.na(isDT)]]
    }
    if(length(dtCol) != 1) stop('date-time or time column should exist')
    timeVar <- 'datetime'
    df[,timeVar] <- as.POSIXct(df[,dtCol], format = isDT[!is.na(isDT)])
    timeCol <- match(timeVar, names(df))
  }
  df_i <- df[!is.na(df[,rateCol]),c(rateCol,timeCol,duraCol)]
  df_c <- df[!is.na(df[,concCol]),c(concCol,timeCol)]
  names(df_i)[1:2] <- c('dose',timeVar)
  if(ncol(df_i) == 3) names(df_i)[3] <- 'duration'
  names(df_c) <- c('dose',timeVar)
  covar <- list(
    wt = df[1,grep('weight', dfn)],
    alb = df[1,grep('alb', dfn)],
    ada = df[1,grep('ada', dfn)]
  )
  list(infList = list(df_i), conList = list(df_c), covariates = covar)
}
