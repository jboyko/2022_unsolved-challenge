require(corHMM)
require(geiger)
require(parallel)
require(ggplotify)
require(phytools)
require(ggplot2)
require(gridExtra)

putRatesinPlace <- function(rates, index.mat){
  index.mat[index.mat==0] <- NA
  Q <- index.mat
  for(i in 1:length(index.mat)){
    Q[i] <- rates[index.mat[i]]
  }
  Q[is.na(Q)] <- 0
  diag(Q) <- -rowSums(Q)
  return(Q)
}

ensureSampling <- function(phy, Q, type){
  meets_criteria_data <- FALSE
  while(!meets_criteria_data){
    sim_dat <- sim.char(phy, Q, 1, "discrete", root = 1)[,,1]
    if(type=="hmm"){
      count <- length(which(sim_dat == 5 | sim_dat == 6 | sim_dat == 7 | sim_dat == 8))
      meets_criteria_data <- length(unique(sim_dat)) == dim(Q)[1] & (count > 0.2*length(phy$tip.label) & count < 0.8*length(phy$tip.label))
    }else{
      meets_criteria_data <- (all(table(sim_dat) >= (0.05 * length(phy$tip.label)))) & (dim(Q)[1] == length(unique(sim_dat)))
    }
  }
  return(sim_dat)
}

organizeHMMData <- function(sim_data){
  sim_data[sim_data == 5] <- 1
  sim_data[sim_data == 6] <- 2
  sim_data[sim_data == 7] <- 3
  sim_data[sim_data == 8] <- 4
  return(sim_data)
}

simulateData <- function(rates=1:5, nTaxa=100){
  phy <- sim.bdtree(b= 1, d = 0.75, "taxa", nTaxa, extinct = FALSE)
  phy <- drop.extinct(phy)
  phy <- ladderize(phy)
  phy$edge.length <- phy$edge.length/sum(phy$edge.length) * 10
  rate.mat <- getFullMat(list(getRateCatMat(2),getRateCatMat(2)), getRateCatMat(2))
  colnames(rate.mat) <- rownames(rate.mat) <- c("0,0", "0,1", "1,0", "1,1")
  rate.mat[rate.mat != 0] <- 1:8
  # independent models
  rate_mat_ind <- equateStateMatPars(rate.mat, list(c(3,8), c(5,7), c(1,6), c(2,4)))
  rate_mat_simp_ind <- equateStateMatPars(rate_mat_ind, list(c(1,2), c(3,4)))
  rate_mat_simp_ind_2 <- getFullMat(list(rate_mat_simp_ind, rate_mat_simp_ind), equateStateMatPars(getRateCatMat(2), 1:2))
  # correlated models
  rate_mat_cor <- rate.mat
  rate_mat_simp_cor <- equateStateMatPars(rate_mat_cor, list(c(1,3, 2,5), c(4,7,6,8)))
  # rates
  Q_ind <- putRatesinPlace(rates[1:2], rate_mat_simp_ind)
  Q_ind_2 <- putRatesinPlace(rates, rate_mat_simp_ind_2)
  Q_cor <- putRatesinPlace(rates[1:2], rate_mat_simp_cor)
  # simulate 
  sim_dat_ind <- ensureSampling(phy, Q_ind, "not_hmm")
  sim_dat_cor <- ensureSampling(phy, Q_cor, "not_hmm")
  sim_dat_ind_2 <- ensureSampling(phy, Q_ind_2, "hmm")

  dat_ind <- data.frame(sp = phy$tip.label, 
                        X = as.numeric(sim_dat_ind == 3 | sim_dat_ind == 4), 
                        Y = as.numeric(sim_dat_ind == 2 | sim_dat_ind == 4))
  dat_cor <- data.frame(sp = phy$tip.label, 
                        X = as.numeric(sim_dat_cor == 3 | sim_dat_cor == 4), 
                        Y = as.numeric(sim_dat_cor == 2 | sim_dat_cor == 4))
  dat_ind_2 <- data.frame(sp = phy$tip.label, 
                          X = as.numeric(sim_dat_ind_2 == 3 | sim_dat_ind_2 == 4 | sim_dat_ind_2 == 7, sim_dat_ind_2 == 8), 
                          Y = as.numeric(sim_dat_ind_2 == 2 | sim_dat_ind_2 == 4 | sim_dat_ind_2 == 6, sim_dat_ind_2 == 8))

  return(list(phy=phy, dat_ind=dat_ind, dat_cor=dat_cor, dat_ind_2=dat_ind_2))
}

fit_model_set <- function(phy, dat){
  # the models we're fitting
  rate.mat <- getFullMat(list(getRateCatMat(2),getRateCatMat(2)), getRateCatMat(2))
  colnames(rate.mat) <- rownames(rate.mat) <- c("0,0", "0,1", "1,0", "1,1")
  rate.mat[rate.mat != 0] <- 1:8
  # independent models
  rate_mat_ind <- equateStateMatPars(rate.mat, list(c(3,8), c(5,7), c(1,6), c(2,4)))
  rate_mat_simp_ind <- equateStateMatPars(rate_mat_ind, list(c(1,2), c(3,4)))
  rate_mat_simp_ind_2 <- getFullMat(list(rate_mat_simp_ind, rate_mat_simp_ind), equateStateMatPars(getRateCatMat(2), 1:2))
  # correlated models
  rate_mat_cor <- rate.mat
  rate_mat_simp_cor <- equateStateMatPars(rate_mat_cor, list(c(1,3, 2,5), c(4,7,6,8)))
  # single rate
  ind_fit <- corHMM(phy, dat, rate.cat = 1, rate.mat = rate_mat_simp_ind, node.states = "none")
  cor_fit <- corHMM(phy, dat, rate.cat = 1, rate.mat = rate_mat_simp_cor, node.states = "none")
  # two rate
  ind_2_fit <- corHMM(phy, dat, rate.cat = 2, rate.mat = rate_mat_simp_ind_2, node.states = "none")
  # output
  out <- list(ind_fit=ind_fit, cor_fit=cor_fit, ind_2_fit=ind_2_fit)
  return(out)
}

# fit the models to an element of datalist
fit_data_list <- function(dat_list){
  ind_dat_fits <- fit_model_set(dat_list$phy, dat_list$dat_ind)
  cor_dat_fits <- fit_model_set(dat_list$phy, dat_list$dat_cor)
  ind_2_dat_fits <- fit_model_set(dat_list$phy, dat_list$dat_ind_2)
  out <- list(ind_dat_fits=ind_dat_fits, cor_dat_fits=cor_dat_fits, ind_2_dat_fits=ind_2_dat_fits)
  return(out)
}

summarizeListElement <- function(maddfitz_simmed_res_element){
  out_table <- rbind(unlist(lapply(maddfitz_simmed_res_element$ind_dat_fits, "[[", "AIC")),
                     unlist(lapply(maddfitz_simmed_res_element$cor_dat_fits, "[[", "AIC")),
                     unlist(lapply(maddfitz_simmed_res_element$ind_2_dat_fits, "[[", "AIC")))
  rownames(out_table) <- c("ind_dat", "cor_dat", "ind_2_dat")
  return(out_table)
}


# generate dataset of n taxa
nTaxa <- 100
dat_list_list <- mclapply(1:100, function(x) simulateData(c(1,5,2,10,4), nTaxa), mc.cores=50)
maddfitz_simmed_res <- mclapply(dat_list_list, fit_data_list, mc.cores=50)

