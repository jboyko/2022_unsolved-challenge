require(corHMM)
require(geiger)
require(parallel)

setwd("2021_unresolved-problem/")

simulateData <- function(ntaxa){
  meets_criteria_phy_index <- integer(0)
  while(length(meets_criteria_phy_index) == 0){
    phy <- sim.bdtree(b= 1, d = 0.5, "taxa", ntaxa, extinct = FALSE)
    phy <- drop.extinct(phy)
    phy <- ladderize(phy)
    subphy <- subtrees(phy)
    meets_criteria_phy_index <- which((unlist(lapply(subphy, function(x) length(x$tip.label))) > (ntaxa * .4)) & (unlist(lapply(subphy, function(x) length(x$tip.label))) < (ntaxa * .6)))
  }
  # unrepliclated bursts 
  # generate character X
  meets_criteria_phy_index <- meets_criteria_phy_index[1]
  dat <- data.frame(sp = phy$tip.label, X = 0, Y = 0)
  dat[dat$sp %in% subphy[[meets_criteria_phy_index]]$tip.label, 2] <- 1
  # generate character Y
  sim_dat <- sim.char(phy, matrix(c(-100,100,100,-100), 2, 2), 1, "discrete")
  dat[dat$sp %in% subphy[[meets_criteria_phy_index]]$tip.label, 3] <- sim_dat[rownames(sim_dat) %in% subphy[[meets_criteria_phy_index]]$tip.label] - 1
  burst_dat <- dat
  # darwins scenario
  # generate character X
  # generate character Y
  dat <- data.frame(sp = phy$tip.label, X = 0, Y = 0)
  dat[dat$sp %in% subphy[[meets_criteria_phy_index]]$tip.label, c(2,3)] <- 1
  darwins_dat_both <- darwins_dat_outside <- darwins_dat_inside <- darwins_dat <- dat
  # set up the modified scenarios
  # outside focal clade
  tips.to.switch.00 <- sample(which(paste0(darwins_dat$X, darwins_dat$Y) == "00"), 2)
  darwins_dat_outside[tips.to.switch.00[1], 2:3] <- c(1,0)
  darwins_dat_outside[tips.to.switch.00[2], 2:3] <- c(0,1)
  # inside focal clade
  tips.to.switch.11 <- sample(which(paste0(burst_dat$X, burst_dat$Y) == "11"), 2)
  darwins_dat_inside[tips.to.switch.11[1], 2:3] <- c(1,0)
  darwins_dat_inside[tips.to.switch.11[2], 2:3] <- c(0,1)
  # both
  darwins_dat_both[tips.to.switch.00[1], 2:3] <- c(1,0)
  darwins_dat_both[tips.to.switch.00[2], 2:3] <- c(0,1)
  darwins_dat_both[tips.to.switch.11[1], 2:3] <- c(1,0)
  darwins_dat_both[tips.to.switch.11[2], 2:3] <- c(0,1)
  return(list(phy=phy, darwin_dat=darwins_dat, burst_dat=burst_dat, darwins_dat_outside=darwins_dat_outside, darwins_dat_inside=darwins_dat_inside, darwins_dat_both=darwins_dat_both))
}

fit_model_set <- function(phy, dat, collapsed=TRUE){
  rate.mat <- getFullMat(list(getRateCatMat(2),getRateCatMat(2)), getRateCatMat(2))
  colnames(rate.mat) <- rownames(rate.mat) <- c("0,0", "0,1", "1,0", "1,1")
  rate.mat[rate.mat != 0] <- 1:8
  # independent models
  rate_mat_ind <- equateStateMatPars(rate.mat, list(c(3,8), c(5,7), c(1,6), c(2,4)))
  rate_mat_ind_2 <- getFullMat(list(rate_mat_ind, rate_mat_ind), equateStateMatPars(getRateCatMat(2), 1:2))
  rate_mat_simp_ind <- equateStateMatPars(rate_mat_ind, list(c(1,2), c(3,4)))
  rate_mat_simp_ind_2 <- getFullMat(list(rate_mat_simp_ind, rate_mat_simp_ind), equateStateMatPars(getRateCatMat(2), 1:2))
  # correlated models
  rate_mat_cor <- rate.mat
  rate_mat_cor_2 <- getFullMat(list(rate_mat_cor, rate_mat_cor), equateStateMatPars(getRateCatMat(2), 1:2))
  rate_mat_simp_cor <- equateStateMatPars(rate_mat_cor, list(c(1,3,2,5), c(4,7,6,8)))
  rate_mat_simp_cor_2 <- getFullMat(list(rate_mat_simp_cor, rate_mat_simp_cor), equateStateMatPars(getRateCatMat(2), 1:2))
  # fit single rate models
  fit_independent <- corHMM(phy, data = dat, rate.cat = 1, rate.mat = rate_mat_ind, collapse = FALSE)
  fit_simple_independent <- corHMM(phy, data = dat, rate.cat = 1, rate.mat = rate_mat_simp_ind, collapse = FALSE)
  fit_correlated <- corHMM(phy, data = dat, rate.cat = 1, rate.mat = rate_mat_cor, collapse = FALSE)
  fit_simple_correlated <- corHMM(phy, data = dat, rate.cat = 1, rate.mat = rate_mat_simp_cor, collapse = FALSE)
  # fit hidden state models
  fit_independent_2 <- corHMM(phy, data = dat, rate.cat = 2, rate.mat = rate_mat_ind_2, collapse = FALSE)
  fit_simple_independent_2 <- corHMM(phy, data = dat, rate.cat = 2, rate.mat = rate_mat_simp_ind_2, collapse = FALSE)
  fit_correlated_2 <- corHMM(phy, data = dat, rate.cat = 2, rate.mat = rate_mat_cor_2, collapse = FALSE)
  fit_simple_correlated_2 <- corHMM(phy, data = dat, rate.cat = 2, rate.mat = rate_mat_simp_cor_2, collapse = FALSE)
  if(collapsed){
    collapsed_mat <- getStateMat4Dat(dat)$rate.mat
    if(dim(collapsed_mat)[1] == 2){
      collapsed_mat <- getRateCatMat(2)
    }
    fit_collapsed <- corHMM(phy, data = dat, rate.cat = 1, rate.mat=collapsed_mat) 
    res <- list(fit_collapsed=fit_collapsed, fit_independent=fit_independent, fit_simple_independent=fit_simple_independent, fit_correlated=fit_correlated,fit_simple_correlated=fit_simple_correlated,fit_independent_2=fit_independent_2,fit_simple_independent_2=fit_simple_independent_2,fit_correlated_2=fit_correlated_2,fit_simple_correlated_2=fit_simple_correlated_2)
  }else{
    res <- list(fit_independent=fit_independent, fit_simple_independent=fit_simple_independent, fit_correlated=fit_correlated,fit_simple_correlated=fit_simple_correlated,fit_independent_2=fit_independent_2,fit_simple_independent_2=fit_simple_independent_2,fit_correlated_2=fit_correlated_2,fit_simple_correlated_2=fit_simple_correlated_2)
  }
  return(res)
}

fit_data <- function(dat_list){
  # darwins
  darwin_res <- fit_model_set(dat_list$phy, dat_list$darwin_dat, TRUE)
  # bursts
  burst_res <- fit_model_set(dat_list$phy, dat_list$burst_dat, TRUE)
  # darwins modified scenario outside
  outside_res <- fit_model_set(dat_list$phy, dat_list$darwins_dat_outside, FALSE)
  # darwins modified scenario inside
  inside_res <- fit_model_set(dat_list$phy, dat_list$darwins_dat_inside, FALSE)
  # darwins modified scenario both
  both_res <- fit_model_set(dat_list$phy, dat_list$darwins_dat_both, FALSE)
  out <- list(burst_res=burst_res, darwin_res=darwin_res, outside_res=outside_res, inside_res=inside_res, both_res=both_res)
  return(out)
}

# generate dataset of 100 taaxa
dat_list <- mclapply(1:100, function(x) simulateData(100), mc.cores=50)

# fit the models to the 100 datasets
maddfitz_res <- mclapply(dat_list, fit_data, mc.cores=50)

