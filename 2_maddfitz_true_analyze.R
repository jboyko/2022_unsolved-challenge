require(corHMM)
require(geiger)
require(parallel)
require(ggplot2)
require(reshape2)
require(RColorBrewer)
require(gridExtra)
require(phytools)
require(ggplotify)

setwd("~/2021_unresolved-problem/")

getModelTable <- function(model_list, type){
  focal_model_list <- model_list[names(model_list) == type][[1]]
  AICs <- getAIC(focal_model_list)
  return(AICs)
}

getAICweight <- function(AICcs){
  dAIC <- AICcs - min(AICcs)
  AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
  return(AICwt)
}

getAIC <- function(model_list){
  return(unlist(lapply(model_list, "[[", "AIC")))
}

# run
load("results/maddfitz_res.Rsave")

# summarize the results
darwin_res <- as.data.frame(do.call(rbind, lapply(maddfitz_res, function(x) getModelTable(x, "darwin_res"))))
burst_res <- as.data.frame(do.call(rbind, lapply(maddfitz_res, function(x) getModelTable(x, "burst_res"))))
outside_res <- as.data.frame(do.call(rbind, lapply(maddfitz_res, function(x) getModelTable(x, "outside_res"))))
inside_res <- as.data.frame(do.call(rbind, lapply(maddfitz_res, function(x) getModelTable(x, "inside_res"))))
both_res <- as.data.frame(do.call(rbind, lapply(maddfitz_res, function(x) getModelTable(x, "both_res"))))

####### PART A: strict darwin's and unreplicated burst scenario
####### PART 1: replicating the original maddfitz result
# strict darwin's scenario
darwin_tmp <- cbind(independent=darwin_res$fit_independent, dependent=darwin_res$fit_correlated)
rowMeans(apply(darwin_tmp, 1, getAICweight)) # avg aic weight
rowSums(apply(darwin_tmp, 1, function(x) (x - min(x)) == 0))
aic_dat <- melt(t(apply(darwin_tmp, 1, getAICweight)))
a <- ggplot(aic_dat, aes(x = Var2, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.8, alpha=0.5) +
  coord_cartesian(ylim=c(0, 1)) + 
  ggtitle("a) Strict Darwin's Scenario") +
  theme_bw() +
  ylab("AIC weight") + 
  xlab("") +
  theme(legend.position = "none")
# strict unreplicated bursts scenario
burst_tmp <- cbind(independent=burst_res$fit_independent, dependent=burst_res$fit_correlated)
rowMeans(apply(burst_tmp, 1, getAICweight))
rowSums(apply(burst_tmp, 1, function(x) (x - min(x)) == 0))
aic_dat <- melt(t(apply(burst_tmp, 1, getAICweight)))
b <- ggplot(aic_dat, aes(x = Var2, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.8, alpha=0.5) +
  coord_cartesian(ylim=c(0, 1)) + 
  ggtitle("b) Strict Unreplicated Bursts Scenario") +
  theme_bw() +
  ylab("AIC weight") + 
  xlab("") +
  theme(legend.position = "none")
grid.arrange(a, b, nrow=1)
####### PART 2: adding the collapsed model to the original maddfitz result
# strict darwin's scenario
darwin_tmp <- cbind(independent=darwin_res$fit_independent, dependent=darwin_res$fit_correlated, collapsed=darwin_res$fit_collapsed)
rowMeans(apply(darwin_tmp, 1, getAICweight))
rowSums(apply(darwin_tmp, 1, function(x) (x - min(x)) == 0))
aic_dat <- melt(t(apply(darwin_tmp, 1, getAICweight)))
a <- ggplot(aic_dat, aes(x = Var2, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.8, alpha=0.5) +
  coord_cartesian(ylim=c(0, 1)) + 
  ggtitle("a) Strict Darwin's Scenario") +
  theme_bw() +
  ylab("AIC weight") + 
  xlab("") +
  theme(legend.position = "none")

# strict unreplicated bursts scenario
burst_tmp <- cbind(independent=burst_res$fit_independent, dependent=burst_res$fit_correlated, collapsed=darwin_res$fit_collapsed)
rowMeans(apply(burst_tmp, 1, getAICweight))
rowSums(apply(burst_tmp, 1, function(x) (x - min(x)) == 0))
aic_dat <- melt(t(apply(burst_tmp, 1, getAICweight)))
b <- ggplot(aic_dat, aes(x = Var2, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.8, alpha=0.5) +
  coord_cartesian(ylim=c(0, 1)) + 
  ggtitle("b) Strict Unreplicated Bursts Scenario") +
  theme_bw() +
  ylab("AIC weight") + 
  xlab("") +
  theme(legend.position = "none")
grid.arrange(a, b, nrow=1)

####### PART B: modified darwin's scenario
outside_tmp <- cbind(independent=outside_res$fit_independent, correlated=outside_res$fit_correlated, independent_2=outside_res$fit_independent_2)
outside_tmp <- as.data.frame(t(apply(outside_tmp, 1, getAICweight)))
inside_tmp <- cbind(independent=inside_res$fit_independent, correlated=inside_res$fit_correlated, independent_2=inside_res$fit_independent_2)
inside_tmp <- as.data.frame(t(apply(inside_tmp, 1, getAICweight)))
both_tmp <- cbind(independent=both_res$fit_independent, correlated=both_res$fit_correlated, independent_2=both_res$fit_independent_2)
both_tmp <- as.data.frame(t(apply(both_tmp, 1, getAICweight)))

####### PART 3: replicating maddfitz result with modified darwin's scenario
# outside scenario
exp(sum(log(outside_tmp$correlated/outside_tmp$independent))/100)
# inside scenario
exp(sum(log(inside_tmp$correlated/inside_tmp$independent))/100)
# both_scenario
exp(sum(log(both_tmp$correlated/both_tmp$independent))/100)

####### PART 4: evidence of rate heterogeneity with modified darwin's scenario
# outside scenario
exp(sum(log(outside_tmp$independent_2/outside_tmp$independent))/100)
# inside scenario
exp(sum(log(inside_tmp$independent_2/inside_tmp$independent))/100)
# both_scenario
exp(sum(log(both_tmp$independent_2/both_tmp$independent))/100)

####### PART 5: remaining evidence of correlation
# outside scenario
exp(sum(log(outside_tmp$correlated/outside_tmp$independent_2))/100)
# inside scenario
exp(sum(log(inside_tmp$correlated/inside_tmp$independent_2))/100)
# both_scenario
exp(sum(log(both_tmp$correlated/both_tmp$independent_2))/100)
# both_scenario
exp(sum(log(both_tmp$independent_2/both_tmp$correlated))/100)

####### PART 6: creating the plot for figure 5 of evidence ratio and modified darwins
cols2 <- setNames(c("#d7191c","#fdae61"), c("1","2"))
cols3 <- setNames(c("#2c7bb6", "#abd9e9"), c("1","2"))
lwd <- 5
types <- c("corr_vs_ind", "corr_vs_hmm")
### outside
phy <- maddfitz_res[[21]]$outside_res$fit_independent$phy
data <- maddfitz_res[[21]]$outside_res$fit_independent$data
er_a <- outside_tmp$correlated/outside_tmp$independent
er_b <- outside_tmp$correlated/outside_tmp$independent_2
map_x <- makeSimmap(phy, data[,c(1,2)], matrix(c(-1e-5,1e-5,1e-5,-1e-5), 2, 2), rate.cat = 1)[[1]]
map_y <- makeSimmap(phy, data[,c(1,3)], matrix(c(-1e-5,1e-5,1e-5,-1e-5), 2, 2), rate.cat = 1)[[1]]
p_x <- as.grob(~plotSimmap(map_x, color=cols2, fsize = 1e-10, lwd = lwd))
p_y <- as.grob(~plotSimmap(map_y, color=cols3, fsize = 1e-10, lwd = lwd, direction="leftwards"))
er <- rbind(data.frame(type = types[1], evidence_ratio = er_a), 
            data.frame(type = types[2], evidence_ratio = er_b))
er$type <- factor(er$type, levels = types)
p_out <- ggplot(er, aes(x = type, y = evidence_ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.8, alpha=0.5) +
  coord_cartesian(ylim=c(0, 50)) + 
  ggtitle("a) Modified Darwin's scenario (outside)") +
  theme_bw() +
  ylab("Support for correlation (evidence ratio)") + 
  xlab("") +
  geom_hline(yintercept = 2.7, size=1, linetype="dashed") + 
  theme(legend.position = "none", plot.title = element_text(size = 21), axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 20)) +
  scale_x_discrete(labels = c('IM','HMIM'))
out_plot <- grid.arrange(p_out, p_x, p_y, nrow=1, widths=c(1.2,1,1))
### inside 
phy <- maddfitz_res[[21]]$inside_res$fit_independent$phy
data <- maddfitz_res[[21]]$inside_res$fit_independent$data
er_a <- inside_tmp$correlated/inside_tmp$independent
er_b <- inside_tmp$correlated/inside_tmp$independent_2
map_x <- makeSimmap(phy, data[,c(1,2)], matrix(c(-1e-5,1e-5,1e-5,-1e-5), 2, 2), rate.cat = 1)[[1]]
map_y <- makeSimmap(phy, data[,c(1,3)], matrix(c(-1e-5,1e-5,1e-5,-1e-5), 2, 2), rate.cat = 1)[[1]]
p_x <- as.grob(~plotSimmap(map_x, color=cols2, fsize = 1e-10, lwd = lwd))
p_y <- as.grob(~plotSimmap(map_y, color=cols3, fsize = 1e-10, lwd = lwd, direction="leftwards"))
er <- rbind(data.frame(type = types[1], evidence_ratio = er_a), 
            data.frame(type = types[2], evidence_ratio = er_b))
er$type <- factor(er$type, levels = types)
p_out <- ggplot(er, aes(x = type, y = evidence_ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.8, alpha=0.5) +
  coord_cartesian(ylim=c(0, 50)) + 
  ggtitle("b) Modified Darwin's scenario (inside)") +
  theme_bw() +
  ylab("Support for correlation (evidence ratio)") + 
  xlab("") +
  geom_hline(yintercept = 2.7, size=1, linetype="dashed") + 
  theme(legend.position = "none", plot.title = element_text(size = 21), axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 20)) +
  scale_x_discrete(labels = c('IM','HMIM'))
ins_plot <- grid.arrange(p_out, p_x, p_y, nrow=1, widths=c(1.2,1,1))
### both 
phy <- maddfitz_res[[21]]$both_res$fit_independent$phy
data <- maddfitz_res[[21]]$both_res$fit_independent$data
er_a <- both_tmp$correlated/both_tmp$independent
er_b <- both_tmp$correlated/both_tmp$independent_2
map_x <- makeSimmap(phy, data[,c(1,2)], matrix(c(-1e-5,1e-5,1e-5,-1e-5), 2, 2), rate.cat = 1)[[1]]
map_y <- makeSimmap(phy, data[,c(1,3)], matrix(c(-1e-5,1e-5,1e-5,-1e-5), 2, 2), rate.cat = 1)[[1]]
p_x <- as.grob(~plotSimmap(map_x, color=cols2, fsize = 1e-10, lwd = lwd))
p_y <- as.grob(~plotSimmap(map_y, color=cols3, fsize = 1e-10, lwd = lwd, direction="leftwards"))
er <- rbind(data.frame(type = types[1], evidence_ratio = er_a), 
            data.frame(type = types[2], evidence_ratio = er_b))
er$type <- factor(er$type, levels = types)
p_out <- ggplot(er, aes(x = type, y = evidence_ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.8, alpha=0.5) +
  coord_cartesian(ylim=c(0, 50)) + 
  ggtitle("c) Modified Darwin's scenario (both)") +
  theme_bw() +
  ylab("Support for correlation (evidence ratio)") + 
  xlab("") +
  geom_hline(yintercept = 2.7, size=1, linetype="dashed") + 
  theme(legend.position = "none", plot.title = element_text(size = 21), axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 20)) +
  scale_x_discrete(labels = c('IM','HMIM'))
bth_plot <- grid.arrange(p_out, p_x, p_y, nrow=1, widths=c(1.2,1,1))

### final plot
evid_ratio_plot <- grid.arrange(out_plot, ins_plot, bth_plot, ncol=1)
ggsave(filename = "figures/in-progress/EvidenceRatios.pdf", plot = evid_ratio_plot, width = 12, height = 20, units = "in")

####### PART 7: introducing the simplified models (full results)
dAICTable <- rbind(rowMeans(apply(darwin_res, 1, function(x) x-min(x))),
                   rowMeans(apply(burst_res, 1, function(x) x-min(x))),
                   c(NA, rowMeans(apply(outside_res, 1, function(x) x-min(x)))),
                   c(NA, rowMeans(apply(inside_res, 1, function(x) x-min(x)))),
                   c(NA, rowMeans(apply(both_res, 1, function(x) x-min(x)))))

dAICTable_sd <- rbind(apply(apply(darwin_res, 1, function(x) x-min(x)), 1, sd),
                   apply(apply(burst_res, 1, function(x) x-min(x)), 1, sd),
                   c(NA, apply(apply(outside_res, 1, function(x) x-min(x)), 1, sd)),
                   c(NA, apply(apply(inside_res, 1, function(x) x-min(x)), 1, sd)),
                   c(NA, apply(apply(both_res, 1, function(x) x-min(x)), 1, sd)))

t(round(dAICTable, 1))
t(round(dAICTable_sd, 1))
rownames(dAICTable) <- c("Darwin's scenario", "Unreplicated bursts scenario", "Modified Darwin's scenario (outside)", "Modified Darwin's scenario (inside)", "Modified Darwin's scenario (both)")
dAICTable <- round(dAICTable, 2)
colnames(dAICTable) <- gsub("fit_", " ", colnames(dAICTable))
colnames(dAICTable) <- gsub("_", " ", colnames(dAICTable))
write.csv(t(dAICTable), file = "dAICtable.csv")

