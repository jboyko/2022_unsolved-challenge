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

summarizeListElement <- function(maddfitz_simmed_res_element){
  out_table <- rbind(unlist(lapply(maddfitz_simmed_res_element$ind_dat_fits, "[[", "AIC")),
                     unlist(lapply(maddfitz_simmed_res_element$cor_dat_fits, "[[", "AIC")),
                     unlist(lapply(maddfitz_simmed_res_element$ind_2_dat_fits, "[[", "AIC")))
  rownames(out_table) <- c("ind_dat", "cor_dat", "ind_2_dat")
  return(out_table)
}

getAICweight <- function(AICcs){
  dAIC <- AICcs - min(AICcs)
  AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
  return(AICwt)
}


# notes
# bind the tree for multiple replicates of darwin's scenario
# strong statistical justification - totally dispraportionate evidence comparing simple independent to correlated
load("maddfitz_simmed_res.Rsave")

table_list <- lapply(maddfitz_simmed_res, summarizeListElement)
plot_data <- as.data.frame(do.call(rbind, lapply(table_list, function(x) melt(apply(x, 1, getAICweight)))))

colnames(plot_data)

ggplot(plot_data, aes(x = Var1, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.8, alpha=0.5, width = 0.3) +
  coord_cartesian(ylim=c(0, 1)) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 15)) + 
  ylab("AICwt") + 
  xlab("") +
  theme(legend.position = "none") +
  facet_wrap(~Var2)


