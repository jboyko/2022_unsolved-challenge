
library(corHMM)
library(interp)
library(lhs)
library(data.table)
library(viridis)
library(phytools)
library(parallel)

######################################################################################################################################
######################################################################################################################################
### FUNCTIONS FOR GETTING CONTOUR COMPONENTS
######################################################################################################################################
######################################################################################################################################

contourSearchPoints <- function(variables, lower, upper, nreps){
    #Creates a latin square design of the parameter space:
    X <- randomLHS(nreps, variables)
    param.points <- matrix(0, nrow=nreps, ncol=variables)
    param.points[,1] <- qunif(X[,1], lower[1], upper[1])
    param.points[,2] <- qunif(X[,2], lower[2], upper[2])
    return(param.points)
}


#Choose 2 to fix, estimate rest
CorHMMSemiFixed <- function(x, phy, model.set.final, fixed.pars, fixed.pair.index, order.test, lewis.asc.bias){
    
    for(pair.index in 1:2){
        model.set.final$index.matrix[fixed.pair.index[pair.index,1], fixed.pair.index[pair.index,2]] <- length(x) + pair.index
    }
    model.set.final$np <- model.set.final$np + 2
    rate.mat <- model.set.final$index.matrix
    rate.mat[is.na(rate.mat)] <- max(rate.mat,na.rm=TRUE)+1
    model.set.final$rate <- rate.mat
    pars <- c(x, log(fixed.pars))
    
    loglik <- corHMM:::dev.corhmm(pars, phy=phy, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p="yang", rate.cat=1, order.test=order.test, lewis.asc.bias=FALSE)
    return(loglik)
}


#Function that calls the semifixed code
contourSearchCorHMM <- function(phy, data, init.val, param.points, collapse=FALSE, fixed.pair.index, n.cores=NULL) {
    
    if(collapse == FALSE){
        rate.mat <- getFullMat(list(getRateCatMat(2),getRateCatMat(2)), getRateCatMat(2))
        colnames(rate.mat) <- rownames(rate.mat) <- c("0,0", "0,1", "1,0", "1,1")
        rate.mat[rate.mat != 0] <- 1:8
    }else{
        rate.mat <- NULL
    }

    input.data <- data
    nCol <- dim(data)[2]
    CorData <- corHMM:::corProcessData(data, collapse = collapse)
    data <- CorData$corData
    if(length(grep("&", CorData$corData[,2])) > 0){
        non_and_chars <- as.numeric(CorData$corData[,2][-grep("&", CorData$corData[,2])])
        and_chars <- as.numeric(unlist(strsplit(CorData$corData[,2][grep("&", CorData$corData[,2])], "&")))
        nObs <- max(c(non_and_chars, and_chars))
    }else{
        nObs <- max(as.numeric(CorData$corData[,2]))
    }
    order.test <- TRUE
    matching <- corHMM:::match.tree.data(phy,data)
    data <- matching$data
    phy <- matching$phy
    
    fixed.pairs.pars.index <- c()
    for(par.index in 1:2){
        fixed.pairs.pars.index <- c(fixed.pairs.pars.index, rate.mat[fixed.pair.index[par.index,1], fixed.pair.index[par.index,2]])
    }
    rate.mat.red <- dropStateMatPars(rate.mat, fixed.pairs.pars.index)
    
    model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy, data=input.data, rate.cat=1, ntraits=nObs, model=NULL, rate.mat=rate.mat.red, collapse=collapse)
    phy <- reorder(phy, "pruningwise")
    phy$node.label <- NULL

    nreps <- dim(param.points)[1]
    res <- matrix(,nreps,3)
    init.vals <- rep(init.val, max(rate.mat.red, na.rm=T))
    lower <- rep(log(1e-9), length(init.vals))
    upper <- rep(log(100), length(init.vals))
    
    if(is.null(n.cores)){
        n.cores=1
    }
    PointEval <- function(nrep.index){
        fixed.pars <- c(param.points[nrep.index,1], param.points[nrep.index,2])
        opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
        out <- nloptr(x0=log(init.vals), eval_f=CorHMMSemiFixed, lb=lower, ub=upper, opts=opts, phy=phy, model.set.final=model.set.final, fixed.pars=fixed.pars, fixed.pair.index=fixed.pair.index, order.test=order.test, lewis.asc.bias=FALSE)
        save(out, file="check.Rsave")
        return(c(-out$objective, fixed.pars[1], fixed.pars[2]))
    }
    res.list <- mclapply(1:nreps, PointEval, mc.cores=n.cores)
    res <- matrix(unlist(res.list), ncol = 3, byrow = TRUE)
    return(res)
}


#Function getting contour data
GetContour <- function(phy, data, init.val, fixed.pair.index = rbind(c("0,0", "1,0"), c("1,0", "0,0")), nreps=10, n.cores=NULL){

    if(any(phy$edge.length<=1e-5)){
        phy$edge.length[phy$edge.length<=1e-5] <- 1e-5
    }
    param.points <- contourSearchPoints(variables=2, lower=c(0,0), upper=c(10,10), nreps=nreps)
    surface.data <- contourSearchCorHMM(phy=phy, data=data, init.val=init.val, param.points=param.points, fixed.pair.index=fixed.pair.index, n.cores=n.cores)
    obj <- list(surface.data=surface.data, focal.params=paste(fixed.pair.index[,1], "->", fixed.pair.index[,2]))
    
    return(obj)
}



#Function that plots contour surface
plot.contour <- function(x, mle.point=NULL, mle.point.col="blue", levels=c(0:20*0.1), xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, col=grey.colors(21, start=0, end=1), ...){
    
    if(is.null(xlab)){
        xlab = x$focal.params[1]
    }
    if(is.null(ylab)){
        ylab = x$focal.params[2]
    }
    if(is.null(xlim)){
        xlim = c(x$focal.params.lower[1], x$focal.params.upper[1], 1)
    }
    if(is.null(ylim)){
        ylim = c(x$focal.params.lower[2], x$focal.params.upper[2], 1)
    }
    
    mydata <- data.frame(x=matrix(x$surface.data[,2],ncol=1),y=matrix(x$surface.data[,3],ncol=1),z=matrix(x$surface.data[,1],ncol=1))
    mydata$z <- (-1)*mydata$z
    mydata$z <- mydata$z-min(mydata$z)
    
    interp.res <- interp(x=mydata$x, y=mydata$y, z=mydata$z, xo=seq(min(mydata$x), max(mydata$x), length = 400), yo=seq(min(mydata$y), max(mydata$y),length = 400), duplicate=FALSE)
    plot(NA,xlab="", ylab="", frame=FALSE, axes=FALSE, xaxs="i", yaxs="i", ylim=ylim[1:2], xlim=xlim[1:2], ...)
    .filled.contour(interp.res$x, interp.res$y, interp.res$z, levels=levels, col=col)
    par(tck=.01)
    axis(2, at = seq(ylim[1], ylim[2], by = ylim[3]), las=1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
    axis(1, at = seq(xlim[1], xlim[2], by = xlim[3]), las=1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
    if(!is.null(mle.point)){
        points(x=mle.point[1], y=mle.point[2], pch=19, col=mle.point.col)
    }
    title(xlab=xlab, line=2.5)
    title(ylab=ylab, line=2)
}


######################################################################################################################################
######################################################################################################################################
### Generate surFace Data -- UNCOMMENT AND IT SHOULD WORK
######################################################################################################################################
######################################################################################################################################

### Getting everything setup ###
load("maddfitz_hmm_res.Rsave")
cor.fit <- maddfitz_hmm_res[[1]]$darwin_res$cor
phy <- cor.fit$phy
data <- cor.fit$data
rate.mat <- getFullMat(list(getRateCatMat(2),getRateCatMat(2)), getRateCatMat(2))
colnames(rate.mat) <- rownames(rate.mat) <- c("0,0", "1,0", "0,1", "1,1")
rate.mat[rate.mat != 0] <- 1
init.val <- corHMM(phy, data, rate.cat=1, rate.mat=rate.mat, collapse=FALSE)
#########################################

##TEST##
#load("check.Rsave")
#rate.mat <- getFullMat(list(getRateCatMat(2),getRateCatMat(2)), getRateCatMat(2))
#colnames(rate.mat) <- rownames(rate.mat) <- c("0,0", "0,1", "1,0", "1,1")
#rate.mat[rate.mat != 0] <- 1:8
#rate.mat.red <- dropStateMatPars(rate.mat, c(3,5))
#rate.mat.red[1,2] <- 7
#rate.mat.red[1,3] <- 8
#init.val <- corHMM(phy, data, rate.cat=1, rate.mat=rate.mat.red, collapse=FALSE, p=c(exp(out$solution), c(67.49949, 42.61204))

############## 00->01 and 00->10 ##############
cont.info1 <- GetContour(phy, data, init.val=init.val$solution[2,1], fixed.pair.index = rbind(c("0,0", "0,1"), c("0,0", "1,0")), nreps=5000, n.cores=50)
cont.info1$surface.data <- cont.info1$surface.data[is.finite(cont.info1$surface.data[,1]),]
save(cont.info1, file="test1_5000points0001_0010.Rsave")
#########################################

############## 11->01 and 11->10 ##############
cont.info2 <- GetContour(phy, data, init.val=init.val$solution[2,1], fixed.pair.index = rbind(c("1,1", "0,1"), c("1,1", "1,0")), nreps=5000, n.cores=50)
cont.info2$surface.data <- cont.info2$surface.data[is.finite(cont.info2$surface.data[,1]),]
save(cont.info2, file="test1_5000points1101_1110.Rsave")
#########################################


############## 01->00 and 01->11 ##############
#cont.info3 <- GetContour(phy, data, init.val=init.val$solution[2,1], fixed.pair.index = rbind(c("0,1", "0,0"), c("0,1", "1,1")), nreps=5000, n.cores=50)
#cont.info3$surface.data <- cont.info3$surface.data[is.finite(cont.info3$surface.data[,1]),]
#save(cont.info3, file="test1_5000points0100_0111.Rsave")
#########################################


############## 10->00 and 10->11 ##############
#cont.info4 <- GetContour(phy, data, init.val=init.val$solution[2,1], fixed.pair.index = rbind(c("1,0", "0,0"), c("1,0", "1,1")), nreps=5000, n.cores=50)
#cont.info4$surface.data <- cont.info4$surface.data[is.finite(cont.info4$surface.data[,1]),]
#save(cont.info4, file="test1_5000points1000_1011.Rsave")
#########################################



######################################################################################################################################
######################################################################################################################################
### Generate Contour 2 Plot
######################################################################################################################################
######################################################################################################################################

#userID  <-  Sys.info()['user']

#set dirs based on userID
#switch(userID,
#"jeremybeaulieu" = { # assume run from dir selon_work/
#    contour.dir <- "/Users/jeremybeaulieu/hisse_fossil_sims/Single_Rate/Contour/Surfaces/";
#    quick.sim.dir <- "/Users/jeremybeaulieu/hisse_fossil_sims/Single_Rate/Contour/QuickStratSim/";
#    out.dir <- "/Users/jeremybeaulieu/hisse_fossil_sims/Tables_Figures/";},
#)

#pdf(paste0(out.dir, "Figure2.pdf"), width=12, height=8)
#par(mfcol=c(2,3),mar=c(4,4.5,0.5,0.5), oma=c(1.5,2,1,1))

#load(paste(contour.dir, "test1_5000points010.Rsave", sep=""))
#cont.info1$focal.params <- c(expression(lambda+mu), expression(mu/lambda))

#plot.contour(cont.info1, mle.point=c(0.4+0.3, 0.3/0.4), mle.point.col="white", levels=c(0, .5, 1, 1.5, 2, 2.5), ylim=c(0,1,.1), xlim=c(0,1.6,.2), col=magma(5))
#mtext("a)",side=3, line=0, adj=0, cex=1)
#abline(v=.7, lty=2)
#abline(h=.75, lty=2)
#text(.1, .95, expression(Extant~only), pos=4)

#plot(0,type='n',axes=FALSE,ann=FALSE)
#load(paste(contour.dir, "simTreeContour100.Rsave", sep=""))

#dev.off()

pdf("Figure1_draft.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4.5,0.5,0.5), oma=c(1.5,2,1,1))

load("test1_5000points0001_0010.Rsave")
cont.info1$focal.params <- c(expression(q["00->01"]), expression(q["00->10"]))
plot.contour(cont.info1, mle.point=c(0.000000001, 0.012090005), mle.point.col="white", levels=c(0, .5, 1, 1.5, 2, 2.5), ylim=c(0,10,2), xlim=c(0,10,2), col=viridis(5))
mtext("a)",side=3, line=0, adj=0, cex=1)
abline(v=0.000000001, lty=2)
abline(h=0.012090005, lty=2)
mtext("a)",side=3, line=0, adj=0, cex=1)

load("test1_5000points0100_0111.Rsave")
cont.info3$focal.params <- c(expression(q["01->00"]), expression(q["01->11"]))
plot.contour(cont.info3, mle.point=c(1.000000e+02, 1e-09), mle.point.col="white", levels=c(0, .5, 1, 1.5, 2, 2.5), ylim=c(0,100,20), xlim=c(0,100,20), col=viridis(5))
abline(v=1.000000e+02, lty=2)
abline(h=1e-09, lty=2)
mtext("c)",side=3, line=0, adj=0, cex=1)

load("test1_5000points1101_1110.Rsave")
cont.info2$focal.params <- c(expression(q["11->01"]), expression(q["11->10"]))
plot.contour(cont.info2, mle.point=c(0.014154138, 0.000000001), mle.point.col="white", levels=c(0, .5, 1, 1.5, 2, 2.5), ylim=c(0,10,2), xlim=c(0,10,2), col=viridis(5))
abline(v=0.014154138, lty=2)
abline(h=0.000000001, lty=2)
mtext("b)",side=3, line=0, adj=0, cex=1)

load("test1_5000points1000_1011.Rsave")
cont.info4$focal.params <- c(expression(q["10->00"]), expression(q["10->11"]))
plot.contour(cont.info4, mle.point=c(0.014154138, 0.000000001), mle.point.col="white", levels=c(0, .5, 1, 1.5, 2, 2.5), ylim=c(0,100,20), xlim=c(0,100,20), col=viridis(5))
abline(v=1.136044e-09, lty=2)
abline(h=1e+02, lty=2)
mtext("d)",side=3, line=0, adj=0, cex=1)

dev.off()



