
## ---- message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------------
library(mixOmics)
data(stemcells)

# The combined data set X
X <- stemcells$gene
dim(X)

# The outcome vector Y:  
Y <- stemcells$celltype 
length(Y) 

summary(Y)


## --------------------------------------------------------------------------------------------------------------------------------------------------
study <- stemcells$study

# Number of samples per study:
summary(study)

# Experimental design
table(Y,study)


## ----MINT-perf,fig.cap='(ref:MINT-perf)'-----------------------------------------------------------------------------------------------------------
mint.plsda.stem <- mint.plsda(X = X, Y = Y, study = study, ncomp = 5)

set.seed(2543) # For reproducible results here, remove for your own analyses
perf.mint.plsda.stem <- perf(mint.plsda.stem) 

plot(perf.mint.plsda.stem)


## ---- results = 'hold'-----------------------------------------------------------------------------------------------------------------------------
perf.mint.plsda.stem$global.error$BER
# Type also:
# perf.mint.plsda.stem$global.error


## ----MINT-plsda-indiv, fig.cap='(ref:MINT-plsda-indiv)'--------------------------------------------------------------------------------------------
final.mint.plsda.stem <- mint.plsda(X = X, Y = Y, study = study, ncomp = 2)

#final.mint.plsda.stem # Lists the different functions

plotIndiv(final.mint.plsda.stem, legend = TRUE, title = 'MINT PLS-DA', 
          subtitle = 'stem cell study', ellipse = T)


## ----stem-plsda-indiv, fig.cap='(ref:stem-plsda-indiv)'--------------------------------------------------------------------------------------------
plsda.stem <- plsda(X = X, Y = Y, ncomp = 2)

plotIndiv(plsda.stem, pch = study,
          legend = TRUE, title = 'Classic PLS-DA',
          legend.title = 'Cell type', legend.title.pch = 'Study')


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(2543)  # For a reproducible result here, remove for your own analyses
tune.mint.splsda.stem <- tune(X = X, Y = Y, study = study, 
                 ncomp = 2, test.keepX = seq(1, 100, 1),
                 method = 'mint.splsda', #Specify the method
                 measure = 'BER',
                 dist = "centroids.dist")

#tune.mint.splsda.stem # Lists the different types of outputs

# Mean error rate per component and per tested keepX value:
#tune.mint.splsda.stem$error.rate[1:5,]


## --------------------------------------------------------------------------------------------------------------------------------------------------
tune.mint.splsda.stem$choice.keepX


## ----MINT-tune, fig.cap='(ref:MINT-tune)'----------------------------------------------------------------------------------------------------------
plot(tune.mint.splsda.stem)


## --------------------------------------------------------------------------------------------------------------------------------------------------
final.mint.splsda.stem <- mint.splsda(X = X, Y = Y, study = study, ncomp = 2,  
                              keepX = tune.mint.splsda.stem$choice.keepX)

#mint.splsda.stem.final # Lists useful functions that can be used with a MINT object


## ----MINT-indiv-global, echo = TRUE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/MINT/'-----------------------------------------------
plotIndiv(final.mint.splsda.stem, study = 'global', legend = TRUE, 
          title = 'Stem cells, MINT sPLS-DA', 
          subtitle = 'Global', ellipse = T)


## ----MINT-indiv-local, echo = TRUE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/MINT/'------------------------------------------------
plotIndiv(final.mint.splsda.stem, study = 'all.partial', legend = TRUE, 
          title = 'Stem cells, MINT sPLS-DA', 
          subtitle = paste("Study",1:4))


## ----MINT-indiv, fig.cap='(ref:MINT-indiv)', echo=FALSE, eval=TRUE, fig.align='center', out.width = '50%', fig.height=4, fig.ncol = 2, fig.subcap=c('', '')----
knitr::include_graphics(c('Figures/MINT/MINT-indiv-global-1.png', 'Figures/MINT/MINT-indiv-local-1.png'))


## ---- echo = TRUE----------------------------------------------------------------------------------------------------------------------------------
plotVar(final.mint.splsda.stem)


## ----MINT-var-col, fig.cap='(ref:MINT-var-col)', echo = FALSE--------------------------------------------------------------------------------------
output <- plotVar(final.mint.splsda.stem, cex = 2, plot = FALSE)
col.var <- rep(color.mixo(1), ncol(X))
names(col.var) = colnames(X)
col.var[rownames(output)[output$x > 0.8]] <- color.mixo(4)
col.var[rownames(output)[output$x < (-0.8)]] <- color.mixo(5)

plotVar(final.mint.splsda.stem, cex = 2, col= list(col.var))


## ----MINT-cim, fig.cap='(ref:MINT-cim)', fig.width=10, fig.height=8--------------------------------------------------------------------------------
# If facing margin issues, use either X11() or save the plot using the
# arguments save and name.save
cim(final.mint.splsda.stem, comp = 1, margins=c(10,5), 
    row.sideColors = color.mixo(as.numeric(Y)), row.names = FALSE,
    title = "MINT sPLS-DA, component 1")


## ----MINT-network, fig.cap='(ref:MINT-network)',fig.width=10, fig.height=8-------------------------------------------------------------------------
# If facing margin issues, use either X11() or save the plot using the
# arguments save and name.save
network(final.mint.splsda.stem, comp = 1,
        color.node = c(color.mixo(1), color.mixo(2)), 
        shape.node = c("rectangle", "circle"))


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Just a head
head(selectVar(final.mint.plsda.stem, comp = 1)$value)


## ----MINT-loading, fig.cap='(ref:MINT-loading)'----------------------------------------------------------------------------------------------------
plotLoadings(final.mint.splsda.stem, contrib = "max", method = 'mean', comp=1, 
             study="all.partial", title="Contribution on comp 1", 
             subtitle = paste("Study",1:4))


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)  # For reproducible results here, remove for your own study
perf.mint.splsda.stem.final <- perf(final.mint.plsda.stem, dist = 'centroids.dist')

perf.mint.splsda.stem.final$global.error


## ----MINT-auc1, echo = TRUE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/MINT/'-------------------------------------------------------
auroc(final.mint.splsda.stem, roc.comp = 1)


## ----MINT-auc2, echo = TRUE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/MINT/'-------------------------------------------------------
auroc(final.mint.splsda.stem, roc.comp = 1, roc.study = '2')


## ----MINT-auc, fig.cap='(ref:MINT-auc)', echo=FALSE, eval=TRUE, fig.align='center', out.width = '50%', fig.height=4, fig.ncol = 2, fig.subcap=c('', '')----
knitr::include_graphics(c('Figures/MINT/MINT-auc1-1.png', 'Figures/MINT/MINT-auc2-1.png'))


## --------------------------------------------------------------------------------------------------------------------------------------------------
# We predict on study 3
indiv.test <- which(study == "3")

# We train on the remaining studies, with pre-tuned parameters
mint.splsda.stem2 <- mint.splsda(X = X[-c(indiv.test), ], 
                                Y = Y[-c(indiv.test)], 
                                study = droplevels(study[-c(indiv.test)]), 
                                ncomp = 1,  
                                keepX = 30)

mint.predict.stem <- predict(mint.splsda.stem2, newdata = X[indiv.test, ], 
                        dist = "centroids.dist",
                        study.test = factor(study[indiv.test]))

# Store class prediction with a model with 1 comp
indiv.prediction <- mint.predict.stem$class$centroids.dist[, 1]

# The confusion matrix compares the real subtypes with the predicted subtypes
conf.mat <- get.confusion_matrix(truth = Y[indiv.test],
                     predicted = indiv.prediction)

conf.mat


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Prediction error rate
(sum(conf.mat) - sum(diag(conf.mat)))/sum(conf.mat)

