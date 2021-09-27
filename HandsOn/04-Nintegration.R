
## ----message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------------------
library(mixOmics)
data(breast.TCGA)

# Extract training data and name each data frame
# Store as list
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna, 
          protein = breast.TCGA$data.train$protein)

# Outcome
Y <- breast.TCGA$data.train$subtype
summary(Y)


## --------------------------------------------------------------------------------------------------------------------------------------------------
design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) <- 0
design 


## ---- results='hold'-------------------------------------------------------------------------------------------------------------------------------
res1.pls.tcga <- pls(X$mRNA, X$protein, ncomp = 1)
cor(res1.pls.tcga$variates$X, res1.pls.tcga$variates$Y)

res2.pls.tcga <- pls(X$mRNA, X$miRNA, ncomp = 1)
cor(res2.pls.tcga$variates$X, res2.pls.tcga$variates$Y)

res3.pls.tcga <- pls(X$protein, X$miRNA, ncomp = 1)
cor(res3.pls.tcga$variates$X, res3.pls.tcga$variates$Y)


## ----diablo-perf, message=FALSE, fig.cap='(ref:diablo-perf)'---------------------------------------------------------------------------------------
diablo.tcga <- block.plsda(X, Y, ncomp = 5, design = design)

set.seed(123) # For reproducibility, remove for your analyses
perf.diablo.tcga = perf(diablo.tcga, validation = 'Mfold', folds = 10, nrepeat = 10)

#perf.diablo.tcga$error.rate  # Lists the different types of error rates

# Plot of the error rates based on weighted vote
plot(perf.diablo.tcga)


## --------------------------------------------------------------------------------------------------------------------------------------------------
perf.diablo.tcga$choice.ncomp$WeightedVote


## --------------------------------------------------------------------------------------------------------------------------------------------------
ncomp <- perf.diablo.tcga$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]


## ---- echo = TRUE, eval = FALSE--------------------------------------------------------------------------------------------------------------------
## # This code may take several min to run, parallelisation is possible
## set.seed(123) # For reproducibility with this handbook, remove otherwise
## test.keepX <- list(mRNA = c(5:9, seq(10, 25, 5)),
##                    miRNA = c(5:9, seq(10, 20, 2)),
##                    proteomics = c(seq(5, 25, 5)))
## 
## tune.diablo.tcga <- tune.block.splsda(X, Y, ncomp = 2,
##                               test.keepX = test.keepX, design = design,
##                               validation = 'Mfold', folds = 10, nrepeat = 1,
##                               dist = "centroids.dist")


## ---- echo = FALSE, eval = TRUE, warning=FALSE-----------------------------------------------------------------------------------------------------
# chunk takes about 2 min to run
set.seed(123) # for reproducibility
test.keepX <- list(mRNA = c(5:9, seq(10, 25, 5)),
                   miRNA = c(5:9, seq(10, 20, 2)),
                   proteomics = c(seq(5, 25, 5)))

tune.diablo.tcga <- tune.block.splsda(X, Y, ncomp = 2, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1, 
                              # use two CPUs for faster computation
                              BPPARAM = BiocParallel::SnowParam(workers = 2),
                              dist = "centroids.dist")


## --------------------------------------------------------------------------------------------------------------------------------------------------
list.keepX <- tune.diablo.tcga$choice.keepX
list.keepX


## ---- eval = FALSE---------------------------------------------------------------------------------------------------------------------------------
## list.keepX <- list( mRNA = c(8, 25), miRNA = c(14,5), protein = c(10, 5))


## ---- message = TRUE-------------------------------------------------------------------------------------------------------------------------------
diablo.tcga <- block.splsda(X, Y, ncomp = ncomp, 
                            keepX = list.keepX, design = design)
#diablo.tcga   # Lists the different functions of interest related to that object


## --------------------------------------------------------------------------------------------------------------------------------------------------
diablo.tcga$design


## ---- eval = FALSE---------------------------------------------------------------------------------------------------------------------------------
## # mRNA variables selected on component 1
## selectVar(diablo.tcga, block = 'mRNA', comp = 1)


## ----plot-diablo, message=FALSE, fig.cap='(ref:plot-diablo)'---------------------------------------------------------------------------------------
plotDiablo(diablo.tcga, ncomp = 1)


## ----diablo-plotindiv, message=FALSE, fig.cap='(ref:diablo-plotindiv)'-----------------------------------------------------------------------------
plotIndiv(diablo.tcga, ind.names = FALSE, legend = TRUE, 
          title = 'TCGA, DIABLO comp 1 - 2')


## ----diablo-plotarrow, message=FALSE, fig.cap='(ref:diablo-plotarrow)'-----------------------------------------------------------------------------
plotArrow(diablo.tcga, ind.names = FALSE, legend = TRUE, 
          title = 'TCGA, DIABLO comp 1 - 2')


## ----diablo-plotvar, message=FALSE, fig.cap='(ref:diablo-plotvar)'---------------------------------------------------------------------------------
plotVar(diablo.tcga, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen'),
        title = 'TCGA, DIABLO comp 1 - 2')


## ----diablo-circos, message=FALSE, fig.cap='(ref:diablo-circos)'-----------------------------------------------------------------------------------
circosPlot(diablo.tcga, cutoff = 0.7, line = TRUE, 
           color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)


## ---- eval = FALSE, message=FALSE------------------------------------------------------------------------------------------------------------------
## # X11()   # Opens a new window
## network(diablo.tcga, blocks = c(1,2,3),
##         cutoff = 0.4,
##         color.node = c('darkorchid', 'brown1', 'lightgreen'),
##         # To save the plot, comment out otherwise
##         save = 'png', name.save = 'diablo-network'
##         )


## ----diablo-network2, fig.cap ='(ref:diablo-network)', echo = FALSE, eval=TRUE---------------------------------------------------------------------
knitr::include_graphics("Figures/DIABLO/diablo-network.png")


## ----eval = FALSE----------------------------------------------------------------------------------------------------------------------------------
## # Not run
## library(igraph)
## myNetwork <- network(diablo.tcga, blocks = c(1,2,3), cutoff = 0.4)
## write.graph(myNetwork$gR, file = "myNetwork.gml", format = "gml")


## ----diablo-loading, message=FALSE, fig.cap='(ref:diablo-loading)'---------------------------------------------------------------------------------
plotLoadings(diablo.tcga, comp = 1, contrib = 'max', method = 'median')


## ----diablo-cim, message=FALSE, out.width = '70%', fig.cap='(ref:diablo-cim)'----------------------------------------------------------------------
cimDiablo(diablo.tcga, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = 1, margin=c(8,20), legend.position = "right")


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123) # For reproducibility with this handbook, remove otherwise
perf.diablo.tcga <- perf(diablo.tcga,  validation = 'Mfold', folds = 10, 
                         nrepeat = 10, dist = 'centroids.dist')

#perf.diablo.tcga  # Lists the different outputs


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Performance with Majority vote
perf.diablo.tcga$MajorityVote.error.rate


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Performance with Weighted vote
perf.diablo.tcga$WeightedVote.error.rate


## ----diablo-auroc, message=FALSE, fig.cap='(ref:diablo-auroc)'-------------------------------------------------------------------------------------
auc.diablo.tcga <- auroc(diablo.tcga, roc.block = "miRNA", roc.comp = 2,
                   print = FALSE)


## ---- message = FALSE------------------------------------------------------------------------------------------------------------------------------
# Prepare test set data: here one block (proteins) is missing
data.test.tcga <- list(mRNA = breast.TCGA$data.test$mrna, 
                      miRNA = breast.TCGA$data.test$mirna)

predict.diablo.tcga <- predict(diablo.tcga, newdata = data.test.tcga)
# The warning message will inform us that one block is missing

#predict.diablo # List the different outputs


## --------------------------------------------------------------------------------------------------------------------------------------------------
confusion.mat.tcga <- get.confusion_matrix(truth = breast.TCGA$data.test$subtype, 
                     predicted = predict.diablo.tcga$WeightedVote$centroids.dist[,2])
confusion.mat.tcga


## --------------------------------------------------------------------------------------------------------------------------------------------------
get.BER(confusion.mat.tcga)

