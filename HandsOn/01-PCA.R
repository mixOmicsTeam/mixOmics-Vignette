
## ----eval=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------
library(mixOmics)
data(multidrug)
X <- multidrug$ABC.trans
dim(X) # Check dimensions of data


## ----screeplot-fig, fig.cap='(ref:screeplot-fig)'--------------------------------------------------------------------------------------------------
tune.pca.multi <- tune.pca(X, ncomp = 10, scale = TRUE)
plot(tune.pca.multi)
# tune.pca.multidrug$cum.var       # Outputs cumulative proportion of variance


## ----echo=TRUE, message=FALSE----------------------------------------------------------------------------------------------------------------------
final.pca.multi <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)
# final.pca.multi  # Lists possible outputs


## --------------------------------------------------------------------------------------------------------------------------------------------------
final.pca.multi$var.tot


## --------------------------------------------------------------------------------------------------------------------------------------------------
final.pca.multi$prop_expl_var$X


## ---- eval = TRUE----------------------------------------------------------------------------------------------------------------------------------
final.pca.multi$cum.var


## ----echo=TRUE, message=FALSE----------------------------------------------------------------------------------------------------------------------
# Top variables on the first component only:
head(selectVar(final.pca.multi, comp = 1)$value)


## ----pca-ABCtrans-group, fig.cap='(ref:pca-ABCtrans-group)'----------------------------------------------------------------------------------------
plotIndiv(final.pca.multi,
          comp = c(1, 2),   # Specify components to plot
          ind.names = TRUE, # Show row names of samples
          group = multidrug$cell.line$Class,
          title = 'ABC transporters, PCA comp 1 - 2',
          legend = TRUE, legend.title = 'Cell line')


## ---- eval = FALSE---------------------------------------------------------------------------------------------------------------------------------
## # Interactive 3D plot will load the rgl library.
## plotIndiv(final.pca.multi, style = '3d',
##            group = multidrug$cell.line$Class,
##           title = 'ABC transporters, PCA comp 1 - 3')


## ---- eval = FALSE---------------------------------------------------------------------------------------------------------------------------------
## plotVar(final.pca.multi, comp = c(1, 2),
##         var.names = TRUE,
##         cex = 3,         # To change the font size
##         # cutoff = 0.5,  # For further cutoff
##         title = 'Multidrug transporter, PCA comp 1 - 2')


## ----pca-var-ABStrans-ccp, eval=TRUE, echo = FALSE, fig.cap='(ref:pca-var-ABStrans-ccp)'-----------------------------------------------------------
col.var <- c(rep(color.mixo(1), ncol(X)))
names(col.var) = colnames(X)
col.var[c('ABCE1', 'ABCB8')] = color.mixo(2)
col.var[c('ABCA8')] = color.mixo(4)
col.var[c('ABCC2')] = color.mixo(5)
col.var[c('ABCC12', 'ABCD2')] = color.mixo(6)

plotVar(final.pca.multi, comp = c(1, 2),
        var.names = TRUE,
        col = list(col.var),
        cex = 3,
        title = 'Multidrug transporter, PCA comp 1 - 2')


## ----pca-ABStrans-biplot, eval=TRUE, fig.cap='(ref:pca-ABStrans-biplot)'---------------------------------------------------------------------------
biplot(final.pca.multi, group = multidrug$cell.line$Class, 
       legend.title = 'Cell line')


## ----pca-ABCtrans-boxplot, eval=TRUE, fig.cap='(ref:pca-ABCtrans-boxplot)'-------------------------------------------------------------------------
ABCC2.scale <- scale(X[, 'ABCC2'], center = TRUE, scale = TRUE)

boxplot(ABCC2.scale ~
        multidrug$cell.line$Class, col = color.mixo(1:9),
        xlab = 'Cell lines', ylab = 'Expression levels, scaled',
        par(cex.axis = 0.5), # Font size
        main = 'ABCC2 transporter')


## --------------------------------------------------------------------------------------------------------------------------------------------------
grid.keepX <- c(seq(5, 30, 5))
# grid.keepX  # To see the grid

set.seed(30) # For reproducibility with this handbook, remove otherwise
tune.spca.result <- tune.spca(X, ncomp = 3, 
                              folds = 5, 
                              test.keepX = grid.keepX, nrepeat = 10) 

# Consider adding up to 50 repeats for more stable results
tune.spca.result$choice.keepX


## ----spca-tuning-ABStrans, fig.cap='(ref:spca-tuning-ABStrans)'------------------------------------------------------------------------------------
plot(tune.spca.result)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# By default center = TRUE, scale = TRUE
keepX.select <- tune.spca.result$choice.keepX[1:2]

final.spca.multi <- spca(X, ncomp = 2, keepX = keepX.select)

# Proportion of explained variance:
final.spca.multi$prop_expl_var$X


## ----spca-ABCtrans-group, fig.cap='(ref:spca-ABCtrans-group)'--------------------------------------------------------------------------------------
plotIndiv(final.spca.multi,
          comp = c(1, 2),   # Specify components to plot
          ind.names = TRUE, # Show row names of samples
          group = multidrug$cell.line$Class,
          title = 'ABC transporters, sPCA comp 1 - 2',
          legend = TRUE, legend.title = 'Cell line')


## ----spca-ABCtrans-biplot, eval=TRUE, fig.cap='(ref:spca-ABCtrans-biplot)'-------------------------------------------------------------------------
biplot(final.spca.multi, group = multidrug$cell.line$Class, 
       legend =FALSE)


## ---- eval = FALSE---------------------------------------------------------------------------------------------------------------------------------
## plotVar(final.spca.multi, comp = c(1, 2), var.names = TRUE,
##         cex = 3, # To change the font size
##         title = 'Multidrug transporter, sPCA comp 1 - 2')


## ----spca-var-ABStrans-ccp, eval=TRUE, fig.cap='(ref:spca-var-ABCtrans-ccp)', echo = FALSE---------------------------------------------------------
col.var <- c(rep(color.mixo(1), ncol(X)))
names(col.var) <- colnames(X)
col.var[c("ABCA9", "ABCB5", "ABCC2","ABCD1")] <- color.mixo(4)

plotVar(final.spca.multi, comp = c(1, 2), var.names = TRUE,
        col = list(col.var), cex = 3, # To change the font size
        title = 'Multidrug transporter, sPCA comp 1 - 2')


## ----echo=TRUE, message=FALSE----------------------------------------------------------------------------------------------------------------------
# On the first component, just a head
head(selectVar(final.spca.multi, comp = 2)$value)


## ----spca-plotLoading, eval=TRUE, fig.cap='(ref:spca-plotLoading)'---------------------------------------------------------------------------------
plotLoadings(final.spca.multi, comp = 2)

