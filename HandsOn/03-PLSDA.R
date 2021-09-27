## ----global_options, include=FALSE----------------------------------------------------------------------------------------------------------------
library(knitr)
# global options
knitr::opts_chunk$set(dpi = 100, echo=TRUE, warning=FALSE, message=FALSE, eval = TRUE,
                      fig.show=TRUE, fig.width= 7,fig.height= 6,fig.align='center', out.width = '50%', message = FALSE,
                      fig.path= 'Figures/PLSDA/')

# The libraries to load
library(kableExtra)



## ---- include=FALSE-------------------------------------------------------------------------------------------------------------------------------
colorize <- function(color, x) {
  if (knitr::is_html_output()) {
    htmlcolor = "black"
    if(color == "blue"){
      htmlcolor = "#388ECC"
    }
    if(color == "orange"){
      htmlcolor = "#F68B33"
    }
    if(color == "grey"){
      htmlcolor = "#585858"
    }
    if(color == "green"){
      htmlcolor = "#009E73"
    }
    if(color == "pink"){
      htmlcolor = "#CC79A7"
    }
    if(color == "yellow"){
      htmlcolor = "#999900"
    }
    if(color == "darkred"){
      htmlcolor = "#CC0000"
    }
    sprintf("<span style='color: %s;'>%s</span>", htmlcolor, x)
  } else {
    sprintf("\\textcolor{%s}{%s}", color, x)
    }
}


## ----results = 'hold', message=FALSE--------------------------------------------------------------------------------------------------------------
library(mixOmics)
data(srbct)
X <- srbct$gene

# Outcome y that will be internally coded as dummy:
Y <- srbct$class 
dim(X); length(Y)


## -------------------------------------------------------------------------------------------------------------------------------------------------
summary(Y)


## ----plsda-pca, fig.cap='(ref:plsda-pca)'---------------------------------------------------------------------------------------------------------
pca.srbct <- pca(X, ncomp = 3, scale = TRUE)

plotIndiv(pca.srbct, group = srbct$class, ind.names = FALSE,
          legend = TRUE, 
          title = 'SRBCT, PCA comp 1 - 2')


## ----plsda-perf, fig.cap='(ref:plsda-perf)'-------------------------------------------------------------------------------------------------------
plsda.srbct <- plsda(X,Y, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda.srbct <- perf(plsda.srbct, validation = 'Mfold', folds = 3, 
                  progressBar = FALSE,  # Set to TRUE to track progress
                  nrepeat = 10)         # We suggest nrepeat = 50

plot(perf.plsda.srbct, sd = TRUE, legend.position = 'horizontal')


## ---- eval = FALSE--------------------------------------------------------------------------------------------------------------------------------
## perf.plsda.srbct


## -------------------------------------------------------------------------------------------------------------------------------------------------
final.plsda.srbct <- plsda(X,Y, ncomp = 3)


## ----plsda-indiv12, echo = TRUE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/PLSDA/'-------------------------------------------------
plotIndiv(final.plsda.srbct, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA on SRBCT comp 1-2',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')


## ----plsda-indiv13, echo = TRUE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/PLSDA/'-------------------------------------------------
plotIndiv(final.plsda.srbct, ind.names = FALSE, legend=TRUE,
          comp=c(2,3), ellipse = TRUE, 
          title = 'PLS-DA on SRBCT comp 2-3',
          X.label = 'PLS-DA comp 2', Y.label = 'PLS-DA comp 3')


## ----plsda-plotindiv, fig.cap='(ref:plsda-plotindiv)', echo=FALSE, eval=TRUE, fig.align='center', out.width = '50%', fig.height=4, fig.ncol = 2, fig.subcap=c('', '')----
knitr::include_graphics(c('Figures/PLSDA/plsda-indiv12-1.png', 'Figures/PLSDA/plsda-indiv13-1.png'))


## -------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.final.plsda.srbct <- perf(final.plsda.srbct, validation = 'Mfold', 
                               folds = 3, 
                               progressBar = FALSE, # TRUE to track progress
                               nrepeat = 50) 


## -------------------------------------------------------------------------------------------------------------------------------------------------
perf.final.plsda.srbct$error.rate$BER[, 'max.dist']


## -------------------------------------------------------------------------------------------------------------------------------------------------
perf.final.plsda.srbct$error.rate.class$max.dist


## ----plsda-background-max, echo = TRUE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/PLSDA/'------------------------------------------
background.max <- background.predict(final.plsda.srbct, 
                                     comp.predicted = 2,
                                     dist = 'max.dist') 

plotIndiv(final.plsda.srbct, comp = 1:2, group = srbct$class,
          ind.names = FALSE, title = 'Maximum distance',
          legend = TRUE,  background = background.max)


## ----plsda-background-cent, echo = FALSE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/PLSDA/'----------------------------------------
background.cent <- background.predict(final.plsda.srbct, 
                                      comp.predicted = 2,
                                      dist = 'centroids.dist') 

plotIndiv(final.plsda.srbct, comp = 1:2, group = srbct$class,
          ind.names = FALSE, title = 'Centroids distance',
          legend = TRUE,  background = background.cent)


## ----plsda-background-mah, echo = FALSE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/PLSDA/'-----------------------------------------
background.mah <- background.predict(final.plsda.srbct, 
                                     comp.predicted = 2,
                                     dist = 'mahalanobis.dist') 

plotIndiv(final.plsda.srbct, comp = 1:2, group = srbct$class,
          ind.names = FALSE, title = 'Mahalanobis distance',
          legend = TRUE,  background = background.mah)


## ----plsda-background, fig.cap='(ref:plsda-background)', echo=FALSE, eval=TRUE, fig.cap='(ref:plsda-background)', fig.align='center', out.width = '30%', fig.height=4, fig.ncol = 3, fig.subcap=c('', '', '')----
knitr::include_graphics(c('Figures/PLSDA/plsda-background-max-1.png', 'Figures/PLSDA/plsda-background-cent-1.png', 'Figures/PLSDA/plsda-background-mah-1.png'))


## -------------------------------------------------------------------------------------------------------------------------------------------------
# Grid of possible keepX values that will be tested for each comp
list.keepX <- c(1:10,  seq(20, 100, 10))
list.keepX


## -------------------------------------------------------------------------------------------------------------------------------------------------
# This chunk takes ~ 2 min to run
# Some convergence issues may arise but it is ok as this is run on CV folds
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, validation = 'Mfold', 
                                 folds = 5, dist = 'max.dist', 
                                 test.keepX = list.keepX, nrepeat = 10)


## -------------------------------------------------------------------------------------------------------------------------------------------------
# Just a head of the classification error rate per keepX (in rows) and comp
head(tune.splsda.srbct$error.rate)


## ----splsda-tune, fig.cap='(ref:splsda-tune)'-----------------------------------------------------------------------------------------------------
# To show the error bars across the repeats:
plot(tune.splsda.srbct, sd = TRUE)


## -------------------------------------------------------------------------------------------------------------------------------------------------
# The optimal number of components according to our one-sided t-tests
tune.splsda.srbct$choice.ncomp$ncomp

# The optimal keepX parameter according to minimal error rate
tune.splsda.srbct$choice.keepX


## -------------------------------------------------------------------------------------------------------------------------------------------------
# Optimal number of components based on t-tests on the error rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp 
ncomp

# Optimal number of variables to select
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  
select.keepX

splsda.srbct <- splsda(X, Y, ncomp = 3, keepX = select.keepX) 


## -------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(34)  # For reproducibility with this handbook, remove otherwise

perf.splsda.srbct <- perf(splsda.srbct, folds = 5, validation = "Mfold", 
                  dist = "max.dist", progressBar = FALSE, nrepeat = 10)

# perf.splsda.srbct  # Lists the different outputs
perf.splsda.srbct$error.rate


## -------------------------------------------------------------------------------------------------------------------------------------------------
perf.splsda.srbct$error.rate.class


## ----splsda-stability, fig.cap='(ref:splsda-stability)', results='hold', fig.show='hold'----------------------------------------------------------
par(mfrow=c(1,2))
# For component 1
stable.comp1 <- perf.splsda.srbct$features$stable$comp1
barplot(stable.comp1, xlab = 'variables selected across CV folds', 
        ylab = 'Stability frequency',
        main = 'Feature stability for comp = 1')

# For component 2
stable.comp2 <- perf.splsda.srbct$features$stable$comp2
barplot(stable.comp2, xlab = 'variables selected across CV folds', 
        ylab = 'Stability frequency',
        main = 'Feature stability for comp = 2')
par(mfrow=c(1,1))


## -------------------------------------------------------------------------------------------------------------------------------------------------
# First extract the name of selected var:
select.name <- selectVar(splsda.srbct, comp = 1)$name

# Then extract the stability values from perf:
stability <- perf.splsda.srbct$features$stable$comp1[select.name]

# Just the head of the stability of the selected var:
head(cbind(selectVar(splsda.srbct, comp = 1)$value, stability))


## ----splsda-indiv12, echo = TRUE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/PLSDA/'------------------------------------------------
plotIndiv(splsda.srbct, comp = c(1,2),
          ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          star = TRUE,
          title = 'SRBCT, sPLS-DA comp 1 - 2')


## ----splsda-indiv23, echo = TRUE, results = 'hide', fig.show = 'hide', fig.path = 'Figures/PLSDA/'------------------------------------------------
plotIndiv(splsda.srbct, comp = c(2,3),
          ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          star = TRUE,
          title = 'SRBCT, sPLS-DA comp 2 - 3')


## ----splsda-indiv, fig.cap='(ref:splsda-indiv)', echo=FALSE, eval=TRUE, fig.cap='(ref:splsda-indiv)', fig.align='center', out.width = '50%', fig.height=4, fig.ncol = 2, fig.subcap=c('(a)', '(b)')----
knitr::include_graphics(c('Figures/PLSDA/splsda-indiv12-1.png', 'Figures/PLSDA/splsda-indiv23-1.png'))


## ----splsda-var, fig.cap='(ref:splsda-var)'-------------------------------------------------------------------------------------------------------
var.name.short <- substr(srbct$gene.name[, 2], 1, 10)
plotVar(splsda.srbct, comp = c(1,2), 
        var.names = list(var.name.short), cex = 3)


## ----splsda-plotloading, fig.cap='(ref:splsda-plotloading)'---------------------------------------------------------------------------------------
plotLoadings(splsda.srbct, comp = 1, method = 'mean', contrib = 'max', 
             name.var = var.name.short)


## ----splsda-cim, fig.width=10, fig.height=8, fig.cap='(ref:splsda-cim)'---------------------------------------------------------------------------
cim(splsda.srbct, row.sideColors = color.mixo(Y))


## ---- results='hold'------------------------------------------------------------------------------------------------------------------------------
set.seed(33) # For reproducibility with this handbook, remove otherwise
train <- sample(1:nrow(X), 50)    # Randomly select 50 samples in training
test <- setdiff(1:nrow(X), train) # Rest is part of the test set

# Store matrices into training and test set:
X.train <- X[train, ]
X.test <- X[test,]
Y.train <- Y[train]
Y.test <- Y[test]

# Check dimensions are OK:
dim(X.train); dim(X.test)


## -------------------------------------------------------------------------------------------------------------------------------------------------
train.splsda.srbct <- splsda(X.train, Y.train, ncomp = 3, keepX = c(20,30,40))


## -------------------------------------------------------------------------------------------------------------------------------------------------
predict.splsda.srbct <- predict(train.splsda.srbct, X.test, 
                                dist = "mahalanobis.dist")


## -------------------------------------------------------------------------------------------------------------------------------------------------
# Just the head:
head(data.frame(predict.splsda.srbct$class, Truth = Y.test))


## -------------------------------------------------------------------------------------------------------------------------------------------------
# Compare prediction on the second component and change as factor
predict.comp2 <- predict.splsda.srbct$class$mahalanobis.dist[,2]
table(factor(predict.comp2, levels = levels(Y)), Y.test)


## -------------------------------------------------------------------------------------------------------------------------------------------------
# Compare prediction on the third component and change as factor
predict.comp3 <- predict.splsda.srbct$class$mahalanobis.dist[,3]
table(factor(predict.comp3, levels = levels(Y)), Y.test)


## -------------------------------------------------------------------------------------------------------------------------------------------------
# On component 3, just the head:
head(predict.splsda.srbct$predict[, , 3])


## ----splsda-roc, fig.cap='(ref:splsda-roc)', results='hold'---------------------------------------------------------------------------------------
auc.srbct <- auroc(splsda.srbct)

