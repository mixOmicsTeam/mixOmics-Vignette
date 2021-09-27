## ----global_options, include=FALSE----------------------------------------------------------------------------------------------------------------
library(knitr)
# global options
knitr::opts_chunk$set(dpi = 100, echo=TRUE, warning=FALSE, message=FALSE, eval = TRUE,
                      fig.show=TRUE, fig.width= 7,fig.height= 6,fig.align='center', out.width = '50%', message = FALSE,
                      fig.path= 'Figures/PLS/')

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


## -------------------------------------------------------------------------------------------------------------------------------------------------
library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic


## -------------------------------------------------------------------------------------------------------------------------------------------------
head(data.frame(rownames(X), rownames(Y)))


## -------------------------------------------------------------------------------------------------------------------------------------------------
y <- liver.toxicity$clinic[, "ALB.g.dL."]


## ----pls1-Q2, eval = TRUE, fig.cap='(ref:pls1-Q2)'------------------------------------------------------------------------------------------------
tune.pls1.liver <- pls(X = X, Y = y, ncomp = 4, mode = 'regression')
set.seed(33)  # For reproducibility with this handbook, remove otherwise
Q2.pls1.liver <- perf(tune.pls1.liver, validation = 'Mfold', 
                      folds = 10, nrepeat = 5)
plot(Q2.pls1.liver, criterion = 'Q2')


## ----spls1-MAE, fig.cap='(ref:spls1-MAE)'---------------------------------------------------------------------------------------------------------
# Set up a grid of values: 
list.keepX <- c(5:10, seq(15, 50, 5))     

# list.keepX  # Inspect the keepX grid
set.seed(33)  # For reproducibility with this handbook, remove otherwise
tune.spls1.MAE <- tune.spls(X, y, ncomp= 2, 
                            test.keepX = list.keepX, 
                            validation = 'Mfold', 
                            folds = 10,
                            nrepeat = 5, 
                            progressBar = FALSE, 
                            measure = 'MAE')
plot(tune.spls1.MAE)


## -------------------------------------------------------------------------------------------------------------------------------------------------
choice.ncomp <- tune.spls1.MAE$choice.ncomp$ncomp
# Optimal number of variables to select in X based on the MAE criterion
# We stop at choice.ncomp
choice.keepX <- tune.spls1.MAE$choice.keepX[1:choice.ncomp]  

choice.ncomp
choice.keepX


## -------------------------------------------------------------------------------------------------------------------------------------------------
spls1.liver <- spls(X, y, ncomp = choice.ncomp, keepX = choice.keepX, 
                    mode = "regression")


## ---- eval = FALSE--------------------------------------------------------------------------------------------------------------------------------
## selectVar(spls1.liver, comp = 1)$X$name


## -------------------------------------------------------------------------------------------------------------------------------------------------
spls1.liver$prop_expl_var$X
tune.pls1.liver$prop_expl_var$X


## ----spls1-ext, fig.cap='(ref:spls1-ext)', message=FALSE------------------------------------------------------------------------------------------
spls1.liver.c2 <- spls(X, y, ncomp = 2, keepX = c(rep(choice.keepX, 2)), 
                   mode = "regression")

plotIndiv(spls1.liver.c2,
          group = liver.toxicity$treatment$Time.Group,
          pch = as.factor(liver.toxicity$treatment$Dose.Group),
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')



## ----spls1-comp1, fig.cap='(ref:spls1-comp1)'-----------------------------------------------------------------------------------------------------
# Define factors for colours matching plotIndiv above
time.liver <- factor(liver.toxicity$treatment$Time.Group, 
                     levels = c('18', '24', '48', '6'))
dose.liver <- factor(liver.toxicity$treatment$Dose.Group, 
                     levels = c('50', '150', '1500', '2000'))
# Set up colours and symbols
col.liver <- color.mixo(time.liver)
pch.liver <- as.numeric(dose.liver)

plot(spls1.liver$variates$X, spls1.liver$variates$Y,
     xlab = 'X component', ylab = 'y component / scaled y',
     col = col.liver, pch = pch.liver)
legend('topleft', col = color.mixo(1:4), legend = levels(time.liver),
       lty = 1, title = 'Time')
legend('bottomright', legend = levels(dose.liver), pch = 1:4,
       title = 'Dose')

cor(spls1.liver$variates$X, spls1.liver$variates$Y)


## ---- eval = TRUE---------------------------------------------------------------------------------------------------------------------------------
set.seed(33)  # For reproducibility with this handbook, remove otherwise

# PLS1 model and performance
pls1.liver <- pls(X, y, ncomp = choice.ncomp, mode = "regression")
perf.pls1.liver <- perf(pls1.liver, validation = "Mfold", folds =10, 
                   nrepeat = 5, progressBar = FALSE)
perf.pls1.liver$measures$MSEP$summary
# To extract values across all repeats:
# perf.pls1.liver$measures$MSEP$values

# sPLS1 performance
perf.spls1.liver <- perf(spls1.liver, validation = "Mfold", folds = 10, 
                   nrepeat = 5, progressBar = FALSE)
perf.spls1.liver$measures$MSEP$summary


## -------------------------------------------------------------------------------------------------------------------------------------------------
dim(Y)


## ----pls2-Q2, fig.cap='(ref:pls2-Q2)'-------------------------------------------------------------------------------------------------------------

tune.pls2.liver <- pls(X = X, Y = Y, ncomp = 5, mode = 'regression')

set.seed(33)  # For reproducibility with this handbook, remove otherwise
Q2.pls2.liver <- perf(tune.pls2.liver, validation = 'Mfold', folds = 10, 
                      nrepeat = 5)
plot(Q2.pls2.liver, criterion = 'Q2.total')


## ---- echo = TRUE, eval = TRUE, warning=FALSE-----------------------------------------------------------------------------------------------------
# This code may take several min to run, parallelisation option is possible
list.keepX <- c(seq(5, 50, 5))
list.keepY <- c(3:10)

set.seed(33)  # For reproducibility with this handbook, remove otherwise
tune.spls.liver <- tune.spls(X, Y, test.keepX = list.keepX, 
                             test.keepY = list.keepY, ncomp = 2, 
                             nrepeat = 1, folds = 10, mode = 'regression', 
                             measure = 'cor', 
                            #   the following uses two CPUs for faster computation
                            # it can be commented out
                            BPPARAM = BiocParallel::SnowParam(workers = 2)
                            )


## ----spls2-tune, fig.cap='(ref:spls2-tune)', out.width = '60%'------------------------------------------------------------------------------------
plot(tune.spls.liver)


## -------------------------------------------------------------------------------------------------------------------------------------------------
#Optimal parameters
choice.keepX <- tune.spls.liver$choice.keepX
choice.keepY <- tune.spls.liver$choice.keepY
choice.ncomp <- length(choice.keepX)

spls2.liver <- spls(X, Y, ncomp = choice.ncomp, 
                   keepX = choice.keepX,
                   keepY = choice.keepY,
                   mode = "regression")


## -------------------------------------------------------------------------------------------------------------------------------------------------
spls2.liver$prop_expl_var


## ---- eval = FALSE--------------------------------------------------------------------------------------------------------------------------------
## selectVar(spls2.liver, comp = 1)$X$value


## -------------------------------------------------------------------------------------------------------------------------------------------------
vip.spls2.liver <- vip(spls2.liver)
# just a head
head(vip.spls2.liver[selectVar(spls2.liver, comp = 1)$X$name,1])


## ---- results = 'hide'----------------------------------------------------------------------------------------------------------------------------
perf.spls2.liver <- perf(spls2.liver, validation = 'Mfold', folds = 10, nrepeat = 5)
# Extract stability
stab.spls2.liver.comp1 <- perf.spls2.liver$features$stability.X$comp1
# Averaged stability of the X selected features across CV runs, as shown in Table
stab.spls2.liver.comp1[1:choice.keepX[1]]

# We extract the stability measures of only the variables selected in spls2.liver
extr.stab.spls2.liver.comp1 <- stab.spls2.liver.comp1[selectVar(spls2.liver, 
                                                                  comp =1)$X$name]


## ----pls2stability, echo = FALSE------------------------------------------------------------------------------------------------------------------
knitr::kable(extr.stab.spls2.liver.comp1, caption = 'Stability measure (occurence of selection) of the variables from X selected with sPLS2 across repeated 10-fold subsampling on component 1.', longtable = TRUE)


## ----spls2-plotIndiv, fig.cap='(ref:spls2-plotIndiv)'---------------------------------------------------------------------------------------------
plotIndiv(spls2.liver, ind.names = FALSE, 
          group = liver.toxicity$treatment$Time.Group, 
          pch = as.factor(liver.toxicity$treatment$Dose.Group), 
          col.per.group = color.mixo(1:4),
          legend = TRUE, legend.title = 'Time', 
          legend.title.pch = 'Dose')


## ----spls2-plotArrow, fig.cap='(ref:spls2-plotArrow)'---------------------------------------------------------------------------------------------
plotArrow(spls2.liver, ind.names = FALSE, 
          group = liver.toxicity$treatment$Time.Group,
          col.per.group = color.mixo(1:4),
          legend.title = 'Time.Group')


## ----spls2-plotVar, fig.cap='(ref:spls2-plotVar)'-------------------------------------------------------------------------------------------------
plotVar(spls2.liver, cex = c(3,4), var.names = c(FALSE, TRUE))


## ----spls2-plotVar2, fig.cap='(ref:spls2-plotVar2)'-----------------------------------------------------------------------------------------------
plotVar(spls2.liver,
        var.names = list(X.label = liver.toxicity$gene.ID[,'geneBank'],
        Y.label = TRUE), cex = c(3,4))


## -------------------------------------------------------------------------------------------------------------------------------------------------
# Define red and green colours for the edges
color.edge <- color.GreenRed(50)

# X11()  # To open a new window for Rstudio
network(spls2.liver, comp = 1:2,
        cutoff = 0.7,
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        # To save the plot, otherwise comment out:
        save = 'pdf', name.save = 'network_liver')


## ----spls2-network, eval=TRUE, echo=FALSE,fig.cap='(ref:spls2-network)'---------------------------------------------------------------------------
knitr::include_graphics("network_liver.pdf")


## -------------------------------------------------------------------------------------------------------------------------------------------------
# X11()  # To open a new window if the graphic is too large
cim(spls2.liver, comp = 1:2, xlab = "clinic", ylab = "genes",
    # To save the plot, otherwise comment out:
    save = 'pdf', name.save = 'cim_liver')


## ----spls2-cim, eval=TRUE, echo=FALSE, fig.cap='(ref:spls2-cim)'----------------------------------------------------------------------------------
knitr::include_graphics("cim_liver.pdf")


## ---- eval = TRUE---------------------------------------------------------------------------------------------------------------------------------
# Comparisons of final models (PLS, sPLS)

## PLS
pls.liver <- pls(X, Y, mode = 'regression', ncomp = 2)
perf.pls.liver <-  perf(pls.liver, validation = 'Mfold', folds = 10, 
                        nrepeat = 5)

## Performance for the sPLS model ran earlier
perf.spls.liver <-  perf(spls2.liver, validation = 'Mfold', folds = 10, 
                         nrepeat = 5)


## ----pls-perf-spls2, fig.cap='(ref:pls-perf-spls2)', out.width = '70%'----------------------------------------------------------------------------
plot(c(1,2), perf.pls.liver$measures$cor.upred$summary$mean, 
     col = 'blue', pch = 16, 
     ylim = c(0.6,1), xaxt = 'n',
     xlab = 'Component', ylab = 't or u Cor', 
     main = 's/PLS performance based on Correlation')
axis(1, 1:2)  # X-axis label
points(perf.pls.liver$measures$cor.tpred$summary$mean, col = 'red', pch = 16)
points(perf.spls.liver$measures$cor.upred$summary$mean, col = 'blue', pch = 17)
points(perf.spls.liver$measures$cor.tpred$summary$mean, col = 'red', pch = 17)
legend('bottomleft', col = c('blue', 'red', 'blue', 'red'), 
       pch = c(16, 16, 17, 17), c('u PLS', 't PLS', 'u sPLS', 't sPLS'))

