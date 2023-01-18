--- 
title: 'mixOmics vignette'
author:
- name: Kim-Anh Le Cao
  email: kimanh.lecao@unimelb.edu.au
  affiliation: Melbourne Integrative Genomics & School of Mathematics and Statistics, The University of Melbourne, Australia

package: mixOmics
site: bookdown::bookdown_site
output: 
  bookdown::gitbook:
    split_bib: no
    includes:
     in_header: header.html
documentclass: book
bibliography: ["bibliography.bib"]
biblio-style: apalike
link-citations: true
github-repo: mixOmicsTeam/mixOmics
description: "Vignette for the R package mixOmics"

vignette: >
  %\VignetteIndexEntry{mixOmics}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---






# Preamble {-}

If you are following our [online course](https://study.unimelb.edu.au/find/short-courses/mixomics-r-essentials-for-biological-data-integration/#course-specifics), the following vignette will be useful as a complementary learning tool. This vignette also covers the essential use cases of various methods in this package for the general `mixOmcis` user. The below methods will be covered:

- (s)PCA, 
- PLS1 and PLS2, 
- (s)PLS-DA,
- N-integration (multi-block sPLS-DA, aka. "DIABLO"), and 
- P-integration (multi-group sPLS-DA, aka "MINT").

As outlined in [1.3](#01:outline), this is not an exhaustive list of all the methods found within `mixOmics`. More information can be found at [our website](http://mixomics.org/) and you can ask questions via our [discourse forum](https://mixomics-users.discourse.group/).


<div class="figure" style="text-align: center">
<img src="InputFigures/MixOmicsAnalysesV2.png" alt="**Different types of analyses with mixOmics** [@mixomics].The biological questions, the number of data sets to integrate, and the type of response variable, whether qualitative (classification), quantitative (regression), one (PLS1) or several (PLS) responses, all drive the choice of analytical method. All methods featured in this diagram include variable selection except rCCA. In N-integration, rCCA and PLS enable the integration of two quantitative data sets, whilst the block PLS methods (that derive from the methods from @Ten11) can integrate more than two data sets. In P-integration, our method MINT is based on multi-group PLS [@Esl14b].The following activities cover some of these methods." width="75%"  />
<p class="caption">(\#fig:00-analyses-diagram)**Different types of analyses with mixOmics** [@mixomics].The biological questions, the number of data sets to integrate, and the type of response variable, whether qualitative (classification), quantitative (regression), one (PLS1) or several (PLS) responses, all drive the choice of analytical method. All methods featured in this diagram include variable selection except rCCA. In N-integration, rCCA and PLS enable the integration of two quantitative data sets, whilst the block PLS methods (that derive from the methods from @Ten11) can integrate more than two data sets. In P-integration, our method MINT is based on multi-group PLS [@Esl14b].The following activities cover some of these methods.</p>
</div>



