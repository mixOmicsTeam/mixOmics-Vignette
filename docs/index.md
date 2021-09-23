--- 
title: 'Hands-on activities'
author:
- name: Kim-Anh Le Cao
  email: kimanh.lecao@unimelb.edu.au
  affiliation: Melbourne Integrative Genomics & School of Mathematics and Statistics, The University of Melbourne, Australia

package: mixOmics
date: '\today'
site: bookdown::bookdown_site
output: 
  bookdown::gitbook:
    split_bib: no
    includes:
     in_header: header.html
documentclass: book
bibliography: ["mybib3.bib"]
biblio-style: apalike
link-citations: true
github-repo: mixOmicsTeam/mixOmics
description: "Hands-on activities using the R package mixOmics"
---





# Preamble {-}

The following vignettes serve as support for our [online course](https://study.unimelb.edu.au/find/short-courses/mixomics-r-essentials-for-biological-data-integration/#course-specifics). More information about the methods can be found at [www.mixOmics.org](www.mixOmics.org)


![(\#fig:methods-fig)**Different types of analyses with mixOmics** [@mixomics].The biological questions, the number of data sets to integrate, and the type of response variable, whether qualitative (classification), quantitative (regression), one (PLS1) or several (PLS) responses, all drive the choice of analytical method. All methods featured in this diagram include variable selection except rCCA. In N-integration, rCCA and PLS enable the integration of two quantitative data sets, whilst the block PLS methods (that derive from the methods from @Ten11) can integrate more than two data sets. In P-integration, our method MINT is based on multi-group PLS [@Esl14b].The following activities cover some of these methods.](XtraFigs/MixOmicsAnalysesV2.pdf) 



