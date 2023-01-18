
pages <- list.files(pattern = "\\.Rmd$") # extract all Rmds in this directory
pages <- pages[pages!="mixOmics Vignette.Rmd"] # remove single Rmd form from list
n_pages <- length(pages)
pages <- c(pages[n_pages], pages[-n_pages]) # ensure "index.Rmd" is first

all_lines <- c() # initialise output

for (page in pages) { # for each page ...
  
  lines <- readLines(page) # extract all lines from RMD
   
  # if the page is one of the following, find where the "fig-path" code chunk is.
  # this indicates the start of that file after the options chunk. Each pages options
  # chunk is not necessary as "index.Rmd" contains the global options chunk. 
  # Hence, remove all lines before the "fig-path" code chunk
  if (!(page %in% c("index.Rmd", "08-SessionInfo.Rmd", "09-References.Rmd"))) { 
    start_line <- which(grepl("fig-path", lines))
    lines <- lines[-(1:(start_line-1))]
  }
  
  if (page != "index.Rmd") { lines <- c(rep("", 10), lines) }
  
  all_lines <- c(all_lines, lines)
  
}

fileConn<-file("mixOmics Vignette.Rmd")
writeLines(all_lines, fileConn)
close(fileConn)
