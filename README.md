# mixOmics vignette files

All Rproject files are stored in `./vignetteRproject`. In addition:

   * `./BiocStyle`:  Contains the vignette files in `BiocStyle` format.
   * `./pdf_book` :  Contains the vignette files in `pdf`          format.
   * `./docs`:       Contains the vignette files in `html`         format.

## How to make the vignettes

   * open the .Rproj file in RStudio
   * **for gitbook** from the `build` pane in RStudio, from `Build Book` dropdown choose `bookdown`.
   * **for BiocStyle** in `index.Rmd` change the author style to `BiocStyle` and then choose `BiocStyle`.
   * **for pdf** ensure appropriate `peamble.tex` is included and choose `pdf`
   * move `Vignette-mixOmics.html` to `./BiocStyle`.
   * move the pdf file vignette to `./pdf_book`.
   * The gitbook/html vignette will now be updated in `./docs`.
   
## Bookdown link
https://mixomicsteam.github.io/Bookdown/


## TODO

   * The `pdf` output needs improvement (if it is actually needed, as it is not needed by Bioc.)
   * html headers tailored to output that include authors and affiliations