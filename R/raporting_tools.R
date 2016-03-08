#' ChIP profiles enrichment of repetitive elements in C. elegans raport
#' 
#' @param infile BigWig file path 
#'   
#' @return HTML file
#' 
#' @author Przemyslaw Stempor
#' 
#' @family raporting
#' @export
#' 
renderReport <- function(infile) {
    message('==> Processing: ', infile)
    require(rmarkdown)
    
    require(Cairo)
    options("bitmapType" = "cairo")
    
    rmarkdown::render(
        "ChIPenrichmentOnRepeats.Rmd",
        params = list(inputfile = infile), 
        output_file = gsub('bw', 'html', basename(infile))
    )
}