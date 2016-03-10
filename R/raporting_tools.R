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
    
    gsub('bw', '_cache', basename(infile))
    
    require(knitr)
    opts_chunk$set(cache.path = gsub('bw', '_cache', basename(infile)))
    
    rmarkdown::render(
        system.file('reporting_templates/ChIPenrichmentOnReporting.Rmd', package = 'JADBtools'),
        params = list(inputfile = infile), 
        output_file = gsub('bw', 'html', basename(infile))
    )
}