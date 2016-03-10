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
#' @examples 
#' #final_out <- lapply(dir(pattern='*.bw'), renderReport)
#' 
renderReport <- function(infile) {
    message('==> Processing: ', infile)
    require(rmarkdown)
    
    require(Cairo)
    options("bitmapType" = "cairo")
    
    gsub('bw', '_cache', basename(infile))
    
    file.copy(
        system.file('reporting_templates/ChIPenrichmentOnReporting.Rmd', package = 'JADBtools'), 
        gsub('bw$', 'Rmd', basename(infile))
    )
    
    rmarkdown::render(
        'ChIPenrichmentOnReporting.Rmd',
        params = list(inputfile = infile), 
        output_file = gsub('bw', 'html', basename(infile))
    )
}