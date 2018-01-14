#' jadb_get_file_paths
#'
#' @param ids ContactExpIDs
#' @param processing raw
#' @param root_prefix /mnt/jadb/DBfile/DBfiles
#'
#' @return paths
#' @export
#'
jadb_get_file_paths <- function(
    ids, processing='raw', root_prefix='/mnt/jadb/DBfile/DBfiles'
) {
    
    if(nchar(root_prefix)) {
        pth <- sapply(ids, getFilePath, processing=processing, url=FALSE)
        pth <- gsub('^files', root_prefix, pth)
    } else {
        pth <- sapply(ids, getFilePath, processing=processing)
    }
    return(pth)
}


#' jadb_download_files
#'
#' @param ids ContactExpIDs
#' @param processing raw
#'
#' @return TRUE if downloaded
#' @export
#'
jadb_download_files <- function(ids, processing='raw', auto=FALSE) {
    
    pth <- jadb_get_file_paths(ids, root_prefix = '', processing=processing)
    
    if(!auto) {
        message('Will download ([y]es/[n]o):\n- ', paste(pth, collapse = '\n- '))
        if(scan(nmax = 1, what = 'character', quiet = TRUE)!='y') return()
    }

    sapply(pth, function(x) {
        if(file.exists(basename(x))) {
            message('Skupping ', basename(x), ' - file exist.')
        } else {
            download.file(x, basename(x))
        }
        
    })
}


#' make_jadb_ids
#'
#' @param ran range
#' @param prefix prefix
#'
#' @return ids
#' @export
#'
make_jadb_ids <- function(ran, prefix='rYD') {
    sprintf('%s%0.3i', prefix, ran)
}
    
