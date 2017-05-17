.onAttach <- function(libname, pkgname) {
    #Fix for annoying BiocGenerics namespace problems
    if(!"package:BiocGenerics" %in% search()) attachNamespace('BiocGenerics')
    
    if(nchar(Sys.getenv('JADB_GROUP'))) {
        GROUP <<- Sys.getenv('JADB_GROUP')
    } else {
        GROUP <<- 'jadb'
    }
    
    if(nchar(Sys.getenv('JADB_GROUP'))) {
        MOUNT <<- Sys.getenv('JADB_MOUNT')
    } else {
        MOUNT <<- '/mnt/jadb/DBfile/DBfiles'
    }
    message('Mount point:\n    ', MOUNT)
    message('    Exists: ', file.exists(MOUNT))
    testConnection()
    
    
    packageStartupMessage('Package loaded.')
}



