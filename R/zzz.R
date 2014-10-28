.onAttach <- function(libname, pkgname) {
    #Fix for annoying BiocGenerics namespace problems
    if(!"package:BiocGenerics" %in% search()) attachNamespace('BiocGenerics')
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb")
    ok <- dbExistsTable(con, 'labfiles')
    dbDisconnect(con)
    if( ok ) {
        packageStartupMessage('Package loaded, DB connection OK.')
    } else {
        stop('DB connection problem. Configure jadb group in $HOME/.my.cnf.')
    }
    
}



