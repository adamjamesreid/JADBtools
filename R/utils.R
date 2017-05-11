#' Test connection to JADB
#'   
#' @return message or error
#' 
#' @author Przemyslaw Stempor
#' 
#' @family dbutils
#' @export
#' 
testConnection <- function() {
    con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf')
    ok <- dbExistsTable(con, 'labfiles')
    message('Connection status:')
    message(paste0('    ', names(dbGetInfo(con)), ': ', dbGetInfo(con), collapse = '\n'))
    dbDisconnect(con)
    if( ok ) {
        message('Connection OK')
    } else {
        stop('DB connection problem. Configure jadb group in $HOME/.my.cnf.')
    }
}

#' Print progress bar JSON
#'   
#' @return message or error
#' 
#' @author Przemyslaw Stempor
#' 
#' @family dbutils
#' @export
#' 
pb <- function(
    mnt, uid, s=TRUE, p=0.01, t='Loading libraries...', tt=''
) {
    cat( toJSON(
        list(success=s, progress=p, text=t, toptext=tt)
    ), file=paste(mnt, '/temp/', uid, '.json', sep='') )
}