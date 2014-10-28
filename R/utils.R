testConnection <- function() {
    con <- dbConnect(dbDriver("MySQL"), group = "jadb")
    ok <- dbExistsTable(con, 'labfiles')
    dbDisconnect(con)
    if( ok ) {
        message('Connection OK')
    } else {
        stop('DB connection problem. Configure jadb group in $HOME/.my.cnf.')
    }
}