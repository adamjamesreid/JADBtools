#' jar_process_exp
#'
#' @param IDS - ContactExpIDs
#' @param genome - Genome version
#'
#' @return NULL
#' @export
#'
#' @examples
#' #jar_process_exp( JADBtools:::make_jadb_ids(106:112, 'CG'), genome = 'cb3ce11' )
#' 
jar_process_exp <- function(IDS, genome='ce11') {
    if( genome=='ce10' ) {
        rcmd <- paste(sprintf(
            'JADBtools::jacl_send_to_cluster(\\"%s\\")', 
            IDS, 'cb3ce11'
        ), collapse = '; ')
    } else {
            rcmd <- paste(sprintf(
                'JADBtools::jacl_send_to_cluster(\\"%s\\", genome=\\"%s\\")', 
                IDS, genome
            ), collapse = '; ')
    }
    cmd <- sprintf('R -e \' library(JADBtools); %s \'', rcmd)
    z <- ssh.utils::run.remote(cmd, remote, verbose = T)
    cat(z$cmd.out, sep='\n')
    
    
}




#' jar_squeue
#'
#' @return NULL
#' @export
#'
jar_squeue <- function(tbl=TRUE, remote='jarun@cb-head2.gurdon.private.cam.ac.uk') {
    cmd <- 'squeue -u jarun'
    z <- ssh.utils::run.remote(cmd, remote, verbose = F)
    if(tbl) {
        read.table(textConnection(z$cmd.out), header = TRUE, stringsAsFactors = FALSE) %>% tbl_df()
    } else {
        cat(z$cmd.out, sep='\n')
    }
}

#' jar_log
#'
#' @param p - pattern to match
#' @param n - number of linest to print
#' @param db - database to search
#'
#' @return NULL
#' @export
#'
jar_log <- function(p='.chip', n=3, db='jadb', remote='jarun@cb-head2.gurdon.private.cam.ac.uk') {
    cmd <- sprintf('cd /mnt/%s/DBfile/DBfiles/_log && tail -n %i *%s*', db, n, p)
    z <- ssh.utils::run.remote(cmd, remote, verbose = F)
    cat(z$cmd.out, sep='\n')
}

#' jar_nlog
#'
#' @param p - pattern to match
#' @param n - number of linest to print
#' @param db - database to search
#'
#' @return NULL
#' @export
#'
jar_nlog <- function(p='.chip', n=3, db='jadb', remote='jarun@cb-head2.gurdon.private.cam.ac.uk') {
    cmd <- sprintf('cd /mnt/%s/DBfile/DBfiles/_log && ls *%s*', db, p)
    z <- ssh.utils::run.remote(cmd, remote, verbose = F)
    gsub(p, '', z$cmd.out)
}
