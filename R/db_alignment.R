#' jadb_addAlignedBAM
#' 
#' @param IDs Vector of JADB ContactExpIDs
#'   
#' @return List 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family tracks
#' @export
#' 
#' @examples
#' #jadb_addAlignedBAM('AA001')
#' #parallel::mclapply(sprintf('REP%.3i', 36:42)[-6], jadb_addAlignedBAM), mc.cores = 8)
jadb_addAlignedBAM <- function(ids) {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "txt.gz", eq=TRUE, processing = "raw", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    prefix <- gsub("\\^[^\\^]+$", '', basename(fls)) 
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    message('Aligning file: ', basename(fls))
    message(getwd())
    
    ALN <- run_bwa(
        basename(fls), 
        file.path(base_dir, '_ref_genomes_/ce10.fa')
    )
    message('Aligned!')
    
    stats <- JADBtools::bamStats(ALN)
    message(stats)
    
    final.path <- file.path('files', exp_dir, gsub('raw\\^NA\\^NA', 'aligned^NA^NA', prefix))
    
    Entry <- addGenericFile(
        ids,
        path = final.path, 
        Processing = 'aligned', 
        Scale = 'NA', 
        Resolution = 'NA',
        filetype_format = 'bam', 
        prefix = 'P',
        comments = stats
    )
    
    out <- file.rename(basename(ALN), basename(Entry$path))
    
    message("Indexing")
    indexBam(basename(Entry$path))
    
    message("Done!")
    
}

run_bwa <- function(file, genome, ncore=8, interperor='bash') {
    
    output <- gsub('\\..+$', '', basename(file))
    #parallel echo ...  ::: *.fastq
    # ref=/mnt/home1/ahringer/ps562/_ref_genomes_/ce10.fa
    cmd <- sprintf(
        "bwa samse %s <(bwa aln -t %i %s %s) %s | samtools view -@ %i -Su - | samtools sort -@ %i - %s",
        genome, ncore, genome, file, file, ncore, ncore, output
    )
    cmd2 <- sprintf('echo "%s" | %s', cmd, interperor)
    message(cmd2)
    system(cmd2)
    return(paste0(output, '.bam'))
}

