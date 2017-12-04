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
jadb_addAlignedBAM <- function(ids, processing='raw', format="txt.gz", genome='ce11') {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = format, eq=TRUE, processing = processing, scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    prefix <- gsub("\\^[^\\^]+$", '', basename(fls)) 
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    message('Aligning file: ', basename(fls))
    message(getwd())
    
    ALN <- run_bwa(
        basename(fls), 
        file.path(base_dir, sprintf('_ref_genomes_/%s/%s.fa', genome, genome))
    )
    message('Aligned!')
    
    stats <- JADBtools::bamStats(ALN)
    message(stats)
    
    final.path <- file.path('files', exp_dir, gsub('raw\\^NA\\^NA', paste0('aligned^NA^', genome), prefix))
    
    Entry <- addGenericFile(
        ids,
        path = final.path, 
        Processing = 'aligned', 
        Scale = 'NA', 
        Resolution = 'NA',
        filetype_format = 'bam', 
        prefix = 'A',
        comments = stats,
        genome = genome,
        uniq = FALSE
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
    # bwa mem -t 24 /mnt/jadb/DBfile/DBfiles/_ref_genomes_/ce10.fa LEM2^Novus48540002_GSM1911031^FALSE^N2^EE_raw^NA^NA_GSM1911031^Cfd37769.txt.gz >
    read_length <- median(nchar(readLines(file, 4000)[seq(2, 4000, by=4)]))
    if(read_length < 70) {
        message('using BWA backtrack')
        cmd <- sprintf(
            #"bwa samse %s <(bwa aln -t %i %s %s) %s | samtools view -@ %i -Su -  | samtools sort -@ %i - %s",
            "bwa samse %s <(bwa aln -t %i %s %s) %s | samtools view -@ %i -bSu - | samtools sort -@ %i -m 44294967296 - -o %s.bam",
            genome, ncore, genome, file, file, ncore, ncore, output
        )
    } else {
        message('using BWA MEM')
        cmd <- sprintf(
            "bwa mem -t %i %s %s | samtools view -@ %i -bSu - | samtools sort -@ %i -m 44294967296 - -o %s.bam",
            ncore, genome, file, ncore, ncore, output
        )
    }
    
    cmd2 <- sprintf('echo "%s" | %s', cmd, interperor)
    message(cmd2)
    system(cmd2)
    return(paste0(output, '.bam'))
}

run_fastqc <- function(file, interperor='bash') {
    
    output <- gsub('\\..+$', '', basename(file))
    #parallel echo ...  ::: *.fastq
    # ref=/mnt/home1/ahringer/ps562/_ref_genomes_/ce10.fa
    cmd <- sprintf(
        "fastqc %s", file
    )
    cmd2 <- sprintf('echo "%s" | %s', cmd, interperor)
    message(cmd2)
    system(cmd2)
    return(paste0(output, '_fastqc.html'))
}

run_fastq_screen <- function(file, interperor='bash') {
    
    output <- gsub('\\..+$', '', basename(file))
    cmd <- sprintf(
        "fastq_screen --aligner bwa %s", file
    )
    cmd2 <- sprintf('echo "%s" | %s', cmd, interperor)
    message(cmd2)
    system(cmd2)
    return(paste0(output, '_screen.html'))
}

run_star <- function(file, interperor='bash') {
    
    output <- gsub('_.{8}.txt.gz', '_', basename(file))
    gndir <- '/mnt/jadb/DBfile/DBfiles/_ref_genomes_/ce11_star'
    
    cmd <- sprintf(
        'STAR --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN 32 --genomeDir %s --readFilesIn %s --outFileNamePrefix ./%s',
        gndir, file, output
    )
    
    cmd2 <- sprintf('echo "%s" | %s', cmd, interperor)
    message(cmd2)
    system(cmd2)
    return(paste0(output, 'Aligned.sortedByCoord.out.bam'))
    
}


#' jadb_processs_sceleton
#' 
#' @param IDs Vector of JADB ContactExpIDs
#'   
#' @return List 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family internal_jadb
#' @export
#' 
#' @examples
#' #jadb_processs_sceleton(
#' #    ids, FUN = run_fastqc, 
#' #    Processing = 'FastQC', Scale = 'NA', Resolution = 'NA', filetype_format = 'HTML',
#' #    format = "bam", processing = "aligned", scale = "NA", resolution='NA'
#' #)
jadb_processs_sceleton <- function(
    ids, FUN, 
    Processing = 'TETS', Scale = 'NA', Resolution = 'NA', filetype_format = 'unknown',
    format = "bam", processing = "aligned", scale = "NA", resolution='NA', 
    genome = "NA", uniq='NA', rename_output=TRUE, prefix_output='P', ...
) {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    #Pre-setup
    base_dir  <- getwd(); on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = format, eq=TRUE, processing = processing, scale = scale, url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 files.')
    
    prefix <- gsub("\\^[^\\^]+$", '', basename(fls)) 
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir); setwd(exp_dir)
    
    message('Processing file: ', basename(fls))
    message(getwd())
    
    final.path <- file.path('files', exp_dir, gsub(
        sprintf('%s.%s.%s.+', processing, scale, resolution), 
        sprintf('%s^%s^%s', Processing, Scale, genome), 
        prefix
    ))
    message('Output file path prefix: ', final.path)

    OUT <- FUN(basename(fls))
    message('Processing done!')

    Entry <- addGenericFile(
        ids,
        path = if(rename_output) final.path else file.path('files', exp_dir, OUT), 
        Processing = Processing, 
        Scale = Scale, 
        Resolution = Resolution,
        filetype_format = filetype_format, 
        prefix = prefix_output,
        comments = '',
        genome = genome,
        uniq = uniq
    )
    if(rename_output) {
        out <- file.rename(basename(OUT), basename(Entry$path))
    } else {
        out <- file.rename(OUT, file.path(dirname(OUT), basename(Entry$path)))
    }
    
    message('Result saved as ', basename(Entry$path))
    message("Done!")
    
}

