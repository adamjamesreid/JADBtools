#' addNonUniqueQ10Beads
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
#' #addNonUniqueQ10Beads('AA001')
addNonUniqueQ10Beads <- function(ids, genome='ce11') {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    outnames <- sapply(basename(fls), rbeads:::reName, proccesing = 'mapq0', scale = 'linear', resolution = genome, ext='.bw')
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    crosslink <- JADBtools::getAnno(ids, anno = 'Crosslinker', EXTABLE = 'labexperimentview')
    input_suffix <- '_HiSeq_nonUNIQ_MAPQ10_200bp_SummedInput_bin25bp.bw'
    if(grepl('^e', crosslink, ignore.case = TRUE)) {
        input <- file.path(base_dir, 'Input/SummedInputs', genome, paste0(genome, '_EGS', input_suffix))
    }  else {
        input <- file.path(base_dir, 'Input/SummedInputs', genome, paste0(genome, '_FRM', input_suffix))
    }
    
    mappability <- file.path(base_dir, paste0("_mappability_files_/", genome, "_gem-mappability_36bp.bw"))
    
    message('File: ', basename(fls), '\n vs. ', basename(input))
    
    message(getwd())
    NRM0 <- beads(
        basename(fls), 
        input, 
        mappability,
        if(genome == 'cb3ce11') file.path(base_dir, '_ref_genomes_/cb3ce11/cb3ce11.fa') else genome, 
        uniq = FALSE, insert = 200L, mapq_cutoff = 10, export = "BEADS", 
        rdata = FALSE, export_er = TRUE, quickMap = TRUE
    )
    
    
    final.path <- file.path('files', exp_dir, gsub('aligned\\^[A-Za-z0-9]+\\^[A-Za-z0-9]+', paste0('BEADSQ10NU^linear^', genome), prefix))
    
    Entry <- addGenericFile(
        ids,
        path = final.path, 
        Processing = 'BEADSQ10NU', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'B',
        genome = genome,
        comments = JADBtools::bamStats(basename(fls)),
        uniq = FALSE
    )
    
    out <- file.rename(basename(path(NRM0)), basename(Entry$path))
    
    
    message("Done!")
    
}


#' addNonUniqueQ10Beads
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
#' #addNonUniqueQ10Beads('AA001')
addNonUniqueNQBeads <- function(ids, genome='ce11') {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    crosslink <- JADBtools::getAnno(ids, anno = 'Crosslinker', EXTABLE = 'labexperimentview')
    input_suffix <- '_HiSeq_nonUNIQ_nonMAPQ_200bp_SummedInput_bin25bp.bw'
    if(grepl('^e', crosslink, ignore.case = TRUE)) {
        input <- file.path(base_dir, 'Input/SummedInputs', genome, paste0(genome, '_EGS', input_suffix))
    }  else {
        input <- file.path(base_dir, 'Input/SummedInputs', genome, paste0(genome, '_FRM', input_suffix))
    }
    
    mappability <- file.path(base_dir, paste0("_mappability_files_/", genome, "_gem-mappability_36bp.bw"))
    
    message('File: ', basename(fls), '\n vs. ', basename(input))
    
    message(getwd())
    NRM0 <- beads(
        basename(fls), 
        input, 
        mappability,
        if(genome == 'cb3ce11') file.path(base_dir, '_ref_genomes_/cb3ce11/cb3ce11.fa') else genome, 
        uniq = FALSE, insert = 200L, mapq_cutoff = 0, export = "BEADS", 
        rdata = FALSE, export_er = TRUE, quickMap = TRUE
    )
    
    
    final.path <- file.path('files', exp_dir, gsub('aligned\\^[A-Za-z0-9]+\\^[A-Za-z0-9]+', paste0('BEADSNQNU^linear^', genome), prefix))
    
    Entry <- addGenericFile(
        ids,
        path = final.path, 
        Processing = 'BEADSNQNU', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'B',
        genome = genome,
        comments = JADBtools::bamStats(basename(fls)),
        uniq = FALSE
    )
    
    out <- file.rename(basename(path(NRM0)), basename(Entry$path))
    
    
    message("Done!")
    
}



#' addUniqueQ10Beads
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
#' #addNonUniqueQ10Beads('AA001')
addUniqueQ10Beads <- function(ids, genome='ce11') {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    crosslink <- JADBtools::getAnno(ids, anno = 'Crosslinker', EXTABLE = 'labexperimentview')
    input_suffix <- '_HiSeq_UNIQ_MAPQ10_200bp_SummedInput_bin25bp.bw'
    if(grepl('^e', crosslink, ignore.case = TRUE)) {
        input <- file.path(base_dir, 'Input/SummedInputs', genome, paste0(genome, '_EGS', input_suffix))
    }  else {
        input <- file.path(base_dir, 'Input/SummedInputs', genome, paste0(genome, '_FRM', input_suffix))
    }
    
    mappability <- file.path(base_dir, paste0("_mappability_files_/", genome, "_gem-mappability_36bp.bw"))
    
    message('File: ', basename(fls), '\n vs. ', basename(input))
    
    message(getwd())
    NRM0 <- beads(
        basename(fls), 
        input, 
        mappability,
        if(genome == 'cb3ce11') file.path(base_dir, '_ref_genomes_/cb3ce11/cb3ce11.fa') else genome, 
        uniq = TRUE, insert = 200L, mapq_cutoff = 10, export = "BEADS", 
        rdata = FALSE, export_er = TRUE, quickMap = TRUE
    )
    
    
    final.path <- file.path('files', exp_dir, gsub('aligned\\^[A-Za-z0-9]+\\^[A-Za-z0-9]+', paste0('BEADSQ10UNIQ^linear^', genome), prefix))
    
    Entry <- addGenericFile(
        ids,
        path = final.path, 
        Processing = 'BEADSQ10UNIQ', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'B',
        genome = genome,
        comments = JADBtools::bamStats(basename(fls)),
        uniq = FALSE
    )
    
    out <- file.rename(basename(path(NRM0)), basename(Entry$path))
    
    
    message("Done!")
    
}

