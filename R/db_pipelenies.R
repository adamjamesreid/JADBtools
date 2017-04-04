jadb_ChIPseq <- function(ids) {
    
    setwd('/mnt/jadb/DBfile/DBfiles')
    library('rbeads')
    
    message('\t => \tPerforming alignment')
    jadb_addAlignedBAM(ids)
    message('\t => \tExporting tracks')
    jadb_addTracksFromBAM(ids)
    message('\t => \tNormalising NU')
    addNonUniqueQ10Beads(ids)
    
}