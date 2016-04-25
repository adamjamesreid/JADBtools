Skipped: LIN13^Q2944_pk10^E^hpl2^L3_aligned^NA^NA_AA315^F3908442.bam



m <- vmatchPattern(reverseComplement(DNAString('TACBGTA')), Celegans, fixed = FALSE)
hits <- resize(m, 1000, fix='center')
cover <- coverage(hits)
export.bw(cover, 'TACBGTA_motif_coverage_1k_ce10.bw')

ss <- getSeq(Celegans, unl[unl  %over%  peak])

m <- vmatchPattern(reverseComplement(DNAString('TACBGTA')), Celegans, fixed = FALSE)
z <- vcountPattern(DNAString('TACBGTA'), ss, fixed = FALSE)

peak <- ChIPseeker::readPeakFile('http://ws190.gurdon.private.cam.ac.uk/REP/motif_peaks_replicates/HPL2%5eQ2324_final_peaks_0.05.narrowPeak', header=FALSE)

ssp <- getSeq(Celegans, peak)
z <- vcountPattern(DNAString('TACBGTA'), ssp, fixed = FALSE)


### GAATWWW or GAATWWWA - Pol2 peak motif


m <- vmatchPattern(DNAString('GAATWWW'), Celegans, fixed = FALSE)
hits <- resize(m, 1000, fix='center')

m2 <- vmatchPattern(DNAString('GAATWWWA'), Celegans, fixed = FALSE)

peaks <- sapply(dir('/Volumes/raid0/_Jcor', pattern = 'narrowPeak', full.names = TRUE), ChIPseeker::readPeakFile)
pp <- sapply(peaks, reduce)

zz <- Reduce(interaction, pp)

seqlengths(small) <- seqlengths(Celegans)[seqlevels(small)]
seqlengths(big) <- seqlengths(Celegans)[seqlevels(big)]

L = 250L
names(small) <- paste(small)
smfa <- getSeq(Celegans, trim(resize(small, L, fix = 'center')))
smfa <- smfa[width(smfa)==L]
#names(smfa) <- paste0('peak', 1:length(smfa))
writeXStringSet(smfa, 'pol2_larva_peaks_small_set.fa')
writeXStringSet(sample(smfa, 60), 'pol2_larva_peaks_small_subset_250bp.fa')

names(big) <- paste(big)
bgfa <- getSeq(Celegans, trim(resize(big, L, fix = 'center')))
bgfa <- bgfa[width(bgfa)==L]
writeXStringSet(bgfa, 'pol2_larva_peaks_big_set_250bp.fa')


parallel meme-chip {} -o memechip_{.} -meme-mod anr -meme-p 4 -time 60 \
-db ~/software/db/motif_databases/EUKARYOTE/jolma2013.meme \
-db ~/software/db/motif_databases/JASPAR/JASPAR_CORE_2014_vertebrates.meme \
-db ~/software/db/motif_databases/MOUSE/uniprobe_mouse.meme \
-db ~/software/db/motif_databases/WORM/uniprobe_worm.meme ::: *.fa


meme-chip pol2_larva_peaks_big_set.fa -o memechip_pol2_larva_peaks_big_set -db /mnt/home3/ahringer/ps562/software/db/motif_databases/EUKARYOTE/jolma2013.meme -db /mnt/home3/ahringer/ps562/software/db/motif_databases/JASPAR/JASPAR_CORE_2014_vertebrates.meme -db /mnt/home3/ahringer/ps562/software/db/motif_databases/MOUSE/uniprobe_mouse.meme -db /mnt/home3/ahringer/ps562/software/db/motif_databases/WORM/uniprobe_worm.meme 


meme-chip -db /mnt/home3/ahringer/ps562/software/db/motif_databases/JASPAR/JASPAR_CORE_2014_vertebrates.meme -o memechip_pol2_larva_peaks_big_set pol2_larva_peaks_big_set.fa