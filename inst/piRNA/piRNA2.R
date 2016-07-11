
'/Volumes/ahringer/Przemyslaw Stempor/PrzemeksAnnotationsCollection/indi_piRNAs_ce10.bed'  %>% import.bed() -> indi

lst <- GRangesList(list(all, dep, indi))
names(lst) <- c("all", "dep", "indi")



clust <- GRanges(rep('chrIV', 2), IRanges(c(4770508, 13259587), c(7108663, 17493793)))
lst   %>% sapply(function(x) sum(x %over% clust))
lst   %>% sapply(function(x) sum(!x %over% clust)) 

lst   %>% sapply(function(x) sum(x %over% clust)) / lengths(lst) -> z


require(Biostrings)
#GSM1195513	N2 5' independent
independent <- readDNAStringSet('http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1195513&format=file&file=GSM1195513%5Fcol%5FtrimcutSX1316%5F5I%2Efasta%2Egz')

#GSM1195514	N2 5' dependent
#GSM1195515	N2 TAP
TAP_EM <-  readDNAStringSet('http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1195515&format=file&file=GSM1195515%5Fcol%5FtrimcutSX1316PP%5FS2%5FL001%5FR1%5F001%2Efasta%2Egz')

#GSM1195516	N2 polyphosphatase


# GSM1195520	prde-1_1 5'independent
# GSM1195521	prde-1_1 5'dependent
# GSM1195522	prde-1 TAP
# GSM1195523	prde-1_2 5' independent
# GSM1195524	prde-1_2 5' dependent


#GSM984436	CIP-TAP, insert 18-40 nt adult

TAP_CM <-  readDNAStringSet('http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM984436&format=file&file=GSM984436%5FCIPTAP0118%5Fprocessed%2Efa%2Egz')

length(TAP_EM)
length(TAP_CM)
sum(TAP_EM  %in% TAP_CM)

#Top1%
TAP_EM_top <- head(TAP_EM, round(length(TAP_EM)*.001))
TAP_CM_top <- head(TAP_CM, round(length(TAP_CM)*.001))

pdict1 <- PDict(TAP_EM_top, tb.start=2, tb.end=18)
TAP_EM_top_locations <- unlist(GRangesList(lapply(names(Celegans), function(chr) {
    message(chr)
    loc <- matchPDict(pdict1, Celegans[[chr]])
    GRanges(chr, unlist(loc))
})))

pdict1 <- PDict(TAP_CM_top, tb.start=1, tb.end=17)
TAP_CM_top_locations <- unlist(GRangesList(lapply(names(Celegans), function(chr) {
    message(chr)
    loc <- matchPDict(pdict1, Celegans[[chr]])
    GRanges(chr, unlist(loc))
})))

require(ChIPseeker)
vennplot(TAP_EM_top_locations, TAP_CM_top_locations)

intersect(TAP_CM_top_locations, TAP_EM_top_locations)
table(TAP_EM_top_locations  %over% TAP_CM_top_locations)
table(TAP_CM_top_locations  %over% TAP_EM_top_locations)


require(BSgenome.Celegans.UCSC.ce10)


pwm <- PWM(DNAStringSet(head(TAP_EM_top, 10), end=21), type='prob') # Computes a PWM for DNA fragments.
library(seqLogo); seqLogo(t(t(pwm) * 1/colSums(pwm)))

pwm <- PWM(DNAStringSet(head(TAP_CM_top, 10), end=17)) # Computes a PWM for DNA fragments.
library(seqLogo); seqLogo(t(t(pwm) * 1/colSums(pwm)))

TAP_CM_top_locations  %>% export.bed('TAP_CM_top_locations001.bed')
TAP_EM_top_locations %>% export.bed('TAP_EM_top_locations001.bed')

TAP_CM_top_locations  %>% export.bed('TAP_CM_top_locations_top01.bed')
TAP_EM_top_locations %>% export.bed('TAP_EM_top_locations_top01.bed')

int <- intersect(TAP_CM_top_locations, TAP_CM_top_locations)
int %>% export.bed('Intersect_top001.bed')
