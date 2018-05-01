blacklist <- import.bed('https://gist.githubusercontent.com/Przemol/ef62ac7ed41d3a84ad6c478132417770/raw/56e98b99e6188c8fb3dfb806ff6f382fe91c27fb/CombinedBlacklists.bed')
nonmappable <- import.bed('https://gist.githubusercontent.com/Przemol/ef62ac7ed41d3a84ad6c478132417770/raw/56e98b99e6188c8fb3dfb806ff6f382fe91c27fb/non_mappable.bed')

if(filter_map) {
    filter <- GenomicRanges::reduce(c(blacklist, nonmappable))
} else {
    filter <- blacklist
}

download.file(
    'https://gist.githubusercontent.com/Przemol/e9f1a3a5053619e69fafbd46759a17e4/raw/c69209270843e2724dc8c6ea9715ddd52930f41c/ce10ToCe11.over.chain', 
    (tempfile() -> chainf), quiet = TRUE
)
chain_ce10ToCe11 <- import.chain(chainf)
blacklist_ce11 <- liftOver(blacklist, chain_ce10ToCe11) %>% GenomicRanges::reduce(min.gapwidth=50) %>%  unlist
filter_ce11 <- liftOver(filter, chain_ce10ToCe11) %>% GenomicRanges::reduce(min.gapwidth=50) %>%  unlist
