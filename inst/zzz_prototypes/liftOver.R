t1 <- import.bed(pp[[1]], extraCols=extraCols_narrowPeak)
liftOver(t1, import.chain(chain)) %>% reduce(min.gapwidth=10) %>% lengths() %>% max
liftOver(t1, import.chain(chain)) %>% reduce(min.gapwidth=10) %>% unlist -> lo
elementMetadata(lo) <- elementMetadata(t1)
export.bed(lo, "SIN3_Q5986_E_N2_L3_BEADSQ10NU_linear_XM001^AA404_IDRconcave_ce10.narrowPeak")

t2 <- import.bed(pp[[2]], extraCols=extraCols_narrowPeak)
liftOver(t2, import.chain(chain)) %>% reduce(min.gapwidth=10) %>% lengths() %>% max
liftOver(t2, import.chain(chain)) %>% reduce(min.gapwidth=10) %>% unlist -> lo
elementMetadata(lo) <- elementMetadata(t2)
export.bed(lo, "SIN3_Q6013_E_N2_L3_BEADSQ10NU_linear_XM002^AA405_IDRconcave_ce10.narrowPeak")