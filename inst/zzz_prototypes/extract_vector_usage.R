require(devtools)
load_all()

fls <- dir()


fsl <- (
"http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/LET418/REP031/LET418_N2_L3_NORM_linear_1bp_IL024^IL025_F5b08121^F0308127.bw"                                       
[721] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/LIN35/REP034/LIN35_N2_L3_NORM_linear_1bp_MC018^AA125_Ffc07678^Fbb07696.bw"                                         
[722] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/LET418/REP030/LET418_N2_L3_NORM_linear_1bp_IL026^IL023_Fdd08103^F8008117.bw"                                       
[723] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/LIN61/REP029/LIN61_N2_L3_NORM_linear_1bp_AA012^AA160_F4407842^Fa507853.bw"                                         
[724] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/LIN13/REP027/LIN13_N2_L3_NORM_linear_1bp_AA207^PK013_Fa007861^F4707870.bw"                                         
[725] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/HPL2/REP022/HPL2_N2_youngAdult_NORM_linear_1bp_AA249^AA250_F9906934^Fc006939.bw"                                   
[726] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/H3K9me2/REP026/H3K9me2_N2_L3_NORM_linear_1bp_AA083^AA082_Fc107578^F2b07590.bw"                                     
[727] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/HPL2/REP024/HPL2_N2_L3_NORM_linear_1bp_AA382^AA381_Fad09213^F6009212.bw"                                           
[728] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/LIN61/REP028/LIN61_N2_L3_NORM_linear_1bp_AA005^AA156_F5a07852^F9107854.bw"                                         
[730] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/MET2/REP025/MET2_N2_L3_NORM_linear_1bp_AA279^AA297_F5a07741^F5107932.bw"                                           
[731] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/LIN35/REP033/LIN35_N2_L3_NORM_linear_1bp_IL016^IL045_Fab07305^Fa707311.bw"                                         
[732] "http://jadb.gurdon.private.cam.ac.uk/db4/files/REPLICATES/HPL2/REP023/HPL2_N2_L3_NORM_linear_1bp_PK011^PK010_F1907331^F6c07325.bw" 

data <- lapply(fls, function(f) {
    message(f)
    try({
        bwf <- BigWigFile(f)
        vec <- extarct_vector(bwf, which = unlist(repeatModel), size = 1L)
        return(vec)
    })
})