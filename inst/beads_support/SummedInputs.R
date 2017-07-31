sumBAMinputs(
    bam.controls = "ce11_FRM_HiSeq_input.bam", 
    bw.mappability = 'mappability/ce11_gem-mappability_36bp.bw', 
    genome = 'ce11',
    out_name = "ce11_FRM_HiSeq_nonUNIQ_200bp",
    uniq = FALSE, 
    insert = 200L, 
    mapq_cutoff = 10L, 
    quickMap = TRUE,
    bin = 25L
)

sumBAMinputs(
    bam.controls = "ce11_EGS_HiSeq_input.bam", 
    bw.mappability = 'mappability/ce11_gem-mappability_36bp.bw', 
    genome = 'ce11',
    out_name = "ce11_EGS_HiSeq_nonUNIQ_MAPQ10_200bp",
    uniq = FALSE, 
    mapq_cutoff = 10L, 
    insert = 200L, 
    quickMap = TRUE,
    bin = 25L
)



sumBAMinputs(
    bam.controls = "ce11_FRM_HiSeq_input.bam", 
    bw.mappability = '../../../_mappability_files_/ce11_gem-mappability_36bp.bw', 
    genome = 'ce11',
    out_name = "ce11_FRM_HiSeq_nonUNIQ_nonMAPQ_200bp",
    uniq = FALSE, 
    insert = 200L, 
    mapq_cutoff = 0, 
    quickMap = TRUE,
    bin = 25L
)

sumBAMinputs(
    bam.controls = "ce11_EGS_HiSeq_input.bam", 
    bw.mappability = '../../../_mappability_files_/ce11_gem-mappability_36bp.bw', 
    genome = 'ce11',
    out_name = "ce11_EGS_HiSeq_nonUNIQ_nonMAPQ_200bp",
    uniq = FALSE, 
    mapq_cutoff = 0, 
    insert = 200L, 
    quickMap = TRUE,
    bin = 25L
)



## ce10 


#ce10_FRM_HiSeq_nonUNIQ_MAPQ10_200bp_SummedInput_bin25bp.bw
sumBAMinputs(
    bam.controls = "ce10_FRM_HiSeq_input.bam", 
    bw.mappability = '../../../_mappability_files_/ce10_gem-mappability_36bp.bw', 
    genome = 'ce10',
    out_name = "ce10_FRM_HiSeq_nonUNIQ_MAPQ10_200bp",
    uniq = FALSE, 
    insert = 200L, 
    mapq_cutoff = 10L, 
    quickMap = TRUE,
    bin = 25L
)
file.rename('ce10_FRM_HiSeq_nonUNIQ_MAPQ10_200bp_SummedInput_linear_25bp.bw', 'ce10_FRM_HiSeq_nonUNIQ_MAPQ10_200bp_SummedInput_bin25bp.bw')


sumBAMinputs(
    bam.controls = "ce10_EGS_HiSeq_input.bam", 
    bw.mappability = '../../../_mappability_files_/ce10_gem-mappability_36bp.bw', 
    genome = 'ce10',
    out_name = "ce10_EGS_HiSeq_nonUNIQ_MAPQ10_200bp",
    uniq = FALSE, 
    insert = 200L, 
    mapq_cutoff = 10L, 
    quickMap = TRUE,
    bin = 25L
)
file.rename('ce10_EGS_HiSeq_nonUNIQ_MAPQ10_200bp_SummedInput_linear_25bp.bw', 'ce10_EGS_HiSeq_nonUNIQ_MAPQ10_200bp_SummedInput_bin25bp.bw')


sumBAMinputs(
    bam.controls = "ce10_FRM_HiSeq_input.bam", 
    bw.mappability = '../../../_mappability_files_/ce10_gem-mappability_36bp.bw', 
    genome = 'ce10',
    out_name = "ce10_FRM_HiSeq_nonUNIQ_nonMAPQ_200bp",
    uniq = FALSE, 
    insert = 200L, 
    mapq_cutoff = 0L, 
    quickMap = TRUE,
    bin = 25L
)
file.rename('ce10_FRM_HiSeq_nonUNIQ_nonMAPQ_200bp_SummedInput_linear_25bp.bw', 'ce10_FRM_HiSeq_nonUNIQ_nonMAPQ_200bp_SummedInput_bin25bp.bw')


sumBAMinputs(
    bam.controls = "ce10_EGS_HiSeq_input.bam", 
    bw.mappability = '../../../_mappability_files_/ce10_gem-mappability_36bp.bw', 
    genome = 'ce10',
    out_name = "ce10_EGS_HiSeq_nonUNIQ_nonMAPQ_200bp",
    uniq = FALSE, 
    insert = 200L, 
    mapq_cutoff = 0L, 
    quickMap = TRUE,
    bin = 25L
)
file.rename('ce10_EGS_HiSeq_nonUNIQ_nonMAPQ_200bp_SummedInput_linear_25bp.bw', 'ce10_EGS_HiSeq_nonUNIQ_nonMAPQ_200bp_SummedInput_bin25bp.bw')
