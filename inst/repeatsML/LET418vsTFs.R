files <- dir('/Volumes/data/_HOT_jointAWG_worm_fly_human/worm_Ce10_noDoubleCounting/source_data')
stage <- strsplit(dir('/Volumes/data/_HOT_jointAWG_worm_fly_human/worm_Ce10_noDoubleCounting/source_data'),'_')  %>% sapply('[[', 3)
factor <- strsplit(dir('/Volumes/data/_HOT_jointAWG_worm_fly_human/worm_Ce10_noDoubleCounting/source_data'),'_')  %>% sapply('[[', 2)

reg <- ChIPseeker::readPeakFile(files[1], header=FALSE)