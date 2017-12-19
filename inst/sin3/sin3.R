tab <- jagui_get_table('labexperimentview')


tab %>% tbl_df %>% filter(Factor=='SIN3') -> sin3

sin3 %>% filter(ExtractID=='ME20') 
sin3 %>% filter(ExtractID=='aa44', Antibody=='Q5986') %>% .$ContactExpID -> g1


JADBtools::jadb_download_files(g1)
g1 %>% sapply(function(x) dir(pattern=x)) -> fls

out <- paste0('cat ', paste(fls, collapse = ' '), ' > SIN3^Q5986_aa44^E^N2^L3_raw^NA^NA_', paste(g1, collapse = '_'), '.txt.gz')
system(out)

jadb_dc('https://docs.google.com/spreadsheets/d/1PA_MRIcb4VM1P23WJQU1LjyikQ3V3rxUO3ORwOJK8G8/edit?usp=sharing', basespace = FALSE)


sin3 %>% filter(ExtractID=='aa44', Antibody=='Q6013') %>% .$ContactExpID  -> g2

JADBtools::jadb_download_files(g2)
g2 %>% sapply(function(x) dir(pattern=x)) -> fls

out <- paste0('cat ', paste(fls, collapse = ' '), ' > SIN3^Q6013_aa44^E^N2^L3_raw^NA^NA_', paste(g2, collapse = '_'), '.txt.gz')
system(out)



sin3$Antibody %>% unique

configr::read.config('~/.my.cnf')$jadb$host