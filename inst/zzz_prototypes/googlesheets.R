require(xlsx)
require(writexl)
load("/Users/przemol/Downloads/dds_ja_fp_cfp1.Rdata")


xlsx <- write_xlsx(all_res, path = "cfp1_sin3_set2_de_res_star_WBcel235e90.xlsx")
gfile <- drive_upload(
    xlsx, 
    'DE/cfp1_sin3_set2_de_res_star_WBcel235e90_auto', 
    type='spreadsheet'
)
drive_browse(gfile)
drive_share

library("googledrive")
require(tidyverse)

devtools::install_github("ropensci/sheetseeR")

install.packages("googlesheets")

require('googlesheets')
require(magrittr)
options(httr_oob_default=TRUE)
gs_auth()



gs_delete()


gs_edit_cells(zz, ws = 1, input = head(all_res[[1]]))

names(all_res)
x=names(all_res)[[1]]

gs <- gs_new('cfp1_sin3_set2_de_res_star_WBcel235e90')
gs %<>% gs_ws_new(ws_title = x)
for(x innames(all_res)) {
    message('==>> ', x)
    gs %<>% gs_ws_new(ws_title = x)
    gs %<>% gs_edit_cells(ws = x, input = all_res[[x]])
}


gs_browse(zz)
gs_delete(zz)