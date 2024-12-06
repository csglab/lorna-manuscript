# setwd('/scratch/asabe/csg/mpaqt/lornash')

library(data.table)
library(stringr)
library(magrittr)

### MPAQT LR+SR
mpaqt_files = list.files('MPAQT/projects/human.c2l2.revio.gencode_v47.GRCh38_p14/samples', '*.MPAQT.LR_SR.tsv', recursive = T, full.names = T)
mpaqt_ids = mpaqt_files %>% basename() %>% str_remove('.MPAQT.LR_SR.tsv')

mpaqt = lapply(mpaqt_files, fread)
names(mpaqt) = mpaqt_ids

mp = lapply(names(mpaqt), function(id) {
  .temp = mpaqt[[id]]
  .temp = .temp[, .(transcript_id, TPM)]
  colnames(.temp) = c('transcript_id', id)
  .temp
})

mps = Reduce(merge, mp)

# cell_lines = mps[, -'transcript_id'] %>% colnames() %>% str_remove('_1$|_2$') %>% unique()
# for(cl in cell_lines) {
#   mps[, (cl) := rowMeans(.SD), .SDcols = patterns(cl)]
# }
# mps = mps[, .SD, .SDcols = c('transcript_id', cell_lines)]

### MPAQT SR
mpaqt_sr_files = list.files('MPAQT/projects/human.c2l2.revio.gencode_v47.GRCh38_p14/samples', '*.MPAQT.SR.tsv', recursive = T, full.names = T)
mpaqt_ids = mpaqt_sr_files %>% basename() %>% str_remove('.MPAQT.SR.tsv')

mpaqt_sr = lapply(mpaqt_sr_files, fread)
names(mpaqt_sr) = mpaqt_ids

mp_sr = lapply(names(mpaqt_sr), function(id) {
  .temp = mpaqt_sr[[id]]
  .temp = .temp[, .(transcript_id, TPM)]
  colnames(.temp) = c('transcript_id', id)
  .temp
})

mps_sr = Reduce(merge, mp_sr)

fwrite(mps, 'MPAQT/projects/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.mpaqt.tpm.csv')

fwrite(mps_sr, 'MPAQT/projects/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.short_reads.tpm.csv')
