library(data.table)
library(stringr)
library(Biostrings)

project = 'encode4'
setwd(str_glue('/scratch/asabe/projects/{project}'))

cell_lines_transcriptome_fasta_files = list.files('data/long-reads/batches-all/encode-pipeline/bulk/transcriptome-fasta', pattern = '*.fasta', full.names = T)

temp = str_split(cell_lines_transcriptome_fasta_files, '/', simplify = T)
cell_lines = temp[, ncol(temp)] %>% str_remove('.fasta')

names(cell_lines_transcriptome_fasta_files) = cell_lines

cl_trs = lapply(cell_lines_transcriptome_fasta_files, readDNAStringSet)

reference_transcriptome_fasta_file = 'data/references/transcriptome/gencode.v43.transcripts.fa'

ref_tr = readDNAStringSet(reference_transcriptome_fasta_file)
names(ref_tr) = str_split_fixed(names(ref_tr), '[|]', n = 2)[, 1]

cl_trs_ref = lapply(cl_trs, function(cl_tr) {
  cl_tr_ref = c(cl_tr, ref_tr)
  cl_tr_ref = cl_tr_ref[!duplicated(names(cl_tr_ref))]
  cl_tr_ref
})

lapply(names(cl_trs_ref), function(cell_line) {
  writeXStringSet(cl_trs_ref[[cell_line]],
                  str_glue('data/long-reads/batches-all/encode-pipeline/bulk/transcriptome-fasta/{cell_line}.gencode_v43.fasta'),
                  format = 'fasta')
})
