setwd('/scratch/asabe/projects/foundation-model/downstream/exon-trap')

library(data.table)
library(stringr)
library(magrittr)
library(Biostrings) %>% suppressPackageStartupMessages()
library(BSgenome.Hsapiens.UCSC.hg38) %>% suppressPackageStartupMessages()

set.seed(42)

ref_genome = BSgenome.Hsapiens.UCSC.hg38
# reporter_vectors_fasta_file = 'data/reporter_vectors.ready.fasta'
vector_nick_fasta_file = 'data/vector_nick_spliceai.fasta'
# query_exons_data_file = 'paper/data/data_zip/5A.txt'
query_exons_data_file = 'paper/data/data_zip/2D.txt'

processed_data_csv_file = '/home/asabe/scratch/projects/foundation-model/preprocess/pre-mrna/data/exon_trap.preprocessed.updated.csv.gz'
metadata_csv_file = '/home/asabe/scratch/projects/foundation-model/preprocess/pre-mrna/data/exon_trap.preprocessed.updated.metadata.csv.gz'

vector = readDNAStringSet(vector_nick_fasta_file)
vector = data.table(id = names(vector), len = width(vector), seq = as.character(vector))
vector = cbind(
  data.table(up_seq = vector[id == 'up_seq', seq]),
  data.table(down_seq = vector[id == 'down_seq', seq]))

query_flank_lens = 2 ** c(4, 6, 8, 9)
intron_start_seq = 'GTAAGTTATGAGGTAGAAAGGTCAACGTCTGA'
intron_end_seq = 'ATGCGTGTTGTGTTGCCTTTCTGTCTTCACAG'

# .listify = function(mat) list(mat[, 1], mat[, 2])

vector[, intron_start := str_locate(up_seq, intron_start_seq) %>% extract(, 'start')]
vector[, intron_end := str_locate(down_seq, intron_end_seq) %>% extract(, 'end')] 

seq_parts = c('exon_left', 'intron_left', 'intron_right', 'exon_right')

vector[, exon_left := str_sub(up_seq, 1, intron_start-1)]
vector[, intron_left := str_sub(up_seq, intron_start, -1)]
vector[, intron_right := str_sub(down_seq, 1, intron_end)]
vector[, exon_right := str_sub(down_seq, intron_end+1, -1)]

vector[, exon_left := str_c('HS', str_to_lower(exon_left))]
vector[, exon_right := str_c(str_to_lower(exon_right), 'E')]
vector[, left_seq := str_c(exon_left, intron_left)]
vector[, right_seq := str_c(intron_right, exon_right)]
vector = vector[, .(left_seq, right_seq)]
vector[, id := 'vecNick']
setnames(vector, 'id', 'vec_id')

exons = fread(query_exons_data_file)
if(!('exon_finder' %in% colnames(exons))) {
  exons[, exon_finder := 'ET']
}

cloned_vectors = 
  lapply(query_flank_lens, function(flank_len) {
    
    exons[, seq := getSeq(ref_genome,
                          names = chrom,
                          start = start - flank_len,
                          end = end + flank_len,
                          strand = strand,
                          as.character = T)]
    
    exons[, query_left := str_sub(seq, 1, flank_len)]
    exons[, query_right := str_sub(seq, nchar(seq)-flank_len+1, -1)]
    exons[, query_seq := str_sub(seq, flank_len+1, nchar(seq)-flank_len)]
    
    exons = exons[, .(exon_id, exon_finder, query_left, query_seq, query_right)]
    
    exon_vectors = exons[, as.list(vector), .(exon_id, exon_finder, query_left, query_seq, query_right)]
    exon_vectors[, query_intronic := str_c(left_seq, query_left, query_seq, query_right, right_seq)]
    exon_vectors[, query_exonic := str_c(left_seq, query_left, str_to_lower(query_seq), query_right, right_seq)]
    exon_vectors[, id_intronic := str_c(exon_id, vec_id, 'intron', flank_len, sep = '_')]
    exon_vectors[, id_exonic := str_c(exon_id, vec_id, 'exon', flank_len, sep = '_')]
    
    vectors_intronic = exon_vectors[, .(id_intronic, query_intronic)]
    vectors_exonic = exon_vectors[, .(id_exonic, query_exonic)]
    colnames(vectors_intronic) = colnames(vectors_exonic) = c('id', 'seq')
    
    vectors_all = rbind(vectors_intronic, vectors_exonic)
    vectors_all
  })

cloned_vectors = rbindlist(cloned_vectors)

setorder(cloned_vectors, id)

cloned_vectors[, seq_len := nchar(seq)]
setnames(cloned_vectors, 'id', 'exon_id')

fwrite(cloned_vectors, processed_data_csv_file)
fwrite(exons, metadata_csv_file)
