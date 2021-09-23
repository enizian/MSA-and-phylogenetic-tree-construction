library(tidyverse)
headers = readLines("blastout.mega.WS11.tsv.gz", n = 4)
headers = headers[4]
headers = headers %>%
  str_remove("# Fields: ") %>%
  str_split(", ") %>% 
  unlist() %>%
  make.names() %>%
  str_replace(fixed(".."), ".") %>%
  str_replace("X.identity", "pct.identity")
megaWS11 = read_tsv("blastout.mega.WS11.tsv.gz", col_names=headers, comment="#")
megaWS28 = read_tsv("blastout.mega.WS28.tsv.gz", col_names=headers, comment="#")
blastnWS11 = read_tsv("blastout.task_blastn.WS11.tsv.gz", col_names=headers, comment="#")
dc_megaWS11 = read_tsv("blastout.task_dc-megablast.WS11.tsv.gz", col_names=headers, comment="#")
tblastx = read_tsv("blastout.tblastx.tsv.gz", col_names=headers, comment="#")

blast.results = bind_rows(list(megaWS11=megaWS11,
                                megaWS28=megaWS28, 
                                blastnWS11=blastnWS11, 
                                dc_megaWS11=dc_megaWS11,
                                tblastx=tblastx), 
                           .id="strategy")

uniq.blast.results = blast.results %>%
  group_by(strategy, subject.acc.ver) %>%
  filter(rank(dplyr::desc(alignment.length), ties.method = "first")==1)

library(UpSetR)
upset.table = uniq.blast.results %>% select(subject.acc.ver, strategy) %>% 
  table() %>%
  as.data.frame.matrix()
upset(upset.table)

uniq.blastn = uniq.blast.results %>%
    ungroup() %>%
    filter(strategy=="blastnWS11",
           str_detect(subject.title, "complete genome")) %>% 
    separate(subject.title,
             into=c("acc", "isolate", "complete", "name", "country", "host"),
             remove = FALSE,
             sep="\\|")

filtered.blastn = uniq.blastn[-c(1462, 1458, 1585),]
filtered.blastn <- filtered.blastn %>%
  mutate_all(function(x) ifelse(x=="", NA, x))
filtered.blastn <- na.omit(filtered.blastn)
filtered.blastn = filtered.blastn %>% 
  group_by(name, country, host) %>% 
  filter(rank(dplyr::desc(alignment.length), ties.method = "first")==1) %>% 
  filter(alignment.length >= 5000)

system("time mafft --maxiterate 100 --thread 2 --reorder --op 0.5 selected_viral_seqs_195.fa  > mafft_maxiter100_195_op.5.fa")
             
library(Biostrings)
ncbi.seqs <- readDNAStringSet("ncbi_virus_110119_2.txt")
selected.seqs = ncbi.seqs[names(ncbi.seqs) %in% filtered.blastn$subject.title]
seq.h = readDNAStringSet("seq_h.txt")
selected.seqs = c(selected.seqs, seq.h)
writeXStringSet(selected.seqs, "selected_viral_seqs.fa")

filtered.blastn = filtered.blastn %>%
  group_by(name, country, host) %>%
  filter(rank(dplyr::desc(alignment.length),ties.method = "min") < 10) %>%
  filter(rank(dplyr::desc(pct.identity), ties.method ="min") <  3)

inpath = "mafft_maxiter100_195_op.5.fa"
outpath = "mafft_maxiter100_195_op.5_trimmed_75pct.fa"
alignment = readDNAMultipleAlignment(inpath)
alignment = DNAMultipleAlignment(alignment,start=1000,end=48250)
alignment = maskGaps(alignment, min.fraction=0.25, min.block.width=1)
alignment = alignment %>% as("DNAStringSet") 
newnames = names(alignment) %>% 
  tibble(name=.) %>%
  mutate(name=str_replace_all(name," ","_")) %>%
  separate(name, 
           into=c("acc", "isolate", "complete", "name", "country", "host"),
           sep="_?\\|_?") %>%
  mutate(name=str_replace(name, "Middle_East_respiratory_syndrome-related","MERS"),
         name=str_replace(name, "Severe_acute_respiratory_syndrome-related", "SARS"),
         newname=paste(acc,name,country,host,sep="|")) %>%
  pull(newname)
names(alignment) = newnames
alignment %>% writeXStringSet(outpath)

system("FastTree -nt -gtr -gamma -out mafft_maxiter100_195_op.5_trimmed.fasttree.tre mafft_maxiter100_195_op.5_trimmed_75pct.fa")
