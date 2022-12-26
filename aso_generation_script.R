

setwd("E:/2022_SysPharm_HSY/00_Projects/2022_20 W1 Reciprocal Regulation")
source('./ASO_Generator.R')
sessionInfo()

gencode.42.gtf <- rtracklayer::import('./gencode.v42.annotation.gtf.gz')
wfdc1.info.all <- data.frame(gencode.42.gtf[gencode.42.gtf$gene_name == 'WFDC1',])


wfdc1.related.fasta <- Biostrings::readDNAStringSet('./Homo_sapiens_WFDC1_sequence.fa')
wfdc1.genomic <- wfdc1.related.fasta[55]
flank.start <- str_split(string = names(wfdc1.genomic), pattern = ':')[[1]][5] %>% as.numeric()
flank.end <- str_split(string = names(wfdc1.genomic), pattern = ':')[[1]][6] %>% as.numeric()

gene.start <- 84294846
gene.end <- 84329844

w1.gene.seq <- wfdc1.genomic
names(w1.gene.seq) <- 'hWFDC1'
w1.gene.seq <- w1.gene.seq$hWFDC1[(gene.start - flank.start +1):(gene.end-flank.start+1)]
wfdc1.genomic.character <- as.character(w1.gene.seq)



# prepare BLAST DB

bl_refseq <- blast(db="E:/2022_SysPharm_HSY/13_Scripts/BLAST/DB/refseq_select_rna/refseq_select_rna")
bl_hg <- blast(db="E:/2022_SysPharm_HSY/13_Scripts/BLAST/DB/human_genome/GCF_000001405.39_top_level")


human.fa.all <- readDNAStringSet('./Orthologues/Homo_sapiens_WFDC1_sequence.fa')
monkey.fa.all <- readDNAStringSet('./Orthologues/Macaca_fascicularis_WFDC1_sequence.fa')
rat.fa.all <- readDNAStringSet('./Orthologues/Rattus_norvegicus_Wfdc1_sequence.fa')
mouse.fa.all <- readDNAStringSet('./Orthologues/Mus_musculus_Wfdc1_sequence.fa')
dog.fa.all <- readDNAStringSet('./Orthologues/Canis_lupus_familiarisboxer_WFDC1_sequence.fa') # boxer 

seq.files.list <- list(human.fa.all, monkey.fa.all, rat.fa.all, mouse.fa.all)
seq.files.all <- DNAStringSet(c(human.fa.all, monkey.fa.all, rat.fa.all,mouse.fa.all,dog.fa.all))


aso.19 <- ASOGeneratorFunction01(wfdc1.info.all, 'WFDC1-201', 19, wfdc1.genomic.character)
aso.19.blast <- ASOGeneratorFunction02(aso.19, genome.blast = TRUE)
aso.19 <- aso.19.blast$aso.df
aso.19 <- ASOGeneratorFunction03(aso.19)

aso.19.filter <- aso.19 %>% filter( (genome_hits < 2)|  (refseq_hit ==1 & (is.na(genome_hits)|genome_hits <2)))  %>% 
    #filter(GC.Contents > 30 & GC.Contents < 70) %>% 
    filter(mouse.perfect.match == 'match'| rat.perfect.match =='match') 


aso.20 <- ASOGeneratorFunction01(wfdc1.info.all, 'WFDC1-201', 20, wfdc1.genomic.character)
aso.20.blast <- ASOGeneratorFunction02(aso.20)
aso.20 <- aso.20.blast$aso.df 
aso.20 <- ASOGeneratorFunction03(aso.20)

aso.20.filter <- aso.20 %>% filter( (genome_hits < 2)|  (refseq_hit ==1 & (is.na(genome_hits)|genome_hits <2)))  %>% 
    #filter(GC.Contents > 30 & GC.Contents < 70) %>% 
    filter(mouse.perfect.match == 'match'| rat.perfect.match =='match') 


aso.18 <- ASOGeneratorFunction01(wfdc1.info.all, 'WFDC1-201', 18, wfdc1.genomic.character)
aso.18.blast <- ASOGeneratorFunction02(aso.18)
aso.18 <- aso.18.blast$aso.df 
aso.18 <- ASOGeneratorFunction03(aso.18)
aso.18.filter <- aso.18 %>% filter( (genome_hits < 2)|  (refseq_hit ==1 & (is.na(genome_hits)|genome_hits <2)))  %>% 
    #filter(GC.Contents > 30 & GC.Contents < 70) %>% 
    filter(mouse.perfect.match == 'match'| rat.perfect.match =='match')


# Filtering Step for All range of ASO
aso.filters <- list(aso.18.filter,aso.19.filter,aso.20.filter)
aso.filters <- rbindlist(aso.filters,use.names = TRUE) %>% as_tibble()

aso.filters <- aso.filters %>% 
    pivot_longer(cols = c('mouse.perfect.match','dog.perfect.match','monkey.perfect.match','rat.perfect.match')) %>% 
    mutate(matched = ifelse(value == 'match',1,0)) %>% 
    group_by(aso_idx,width) %>% 
    mutate(match_count = sum(matched)) %>% 
    ungroup() %>% dplyr::select(-matched) %>%
    pivot_wider(names_from = 'name', values_from = 'value')

aso.filters.final <- aso.filters %>% filter(mouse.perfect.match == 'match') %>% 
    arrange(chr_start,desc(match_count), desc(width)) %>% 
    distinct(chr_start, .keep_all = TRUE)  

write_clip(aso.filters.final)

xlsx::write.xlsx(aso.18, file = "W1 ASO 18-20mer Design.xlsx", row.names = FALSE)

write_csv(aso.18, "W1 ASO 18-20mer Design.csv", progress = TRUE)
write_csv(aso.19, "W1 ASO 19mer Design.csv", progress = TRUE)
write_csv(aso.20, "W1 ASO 20mer Design.csv", progress = TRUE)

#### 17mer ####


aso.17 <- ASOGeneratorFunction01(wfdc1.info.all, 'WFDC1-201', 17, wfdc1.genomic.character)
aso.17.blast <- ASOGeneratorFunction02(aso.17)
aso.17 <- aso.17.blast$aso.df 
aso.17 <- ASOGeneratorFunction03(aso.17)
aso.17.filter <- aso.17 %>% filter( (genome_hits < 2)|  (refseq_hit ==1 & (is.na(genome_hits)|genome_hits <2)))  %>% 
    #filter(GC.Contents > 30 & GC.Contents < 70) %>% 
    filter(mouse.perfect.match == 'match'| rat.perfect.match =='match')

aso.17.filter


