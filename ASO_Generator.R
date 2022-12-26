### 
#
#  ASO Generator Functions 
# 2022-12-16 HSY 
# 
# 
# ASO generation via...
# ASOGeneratorFunction01- initial enumeration 
# ASOGeneratorFunction02 - BLAST
# ASOGeneratorFunction03 - Ortholog Search

# ASOGeneratorFunction01 - enumerate ASO and add basic informations 
# Inputs : 
# gene.gtf.df : gene_information from Gencode GTF, 
# main_transcripts : a character vector of transcript_name of main transcript to make positional information 
# gene.seq.character, gene.start.position

# aso width,  ( to indicate genomic position of ASO)
# Outputs : ASO design with chromosome position, antisense, sense sequence, 

library(rtracklayer)
library(tidyverse)
library(glue)
library(progress)
library(ggtranscript)
library(clipr)
library(data.table)


library(GenomicRanges)
library(Biostrings)
library(BSgenome)
library(rBLAST)
library(org.Hs.eg.db)


source('./ASO_Generator_Help.R')


ASOGeneratorFunction01 <- function(gene.gtf.df, main_transcript, width, gene.seq.character){
    
    # Enumerate ASO with fixed width 
    # Add sequence and antisense information 
    # Add Position & Sequence Information 
    
    # Gene information process

    
    gene.gtf.df.gene <- gene.gtf.df %>% dplyr::filter(type =='gene')
    gene.gtf.df.transcript <- gene.gtf.df %>% dplyr::filter(type == 'transcript') %>% arrange(transcript_name)
    
    gene.gtf.df.main.transcript <- gene.gtf.df %>% dplyr::filter(transcript_name == main_transcript & type !='transcript') %>% 
        dplyr::select(start,end, width, type, exon_number) %>% 
        mutate(exon_idx = paste0('exon_', exon_number)) %>% mutate(subtype = ifelse(type=='exon',exon_idx, as.character(type)))
    
    gene.gtf.df.main.transcript.intron <- ggtranscript::to_intron(gene.gtf.df.main.transcript %>% filter(type == 'exon')) %>% 
        mutate(start = start +1, end = end-1, exon_number = as.numeric(exon_number)-1) %>% 
        mutate(width = end-start + 1 , subtype = paste0(type,'_',exon_number), 
               exon_idx = paste0('exon_',exon_number), type = as.character(type), exon_number = as.character(exon_number)) %>% 
        dplyr::select(start, end, width, type, exon_number, exon_idx, subtype)
    
    gene.gtf.df.main.transcript <- gene.gtf.df.main.transcript %>% mutate(type=as.character(type)) %>% rows_append(gene.gtf.df.main.transcript.intron)
    
    gene.genomic.start <- gene.gtf.df.gene$start

    # make aso df 
    

    aso.df <- data.frame(aso_idx = 1:(gene.gtf.df.gene$width-width+1),
                            chr_start = gene.gtf.df.gene$start:(gene.gtf.df.gene$end-width + 1) ,
                            chr_end = (gene.gtf.df.gene$start+width - 1) :(gene.gtf.df.gene$end), 
                            width = width ) 
    aso.df$seq <- sapply(aso.df$chr_start, add.seq, seq=gene.seq.character, seq.start.idx = gene.genomic.start, aso.width = width)

    # add position information
    aso.df$position <- ''
    print('add position information')
    pb <- progress_bar$new(format = "  Processing [:bar] :current/:total :percent eta: :eta",
                           total = nrow(aso.df), 
                           clear = FALSE)
    for(i in 1:nrow(aso.df)){
        pb$tick()
        position.names <- asoModuleFinder.onlynames(aso.df$chr_start[i], aso_end = aso.df$chr_end[i],
                                                 gene.gtf.df.main.transcript$start, 
                                                 gene.gtf.df.main.transcript$end, 
                                                 gene.gtf.df.main.transcript$subtype)
        aso.df$position[i] <- paste0(unique(position.names[position.names !='']), collapse = '_')
    }
    
    
    # make junctional aso ( should change add.seq )
    
    aso.df.junction <- createJunctionalASO(gene.gtf.df, 
                                           gene.gtf.df.transcript %>% dplyr::filter(transcript_name == main_transcript), 
                                           aso.width = width, 
                                           gene.seq.character, 
                                           gene.genomic.start)
    
    
    
    # merge ASO and Junctional DF 
    
    aso.df.junction <- aso.df.junction %>% mutate(chr_start = j1.start, chr_end = j2.end)
    aso.df.cols <- c("aso_idx","width","position","seq","chr_start","chr_end") 
    
    aso.df$aso_idx <- as.character(aso.df$aso_idx)
    
    aso.df <- aso.df %>% dplyr::select(aso.df.cols) %>%  bind_rows(aso.df.junction %>% dplyr::select(aso.df.cols))
    
    
    # Calculate GC contents
    
    aso.df$GC.Contents <- sapply( aso.df$seq, x <- function(x){round( str_count(x, pattern = "G|C") / width * 100 ) })
    
    # Calculate Duplicated  genes 
    
    duplicated.seqs <- aso.df %>% group_by(seq) %>% summarise(duplicated.in.gene = n()) %>% arrange(desc(duplicated.in.gene))
    aso.df <- aso.df %>% left_join(duplicated.seqs)
    
    # Calculate Antisense sequence 
    aso.df$antisense <- as.character(reverseComplement(DNAStringSet(aso.df$seq)))
    
    aso.df <- aso.df %>% dplyr::select(aso_idx, width, position, antisense, seq, chr_start, chr_end, everything())
    return(aso.df)
    
    
}


# 2. BLAST 

ASOGeneratorFunction02 <- function(aso.df, blast.db, Perc.Ident = 100, batch = 100, genome.blast = TRUE){
    
    # Function call for BLAST 
    
    seq <- aso.df$seq
    names(seq) <- aso.df$aso_idx
    seq <- DNAStringSet(x = seq, use.names = TRUE)
    width <- aso.df$width[1]
    
    blast.result <- list()
    
    # Refseq BLAST 
    print('refSeq BLAST')
    refseq.blast <- predict(object = bl_refseq, 
                            newdata = seq, 
                            BLAST_args = glue('-taxids 9606 \\
                                              -task megablast \\
                                              -word_size {width}'))
    
    refseq.blast.perfect <- refseq.blast %>% 
        filter(Perc.Ident == Perc.Ident & Alignment.Length == width & Q.end == width &  Q.start == 1) 
    
    
    refseq.blast.perfect$Symbol <- mapIds(x =org.Hs.eg.db, 
                                          keys = gsub(refseq.blast.perfect$SubjectID, pattern = '\\..*$', replacement = ''),
                                          keytype = 'REFSEQ', 
                                          column = 'SYMBOL', 
                                          multiVals = function(x){paste0(x, collapse = ',')}) 
    
    refseq.hit.df <- refseq.blast.perfect %>% 
        dplyr::select(QueryID, Symbol) %>% 
        group_by(QueryID) %>% 
        mutate(refseq_hit=n()) %>% ungroup() %>% 
        dplyr::select(QueryID,refseq_hit) %>% 
        mutate(QueryID = as.character(QueryID)) %>% 
        distinct()
    
    aso.df <- aso.df %>% left_join(refseq.hit.df, by = c('aso_idx' = 'QueryID'))
    
    blast.result[[1]] <- refseq.blast
    
    # genome blast 
    if(genome.blast){
        
        aso.genome.blast <- data.frame()
        b <- batch
        
        pb <- progress_bar$new(format = "  Processing [:bar] :current/:total :percent eta: :eta",
                               total = ceiling(length(seq)/b), 
                               clear = FALSE)
        
        for(i in 1:(ceiling(length(seq)/b)) ){
            
            pb$tick()
            # 1~10 11~20 21 ~ 25
            
            t.end <- min( b*i, length(seq))
            
            gb.i <- predict(object = bl_hg, 
                            newdata = seq[( b*(i-1) + 1 ) : t.end], 
                            BLAST_args = glue('-taxids 9606 \\
                                          -perc_identity 100 \\
                                          -task "megablast" \\
                                          -num_threads 19 \\
                                          -max_hsps 100 \\
                                          -word_size {width}'))
            
            gb.i <- gb.i %>% dplyr::filter(Alignment.Length == width & Perc.Ident == Perc.Ident)
            
            aso.genome.blast <- rbind(aso.genome.blast, gb.i)
            gb.i <- NULL
            
        }
        
        human.genome.hit.df <- aso.genome.blast %>% 
            filter(grepl(pattern = "^NC_", x = SubjectID)) %>% # Only use chromosomal sequences 
            group_by(QueryID) %>% summarise(genome_hits =n()) %>% 
            mutate(QueryID = as.character(QueryID)) %>% 
            arrange(desc(genome_hits))
        
        aso.df <- aso.df %>% left_join(human.genome.hit.df, by = c('aso_idx' = 'QueryID'))
        blast.result[[2]] <- aso.genome.blast
    }
    
    
    return(list(aso.df = aso.df, blast.result = blast.result))
    
}


# 3. find ortholog 

ASOGeneratorFunction03 <- function(aso.df, orghologs = c('monkey.fa','mouse.fa','rat.fa','dog.fa')){
    
    mp.jun.l <- list()
    pb <- progress_bar$new(format = "  Processing [:bar] :current/:total :percent eta: :eta",
                           total = nrow(aso.df), 
                           clear = FALSE)
    
    for(i in 1:nrow(aso.df)){
        pb$tick()
        mp <- vmatchPattern(pattern = aso.df$seq[i] , subject = seq.files.all)
        if(length(unlist(mp)) == 0){ next }
        mp.jun.df <- as.data.frame(unlist(mp))
        mp.jun.df$aso_idx <- aso.df$aso_idx[i]
        mp.jun.l[[i]] <- mp.jun.df    
    }
    mp.jun.l.df <- rbindlist(mp.jun.l) %>% as_tibble()
    
    
    aso.df$mouse.perfect.match <- ifelse(aso.df$aso_idx %in% (mp.jun.l.df %>% filter(names %in% names(mouse.fa.all)) %>% distinct(aso_idx) %>% pull(aso_idx)), 'match','mismatch')
    aso.df$monkey.perfect.match <- ifelse(aso.df$aso_idx %in% (mp.jun.l.df %>% filter(names %in% names(monkey.fa.all)) %>% distinct(aso_idx) %>% pull(aso_idx)), 'match','mismatch')
    aso.df$rat.perfect.match <- ifelse(aso.df$aso_idx %in% (mp.jun.l.df %>% filter(names %in% names(rat.fa.all)) %>% distinct(aso_idx) %>% pull(aso_idx)), 'match','mismatch')
    aso.df$dog.perfect.match <- ifelse(aso.df$aso_idx %in% (mp.jun.l.df %>% filter(names %in% names(dog.fa.all)) %>% distinct(aso_idx) %>% pull(aso_idx)), 'match','mismatch')
    
    
    return(aso.df)

    
}





system("blastn 
       -version")










