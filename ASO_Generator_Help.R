#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# ASO Generator Help Functions
# 2022-12-17 HSY
#




#### Search Functions ####

asoModuleFinder <-  function(aso_start,aso_end, module_starts, module_ends, module_names){
    # Search aso within module range and indicate overlap 
    module_list <- ifelse(aso_start >= module_starts,
                          ifelse(aso_start <=  module_ends, 
                                 ifelse(aso_end <= module_ends, module_names , 
                                        paste("Start_with",module_names)), 
                                 ""), 
                          ifelse(aso_end >= module_starts, 
                                 ifelse(aso_end > module_ends, 
                                        glue("aso include {module_names}"), 
                                        paste("Ends_with",module_names) ), 
                                 ""))
}

asoModuleFinder.onlynames <-  function(aso_start,aso_end, module_starts, module_ends, module_names){
    
    module_list <- ifelse(aso_start >= module_starts,
                          ifelse(aso_start <=  module_ends, 
                                 ifelse(aso_end <= module_ends, module_names , 
                                        module_names), 
                                 ""), 
                          ifelse(aso_end >= module_starts, 
                                 ifelse(aso_end > module_ends, 
                                        module_names, 
                                        module_names ), 
                                 ""))
    module_list
}


add.seq <- function(seq.start.chr.position, seq = gene.seq.character, seq.start.idx= gene.genomic.start, aso.width = width ){
    substr(seq,(seq.start.chr.position - seq.start.idx + 1),(seq.start.chr.position - seq.start.idx + aso.width))
}



addJunctionalASO <- function(gene.gtf.df, transcript, aso.width= width){
    
    # Create Junctional aso, 
    
    exons <- gene.gtf.df %>% filter(transcript_name == transcript & type =='exon')
    exons$start <- as.numeric(exons$start)
    exons$end <- as.numeric(exons$end)
    exons$exon_number <- as.numeric(exons$exon_number)
    junction.aso <- data.frame()
    if(max(exons$exon_number) > 1){
        for(i in 1:(max(exons$exon_number)-1) ){
            j1_start = exons$end[i] - aso.width + 2
            j2_end = exons$start[i+1] + aso.width - 2
            junction.aso.i <- data.frame(j1.start = j1_start:exons$end[i], 
                                         j1.end = exons$end[i],
                                         j2.start = exons$start[i+1],
                                         j2.end =  exons$start[i+1]:j2_end,
                                         position = paste0(exons$transcript_name[i],'_',exons$exon_number[i], '_', exons$exon_number[i+1])
            )
            junction.aso.i <- junction.aso.i %>% mutate(width = j1.end - j1.start + j2.end - j2.start + 2) %>% filter(width == aso.width)
            junction.aso <- rbind(junction.aso, junction.aso.i)
        }
        
    }
    return(junction.aso)
}

createJunctionalASO <- function(gene.gtf.df, gene.gtf.df.transcript, aso.width, gene.seq.character, gene.genomic.start){
    
    aso.df.junction <- data.frame()
    for(i in 1:nrow(gene.gtf.df.transcript)){
        junctional.aso.i <- addJunctionalASO(gene.gtf.df, 
                                             gene.gtf.df.transcript$transcript_name[i], 
                                             aso.width = aso.width)
        aso.df.junction <- rbind(aso.df.junction, junctional.aso.i)
        
    }
    
    aso.df.junction$j1.seq <- mapply(add.seq, seq.start.chr.position = aso.df.junction$j1.start, 
                                     seq = gene.seq.character, seq.start.idx= gene.genomic.start,
                                     aso.width =  (aso.df.junction$j1.end - aso.df.junction$j1.start +1 ) )
    
    aso.df.junction$j2.seq <- mapply(add.seq, seq.start.chr.position = aso.df.junction$j2.start, 
                                     seq = gene.seq.character, seq.start.idx= gene.genomic.start,
                                     aso.width =  (aso.df.junction$j2.end - aso.df.junction$j2.start +1 ) )
    
    aso.df.junction$seq <- paste0(aso.df.junction$j1.seq,aso.df.junction$j2.seq)
    
    aso.df.junction <- aso.df.junction %>% group_by(j1.start, j2.start, j1.end, j2.end, j1.seq, j2.seq, seq) %>% mutate(position = paste0(position, collapse = ',')) %>% distinct() %>% ungroup()
    
    aso.df.junction$aso_idx <- paste0('J',1:nrow(aso.df.junction))
    return(aso.df.junction)
    
}






