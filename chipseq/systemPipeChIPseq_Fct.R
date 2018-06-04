f##############################                                                                                                                                               
## PWM matching in sequences ##                                                                                                                                               
###############################                                                                                                                                               
                                                                                                                                                                              
library (MotifDb)                                                                                                                                                             
as.list(query(MotifDb, 'athaliana'))                                                                                                                                          
as.list(query(MotifDb, 'hsapiens'))                                                                                                                                           
chr <- DNAStringSet(c(seq1="AAAGCTAAAGGTAAAGCAAAACGCCGCCG", seq2="CGCCGCCG"))                                                                                                 
pwm <- as.list(query(MotifDb, 'athaliana'))[[1]]                                                                                                                              
sapply(chr, function(x) matchPWM(pwm, x, min.score="90%"))                                                                                                                    
                                                                                                                                                                              
##########################                                                                                                                                                    
## Matches among motifs ##                                                                                                                                                    
##########################                                                                                                                                                    
                                                                                                                                                                              
library (MotIV)                                                                                                                                                               
?motifMatch                                                                                                                                                                   
egr1.motif <- MotifDb [egr1.mouse.jaspar.rows]                                                                                                                                
ts <- motifMatch (as.list (pwm), as.list (pwm), top=11)   

#####################################
## Motif enrichment with PWMEnrich ##
#####################################
library(MotifDb)
library(Biostrings) 
library(PWMEnrich)
library(GenomicFeatures)

peak_seq <- readDNAStringSet('results/SRR038848_1.fastq_trim.gz_C.bam_macs2_peaks.xls.fasta',format='fasta')
athaliana_motifs <- as.list(query(MotifDb, 'athaliana'))

# PWMEnrich needs pfm as integer matrix 
pfm_int_mat <- function(pfm){
	pfm <- as.integer(as.matrix(pfm)*10000000)
	pfm <- matrix(pfm, nrow=4)
	rownames(pfm) <- c('A', 'C', 'G', 'T')
	colnames(pfm) <- 1:dim(pfm)[2]
	return(pfm)
}

# Convert MotifDb athaliana motifs into integer matrix 
pfm <- sapply(athaliana_motifs, pfm_int_mat)

# Ned promoter sequence for calculating background frequencyi
txdb <- loadDb("./data/TAIR10.sqlite")
promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)
genomeFa <- readDNAStringSet('data/tair10.fasta')
promoter_seq<- getSeq(genomeFa, promoters[seqnames(promoters)=='Chr1' | seqnames(promoters)=='Chr2' | seqnames(promoters)=='Chr3' | seqnames(promoters)=='Chr4' | seqnames(promoters)=='Chr5' | seqnames(promoters)=='ChrC' | seqnames(promoters)=='ChrC'])

# Motif enrichment needs LogN background
PWM_bg <- makePWMLognBackground(clean(promoter_seq), pfm)

# Motif enrichment
res <- motifEnrichment(peak_seq[1], toy_bg)
#


