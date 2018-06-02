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
library('PWMEnrich')

peak_seq <- readDNAStringSet('results/SRR038848_1.fastq_trim.gz_C.bam_macs2_peaks.xls.fasta',format='fasta')
athaliana_motifs <- as.list(query(MotifDb, 'athaliana'))

# I need to check the pwm to pfm conversion
# Is it mathmetically sound?
pwm_to_pfm <- function(pwm){
	pfm <- as.integer(as.matrix(pwm)*100000)
	pfm <- matrix(pfm, nrow=4)
	rownames(pfm) <- c('A', 'C', 'G', 'T')
	colnames(pfm) <- 1:dim(pfm)[2]
	return(pfm)
}


############################### 
pfm_toy <- sapply(athaliana_motifs[1:5], pwm_to_pfm)
#
toy_bg <- makePWMLognBackground(peak_seq, pfm_toy)
#
res <- motifEnrichment(peak_seq[1], toy_bg)
#


