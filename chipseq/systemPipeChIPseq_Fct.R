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
library(ShortRead)
library(BSgenome)

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

# Need promoter sequence for calculating background frequency
txdb <- loadDb("./data/TAIR10.sqlite")
promoters <- promoters(genes(txdb), upstream = 1000, downstream = 0)
#promoters <- sample(promoters, 5000, replace=FALSE)
genomeFa <- readDNAStringSet('data/tair10.fasta')
#promoter_seq <- getSeq(genomeFa, promoters)
promoter_seq<- getSeq(genomeFa, promoters[seqnames(promoters)=='Chr1' | 
		      seqnames(promoters)=='Chr2' | seqnames(promoters)=='Chr3' | 
		      seqnames(promoters)=='Chr4' | seqnames(promoters)=='Chr5'])

# Motif enrichment needs LogN background
sample_promoter_seqs <- sample(promoter_seq, 1000, replace=FALSE)
PWMLogn.tair10.MotifDb.Athal<- makePWMLognBackground(clean(promoter_seq), pfm)
save(PWM_bg, file='PWMLogn.tair10.MotifDb.Athal')
load('Athal.TAIR10.PWM.bg')

# Motif enrichment for single peak sequences
res <- sapply(peak_seq, function(x) motifEnrichment(x, PWM_bg))
report <- sequenceReport(res[[1]], 1)
plot(report[report$p.value < 0.05], fontsize=7, id.fontsize=6)

# Motif enrichment in multiple sequences
res <- motifEnrichment(peak_seq, PWM_bg)
report <- groupReport(res)
report

png(filename='group_enriched_motif.png')
plot(report[1:10], fontsize=7, id.fontsize=5)
dev.off()
