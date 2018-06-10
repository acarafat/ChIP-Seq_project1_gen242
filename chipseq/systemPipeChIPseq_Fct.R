###############################                                                                                                                                               
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
## see here: http://bioconductor.org/packages/release/bioc/vignettes/PWMEnrich/inst/doc/PWMEnrich.pdf 

