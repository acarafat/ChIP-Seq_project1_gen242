## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

#Loading a number of the needed libraries at the start
## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(systemPipeR)
    library(BiocParallel)
    library(Biostrings)
    library(Rsamtools)
    library(GenomicRanges)
    library(ggplot2)
    library(GenomicAlignments)
    library(ShortRead)
    library(ape)
    library(systemPipeRdata)
    library(systemPipeR)
    
  
})

## ----genChip_workflow, eval=FALSE----------------------------------------

# to ensure we have all the required and supporting files, git clone the group git directory to generate the work environment
system("git clone https://github.com/girke-class/gen242_2018_ChIP-Seq1") ##cheating?????????/





#I dont want to make a new bunch of files each time i run this R file, I have already used git clone
#genWorkenvir(workflow="chipseq")
#setwd("chipseq_proj")


### ----r_environment, eval=FALSE-------------------------------------------
###system("hostname") # should return name of a compute node starting with i or c
### getwd() # checks current working directory of R session
### dir() # returns content of current working directory



## ----load_systempiper, eval=TRUE-----------------------------------------
library(systemPipeR)

#There are no functions specifically referenced from the Fct file
### ----load_custom_fct, eval=FALSE-----------------------------------------
#source("systemPipeChIPseq_Fct.R")



#Specifying and displaying some of the targets file holding the ChIP seq reads
## ----load_targets_file, eval=TRUE----------------------------------------
targetspath <- system.file("./targets_chip.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4,-c(5,6)]



#Quality controlling the reads, dropping any read with a base below the specified PHRED quality score
#We decided to rerun the workflow without filtering due to losing so much data to filtering
## ----proprocess_reads, eval=FALSE, messages=FALSE, warning=FALSE, cache=TRUE----
# args <- systemArgs(sysma="param/trim.param", mytargets="targets_chip.txt")
# system("cat ")
# 
# #### quality control Make more permissive! 
# filterFct <- function(fq, cutoff=20, Nexceptions=0) {
#   qcount <- rowSums(as(quality(fq), "matrix") <= cutoff)
#   fq[qcount <= Nexceptions] # Retains reads where Phred scores are >= cutoff with N exceptions
#   }
#   preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
#   writeTargetsout(x=args, file="targets_chip_trim.txt", overwrite=TRUE)


#####Tophat aligner, we chose instead to use Bowtie 2
## ----fastq_report, eval=FALSE--------------------------------------------
#args <- systemArgs(sysma="param/tophat.param", mytargets="targets_chip.txt")
# library(BiocParallel); library(BatchJobs)
# f <- function(x) {
#     library(systemPipeR)
#     args <- systemArgs(sysma="param/tophat.param", mytargets="targets_chip.txt")
#     seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
# }
# funs <- makeClusterFunctionsSLURM("slurm.tmpl")
# param <- BatchJobsParam(length(args), resources=list(walltime="00:20:00", ntasks=1, ncpus=1, memory="2G"), cluster.functions=funs)
# register(param)
# fqlist <- bplapply(seq(along=args), f)
# pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
# seeFastqPlot(unlist(fqlist, recursive=FALSE))
#dev.off()


############Bowtie 2 aligner, our prefered aligner. Ran on a single system due to trouble with slurm.
## ----bowtie2_align, eval=FALSE-------------------------------------------
args <- systemArgs(sysma="param/bowtieSE.param", mytargets="targets_chip.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file
# moduleload(modules(args)) # Skip if a module system is not used
# 
# system("module load bowtie2")#trying to fix error "sh: bowtie2-build: command not found," get error:sh: module: command not found
# 
# system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta") # Indexes reference genome
# resources <- list(walltime="1:00:00", ntasks=1, ncpus=cores(args), memory="10G")
# reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01",
#                  resourceList=resources)
# waitForJobs(reg)


## ----bowtie2_align_seq, eval=FALSE---------------------------------------
runCommandline(args)
writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)#**need to move this down below runCommandline???


######check if the .bam files have been created
## ----check_files_exist, eval=FALSE---------------------------------------
file.exists(outpaths(args))


#############Write out an Excel file with alignent stats. compare to  alignStatsTheo.xls which was the result of alignment on the filtered data
## ----align_stats, eval=FALSE---------------------------------------------
read_statsDF <- alignStats(args=args)
write.table(read_statsDF, "./results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
read.delim("results/alignStats.xls")

read.delim("results/alignStatsTheo.xls")

###########creates the .bam files links for graphing  
## ----symbol_links, eval=FALSE--------------------------------------------
symLink2bam(sysargs=args, htmldir=c("~/.html/", "gen242/"),
            urlbase="http://biocluster.ucr.edu/~tkata002/", urlfile="results/IGVurl.txt")

####creating objects to store coverage information**only using one of the 7 outpaths args??
## ----rle_object, eval=FALSE----------------------------------------------
library(rtracklayer); library(GenomicRanges); library(Rsamtools); library(GenomicAlignments)
aligns <- readGAlignments(outpaths(args)[1])
cov <- coverage(aligns)
head(cov)

####plots coverage for all chromosomes using the merged bam 
## ----plot_coverage, eval=FALSE-------------------------------------------

#the par for the for loop isnt working, the par command seems to have no effect on autoplot output. So, I have to export these graphs one by one from Rstudio.
#par(mfrow=c(2,4))

#chr lengths for graphing obtained from: https://www.arabidopsis.org/portals/genAnnotation/gene_structural_annotation/agicomplete.jsp
#to get this to run, had to excise and run non-interactively on the high memory computational node. **
library(ggbio)
chr <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chr_length <- c(34964571,22000000,26000000,21000000,31000000)
for(i in chr){
  z <- as.character(i)
  pdf(z)
  q = 1
  myloc <- c(i, 1, chr_length[q])
  a <- readGAlignments(outpaths(args)[1], use.names=TRUE, param=ScanBamParam(which=GRanges(myloc[1], IRanges(as.numeric(myloc[2]), as.numeric(myloc[3])))))
  autoplot(a, aes(color = strand, fill = strand), facets = strand ~ seqnames, stat = "coverage")
  dev.off()
  q <- q+1
}



####combines the bam files from our multiple fastq files
## ----merge_bams, eval=FALSE----------------------------------------------
 args <- systemArgs(sysma=NULL, mytargets="targets_bam.txt")
 outpaths(args)
 args_merge <- mergeBamByFactor(args, overwrite=TRUE)
 writeTargetsout(x=args_merge, file="targets_mergeBamByFactor.txt", overwrite=TRUE)

 file.exists(outpaths(args_merge)) 

 #####peak calling without referenes, Had to do this in R on cluster directly line by line. This uses the control sample as reference to normalize by the random background
## ----call_peaks_macs_noref, eval=FALSE-----------------------------------
args <- systemArgs(sysma="param/macs2_noinput.param", mytargets="targets_mergeBamByFactor.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file



moduleload(module="python", envir=c("PATH", "LD_LIBRARY_PATH", "PYTHONPATH")) # Temp solution due to Python path change
  moduleload(module="python")

#system("module load python")
 
runCommandline(args)
file.exists(outpaths(args))
writeTargetsout(x=args, file="targets_macs.txt", overwrite=TRUE)


#we do not use the peak calling with reference for analysis, as we want to compare the AP1 and control peaks.
###############peak calling with reference#######################
#to get this to run, had to go directly to the culster and ran MACS2 directly with the code from the args_input file 
#macs2 callpeak   -t /bigdata/gen242/shared/ChIP-Seq1/results/SRR038848_1.fastq_trim.gz_C.bam -c /bigdata/gen242/shared/ChIP-Seq1/results/SRR038845_1.fastq_trim.gz_AP1.bam -n /bigdata/gen242/tkata002/5_30/gen242_2018_ChIP-Seq1/chipseq/results/SRR038848_1.fastq_trim.gz_C.bam_macs2 -f BAM -g 1.2e8 -B -q 0.01
###############################################################3
# ## ----call_peaks_macs_withref, eval=FALSE---------------------------------
# writeTargetsRef(infile="targets_mergeBamByFactor.txt", outfile="targets_bam_ref.txt", silent=FALSE, overwrite=TRUE)
# args_input <- systemArgs(sysma="param/macs2.param", mytargets="targets_bam_ref.txt")
# sysargs(args_input)[1] # Command-line parameters for first FASTQ file
# #unlink(outpaths(args_input)) # Note: if output exists then next line will not be run
# runCommandline(args_input)
# file.exists(outpaths(args_input))
# writeTargetsout(x=args_input, file="targets_macs.txt", overwrite=TRUE) #XX


###########consensus peaks, dont understand currently we have AP1_1 and C_1A, not M1A and A1A???
####################this calls peaks.annotated.xls but we have not annotated yet????
# ## ----consensus_peaks, eval=FALSE-----------------------------------------
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/rangeoverlapper.R")
peak_C <- outpaths(args)["C_1A"] ####currently only has "C_1A" now!!!!
peak_C <- as(read.delim(peak_C, comment="#")[,1:3], "GRanges") ##!!! why 1:3????



peak_A1A <- outpaths(args)["A1A"]
peak_A1A <- as(read.delim(peak_A1A, comment="#")[,1:3], "GRanges")

(myol1 <- subsetByOverlaps(peak_M1A, peak_A1A, minoverlap=1)) # Returns any overlap

myol2 <- olRanges(query=peak_C, subject=peak_C, output="gr") # Returns any overlap with OL length information
myol2[values(myol2)["OLpercQ"][,1]>=50] # Returns only query peaks with a minimum overlap of 50%






##############
# ## ----chip_peak_anno, eval=FALSE------------------------------------------
library(ChIPpeakAnno); library(GenomicFeatures)
args <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
#txdb <- loadDb("./data/tair10.sqlite")
txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana")
ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type"))
for(i in seq(along=args)) {
    peaksGR <- as(read.delim(infile1(args)[i], comment="#"), "GRanges")
    annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData=genes(txdb))
    df <- data.frame(as.data.frame(annotatedPeak), as.data.frame(values(ge[values(annotatedPeak)$feature,])))
    write.table(df, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")
}
writeTargetsout(x=args, file="targets_peakanno.txt", overwrite=TRUE) #XX
# 


############### Peak annotation using Biomart, IS THIS NECCESSARRY???????????????????
# ## ----chip_peak_anno_full_annotation, include=FALSE, eval=FALSE-----------
# Perform previous step with full genome annotation from Biomart
txdb <- makeTxDbFromBiomart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org")
tx <- transcripts(txdb, columns=c("tx_name", "gene_id", "tx_type"))
ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type")) # works as well
seqlevels(ge) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM")
table(mcols(tx)$tx_type)
tx <- tx[!duplicated(unstrsplit(values(tx)$gene_id, sep=","))] # Keeps only first transcript model for each gene]
seqlevels(tx) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM")
annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData = tx)

head(annotatedPeak)

##Did this, but there doesn't seem to be a write.table step so i have added it here:
df2 <- data.frame(as.data.frame(annotatedPeak))    #, as.data.frame(values(ge[values(annotatedPeak)$feature,])))
write.table(df2, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")



# ## ----chip_peak_seeker, eval=FALSE----------------------------------------
library(ChIPseeker)
for(i in seq(along=args)) {
    peakAnno <- annotatePeak(infile1(args)[i], TxDb=txdb, verbose=FALSE)
    df <- as.data.frame(peakAnno)
    write.table(df, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")
}
writeTargetsout(x=args, file="targets_peakanno.txt", overwrite=TRUE)#XX

# 
# ## ----chip_peak_seeker_plots, eval=FALSE----------------------------------
peak <- readPeakFile(infile1(args)[1])
covplot(peak, weightCol="X.log10.pvalue.")
peakHeatmap(outpaths(args)[1], TxDb=txdb, upstream=1000, downstream=1000, color="red")

pdf("peakHeatmap.pdf")
dev.off()

plotAvgProf2(outpaths(args)[1], TxDb=txdb, upstream=1000, downstream=1000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
# 
## ----count_peak_ranges, eval=FALSE---------------------------------------
library(GenomicRanges)
args <- systemArgs(sysma="param/count_rangesets.param", mytargets="targets_macs.txt")
args_bam <- systemArgs(sysma=NULL, mytargets="targets_bam.txt")
bfl <- BamFileList(outpaths(args_bam), yieldSize=50000, index=character())
countDFnames <- countRangeset(bfl, args, mode="Union", ignore.strand=TRUE)
writeTargetsout(x=args, file="targets_countDF.txt", overwrite=TRUE) #XX


########creates edgeR.xls files which contain differential expression infomration as well as other important values such as gene ids
## ----diff_bind_analysis, eval=FALSE--------------------------------------
args_diff <- systemArgs(sysma="param/rundiff.param", mytargets="targets_countDF.txt")
cmp <- readComp(file=args_bam, format="matrix")
dbrlist <- runDiff(args=args_diff, diffFct=run_edgeR, targets=targetsin(args_bam),
                    cmp=cmp[[1]], independent=TRUE, dbrfilter=c(Fold=.2, FDR=.05))
writeTargetsout(x=args_diff, file="targets_rundiff.txt", overwrite=TRUE) #XX
# 


##?????? i dont think this is working for me
#####looking for this file: SRR038848_1.fastq_trim.gz_C.bam_macs2_peaks.annotated.xls
# ## ----go_enrich, eval=FALSE-----------------------------------------------
args <- systemArgs(sysma="param/macs2.param", mytargets="targets_bam_ref.txt")
args_anno <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
annofiles <- outpaths(args_anno)

file.exists(outpaths(args_anno)) ##problem, this file does exist, bit ive been unable to make it work myself (SRR038848_1.fastq_trim.gz_C.bam_macs2_peaks.annotated.xls)

his

#gene_ids <- sapply(names(annofiles), function(x) unique(as.character(read.delim(annofiles[x])[,"geneId"])), simplify=FALSE)#XX

gene_ids <- unique(as.character(genes))#using  my own gene ids which are supplied from the annotated xls file created by the Parse_peak_sequences step


load("data/GO/catdb.RData")
BatchResult <- GOCluster_Report(catdb=catdb, setlist=gene_ids, method="all", id_type="gene", CLSZ=2, cutoff=0.9,  recordSpecGO=NULL) #gocats=c("MF", "BP", "CC"),

## ----parse_peak_sequences, eval=FALSE------------------------------------
library(Biostrings); library(seqLogo); library(BCRANK)
args <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
rangefiles <- infile1(args)
for(i in seq(along=rangefiles)) {
  
  i = 1
    df <- read.delim(rangefiles[i], comment="#")
    peaks <- as(df, "GRanges")
    names(peaks) <- paste0(as.character(seqnames(peaks)), "_", start(peaks), "-", end(peaks))
    #peaks <- peaks[order(values(peaks)$X.log10.pvalue.,values(peaks)$, decreasing=TRUE)] this is different and I dont know why
    peaks <- peaks[order(values(peaks)$X.log10.pvalue., decreasing=TRUE)]
    pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
    names(pseq) <- names(peaks)
    writeXStringSet(pseq, paste0(rangefiles[i], "Theo.fasta"))
}
head(peaks)


############ de novo motif descovery, not our project
# ## ----bcrank_enrich, eval=FALSE-------------------------------------------
# set.seed(0)
# BCRANKout <- bcrank(paste0(rangefiles[1], ".fasta"), restarts=25, use.P1=TRUE, use.P2=TRUE)
# toptable(BCRANKout)
# topMotif <- toptable(BCRANKout, 1)
# weightMatrix <- pwm(topMotif, normalize = FALSE)
# weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
# pdf("results/seqlogo.pdf")
# seqLogo(weightMatrixNormalized)
# dev.off()
# 
# ## ----sessionInfo---------------------------------------------------------
# sessionInfo()



##############Challenge project part 1, FDR "Prioritize/rank peaks by FDR from differential binding analysis"
edgeR_table <- read.table("./results/SRR038848_1.fastq_trim.gz_C.bam_macs2_peaks.edgeR.xls")
####sorting by smallest FDR and then smallest p value
edgeR_table_sorted <- edgeR_table[order(edgeR_table$AP1.C_FDR, edgeR_table$AP1.C_PValue),]
head(edgeR_table_sorted)


# hist(edgeR_table_sorted$AP1.C_FDR)
# hist(edgeR_table_sorted$AP1.C_PValue)
# 
# scatter.smooth(edgeR_table_sorted$AP1.C_FDR, edgeR_table_sorted$AP1.C_PValue)# looks bad, many very high FDR
# 
# #sorting out fdr greater than .9
# edgeR_table_sorted_FDR_cutoff <- edgeR_table_sorted[edgeR_table_sorted$AP1.C_FDR<= 0.9,]
# dim(edgeR_table_sorted_FDR_cutoff)
# 
# scatter.smooth(edgeR_table_sorted_FDR_cutoff$AP1.C_FDR, edgeR_table_sorted_FDR_cutoff$AP1.C_PValue)
# 
# 
# FDR_sorted <- peaks[rownames(edgeR_table_sorted)]
# FDR_seq <- getSeq(FaFile('./data/tair10.fasta'), FDR_sorted)
# 
# peaks[rownames(edgeR_table_sorted)[1:5]]
# tail(edgeR_table_sorted)
# 
# min(edgeR_table$AP1.C_PValue)
# 
# unique(edgeR_table$AP1.C_FDR)
# sort1.hsb2 <- hsb2[order(read) , ]
# sort()


###Presenting the peak sequences from  Macs2 peak calling. These are the sequences which were highly enriched within the set of sequences returned from the Illumina sequencing
#the sequences are saved as character strings with names being the chromosome and range designation.
###################Challenge project part 2 "Parse peak sequences from genome"
chr_test <- readDNAStringSet("./results/SRR038848_1.fastq_trim.gz_C.bam_macs2_peaks.xls.fasta")

peak_seqs <- as.character(chr_test)

#class(peak_seqs)
head(peak_seqs)
ggplot(data = peak_seqs)
autoplot(chr_test)

##################### Challenge project part 3 "Determine which peak prioritization method shows the highest enrichment for any of the motifs available in the Jaspar database (motifDB)."
ann_table <- read.table("./results/SRR038848_1.fastq_trim.gz_C.bam_macs2_peaks.annotated.xls", sep="\t", header = T)
#head(ann_table)
genes <- ann_table$geneId
venn_targ <- (unique(genes[ann_table$fold_enrichment >= 1.8]))
length(venn_targ)

#these are the table 3 data from paper, the differentially expressed genes "249 genes (which we refer to as high-confidence targets) showed robust (>1.8-fold) differential expression" from microarray data
xl <- read.table("/bigdata/gen242/tkata002/5_30/gen242_2018_ChIP-Seq1/chipseq/genes.csv")
head(xl)
dim(xl)
head(ann_table$gene_id)

#finding out how many genes overlap between paper and our own workflow output from EdgeR
l1 <-levels(xl$V1)
l2 <-levels(ann_table$geneId)

intersect(l1,l2)
sum(l1%in%l2)
length(xl$V1)
length(ann_table$geneId)

intersect(l1,venn_targ)

intersect(a,b)
library (MotifDb)                                                                                                                                                             

as.list(query(MotifDb, 'athaliana'))                                                                                                                                          

#as.list(query(MotifDb, 'hsapiens'))

chr_test <- readDNAStringSet("./results/SRR038848_1.fastq_trim.gz_C.bam_macs2_peaks.xls.fasta")

chr <- DNAStringSet(c(seq1="TTTCCGAAAAGGTATCACATGCCAAGTTTGGCCTCAC", seq2="GGATGCAACACGAGGACTTCCCGGGAGGTCACCCATCCTAGTACTACTCTC"))##does this only search within individual sequences???                                                                                               
chr <- DNAStringSet(chr_test)

pwm <- as.list(query(MotifDb, 'athaliana'))[[1]]                                                                                                                              
pwmall <- as.list(query(MotifDb, 'athaliana'))


pwm*100000

y <- sapply(chr, function(x) matchPWM(pwm, x, min.score="70%", with.score = T))### finds motifs from our peak sequences



install.packages("MotIV")
library("MotIV")                                                                                                                                                               
biocLite("MotIV")


?motifMatch                                                                                                                                                                   
egr1.motif <- MotifDb [egr1.mouse.jaspar.rows]                                                                                                                                
ts <- motifMatch (as.list (pwm), as.list (pwm), top=11)  



library(MotifDb)
mdb <- MotifDb
matrices <- subset(mdb, dataSource=='JASPAR')

for(i in 1:length(genes)){
  
  
}

#matrices <- subset (mdb, dataSource=='UniPROBE')

jaspar.test.matrices <- subset(mdb, geneSymbol== as.character(genes[1]) & dataSource == 'JASPAR')

#all.egr1.matrices <- subset (mdb, geneId=='13653')



# tbl.motifDbExample <- data.frame(motifName=c("Mmusculus-jaspar2016-Ahr::Arnt-MA0006.1",
#                                               "Hsapiens-jaspar2016-FOXI1-MA0042.2",
#                                               "Hsapiens-jaspar2016-HLF-MA0043.2"),
#                                   chrom=c("chr1", "chr1", "chr1"),
#                                   start=c(1000005, 1000085, 1000105),
#                                   start=c(1000013, 1000092, 1000123),
#                                   score=c(0.85, 0.92, 0.98),
#                                   stringsAsFactors=FALSE)
# tbl.out <- associateTranscriptionFactors(MotifDb, tbl.motifDbExample, source="MotifDb",
#                                          expand.rows=TRUE)
# dim(tbl.out) # one new column ("geneSymbol"), no new rows



#########PWMEnrich
source("https://bioconductor.org/biocLite.R")
biocLite("PWMEnrich")
library("PWMEnrich")
library(Biostrings)
##########this supposedly loads what getBackgroundFrequencies is internally calling for arabadopsis

genomic.acgt = getBackgroundFrequencies("BSgenome.Athaliana.TAIR.TAIR10")

Prom1000bp.genomic.acgt = makePriors(Prom1000bp, 1)


library(PWMEnrich.Athaliana.background)
# load the pre-compiled lognormal background
data(PWMLogn.dm3.MotifDb.Dmel)
sequences = readDNAStringSet(system.file(package="PWMEnrich",
                                         dir="extdata", file="tinman-early-top20.fa"))
res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)


###########################generating graphs:


##FDR vs p val
hist(edgeR_table_sorted$AP1.C_FDR)
hist(edgeR_table_sorted$AP1.C_PValue)

hist(edgeR_table_sorted$AP1.C_logFC)

###Could make a variable with positive or negative for logFC and do t test to determine if there is significant diff between counts of pos, vs negative

v <- log(edgeR_table_sorted$AP1.C_FDR)

scatter.smooth(edgeR_table_sorted$AP1.C_FDR, edgeR_table_sorted$AP1.C_PValue)# looks bad, many very high FDR


#sorting out fdr greater than .9
edgeR_table_sorted_FDR_cutoff <- edgeR_table_sorted[edgeR_table_sorted$AP1.C_FDR<= 0.9,]
dim(edgeR_table_sorted_FDR_cutoff)

scatter.smooth(edgeR_table_sorted_FDR_cutoff$AP1.C_FDR, edgeR_table_sorted_FDR_cutoff$AP1.C_PValue)


#### transcript type breakdown
names(ann_table)
tx_breakdown <- barplot(prop.table(table(ann_table$tx_type)))
  
  

#####
hist(ann_table$width)
levels(ann_table$feature)
ann_table$
###### plot coverage
library(GenomicRanges)

set.seed(1); 
N <- 100; 
gr <- GRanges(seqnames = sample(c("chr1", "chr2", "chr3"), size = N, replace = TRUE), IRanges(start = sample(1:300, size = N, replace = TRUE), width = sample(70:75, size = N,replace = TRUE)), strand = sample(c("+", "-"), size = N, replace = TRUE), value = rnorm(N, 10, 3), score = rnorm(N, 100, 30), sample = sample(c("Normal", "Tumor"), size = N, replace = TRUE), pair = sample(letters, size = N, replace = TRUE))

#testing, our data doesnt appear to be strand specific, but this shows peaks on all chr
autoplot(peaks[Rle=="Chr1"], aes(color = strand, fill = strand), facets = strand ~ seqnames) # doesnt work
test <- autoplot(peaks, aes(color = strand, fill = strand), binwidth = 500) #makes the peaks by chrom graph
str(test)



autoplot(peaks, geom = "bar", aes(fill = peaks$fold_enrichment))

autoplot(a, layout = 'circle') # doenst work but a should be all reads for coverage

autoplot(peaks, aes(color = strand, fill = strand), facets = strand ~ seqnames)




###############getting untrimmed read numbers
system("grep  -c \"^>\" ./data/SRR038846_1.fastq.gz")
