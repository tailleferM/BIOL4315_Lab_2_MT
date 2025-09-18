library(Rqc)
folder <- system.file(package="ShortRead","extdata/E-MTAB-1147")
qcRes <- rqc(path = folder, pattern = ".fastq.gz", openBrowser = FALSE)
rqcCycleQualityBoxPlot(qcRes)
rqcCycleBaseCallsLinePlot(qcRes)
rqcReadFrequencyPlot(qcRes)

library(fastqcr)

#Aggregating Multiple FastQC Reports into a Data Frame

#Demo QC directory containing zipped FASTQC reports
qc.dir <- system.file("fastqc_results", package = "fastqcr")
qc <- qc_aggregate(qc.dir)
qc

# Inspecting QC Problems

# See which modules failed in the most samples
qc_fails(qc, 'module')
# Or, see which samples failed the most
qc_fails(qc, 'sample')

# Building Multi QC Reports
qc_report(qc.dir, result.file = 'multi-qc-report')

# Building One-Sample QC Reports (+ Interpretation)
qc.file <- system.file('fastqc_results', 'S1_fastqc.zip', package = 'fastqcr')

#view the report rendered by R functions
qc_report(qc.file, result.file = 'one-sample-report', interpret = TRUE)


#2.6
library(QuasR)
# obtain a list of fastq file paths
fastqFiles <- system.file(package='ShortRead', 'extdata/E-MTAB-1147',
                          c("ERR127302_1_subset.fastq.gz",
                            "ERR127302_2_subset.fastq.gz"))
#defined processed fastq file names
outfiles <- paste(tempfile(pattern=c("processed_1_",
                                     "processed_2_")),".fastq",sep="")
# process fastq files
# remove reads that have more than 1 N, (nBases)
# trim 3 bases from the end of reads (truncateEndBases)
# Remove ACCCGGGA pattern if it occurs at the start (Lpattern)
# remove reads shorted than 40 base-pairs (minLength)
preprocessReads(fastqFiles, outfiles,
                nBases=1,
                truncateEndBases=3,
                Lpattern='ACCCGGGA',
                minLength=40)

library(ShortRead)

#obtain a list pf fastq file paths 
fastqFile <- system.file(package="ShortRead",
                         "extdata/E-MTAB-1147",
                         "ERR127302_1_subset.fastq.gz")
#read fastq file
fq = readFastq(fastqFile)
#get quality scores per base as a matrix
qPerBase = as(quality(fq),"matrix")
#get number of bases per read that have quality score below 20
#we use this 
qcount = rowSums(qPerBase <= 20)
#Number of reas where all Phred scores >= 20
fq[which(qcount == 0)]
#We can finally write out the filtered fastq file with the ShortRead::writeFastq() function
#mode = 'a' allows you to rewrite files that are already written
ShortRead::writeFastq(fq[which(qcount == 0 )], 
                      paste(getwd(), "Qfiltered3.fastq", sep='/'),mode ='a')

#set up streaming with block size 1000
# every time we call the yield() function 1000 read portion 
# of the file will be read successively
f <- FastqStreamer(fastqFile,readerBlockSize = 1000)
#we set up a while loop to call yield() function to go through the file
while(length(fq <- yield(f))){
  #remove reads where all quality scores are <20
  #get quality scores per base as a matrix
  qPerBase = as(quality(fq),'matrix')
  #get number of bases per read that have Q score < 20
  qcount = rowSums(qPerBase <= 20)
  #write fastq file with mode 'a', so every new block is
  #written to the same file
  writeFastq(fq[which(qcount == 0)],
             paste(fastqFile,'Qfiltered', sep='_', mode='a'),mode= 'a')
}

# working on question 9 
library(ShortRead)
fastqFiles9 <- c('/Users/mtaillefer00/ERR11203340_2.fastq.gz','/Users/mtaillefer00/ERR11203340_1.fastq.gz')
cJujuni = readFastq(fastqFiles9)
library(Rqc)
QCJujuni <- rqc(path = '/Users/mtaillefer00/lab2q9',pattern = '.fastq.gz', openBrowser = FALSE)
rqcCycleQualityBoxPlot(QCJujuni)
#The quality score drops below 20 at cycle 101
rqcCycleBaseCallsLinePlot(QCJujuni)
#the proportion of Ns is constant, at 0
rqcReadFrequencyPlot(QCJujuni)
#~45% of reads occur once, ~25% occur twice and 30% occur more than twice

# working on question 10
# truncateEndBases()
