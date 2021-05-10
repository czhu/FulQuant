library(nanopore)
library(doMC)

folder = "."

outfolder = file.path(folder, "fastq_trimmed")
if(!file.exists(outfolder)) dir.create(outfolder)

infolder1 = file.path(folder, "fastq", "demultiplexed")
infolder2 = file.path(folder, "fastq", "trimmed")

infiles = dir(infolder1, "*fastq.gz", full.names=TRUE)

### filtering
minimalResidualLength = 200
minimalMeanQscore = 6
saveLog=TRUE
readsChunkSize=1e6


registerDoMC(cores=round(length(infiles)/2))

foreach(readsFile=iter(infiles)) %dopar% {
    bn = sub(".fastq.gz","", basename(readsFile))
    message("Working on ", bn)
    ssFile = file.path(infolder2, paste0(bn, "_tab_for_trimming.txt.gz"))

    outTrimmedReads = file.path(outfolder,paste0(bn, ".fastq"))
    outDiscardedReads = file.path(outfolder,paste0(bn, "_reads_dicarded.fastq"))
    processlogFile = file.path(outfolder,  paste0(bn, "_trimming_pipeline.log"))
    unlink(c(outTrimmedReads,outDiscardedReads,processlogFile))

    tsv = read_tsv(ssFile,col_types=list( "end.right"="d","start.right"="d", "end.left"="d",
        "start.left"="d", start="d",end="d"))

    trimTab = tsv %>% mutate(
        ## use polyT if possible
        trim_start = ifelse(
            grepl("^.-T-.$",conf), end + 1,
            ifelse(grepl("^P-.-.$",conf), end.left + 1, 1)),
        trim_end = ifelse(
            grepl("^.-A-.$",conf), start -1,
            ifelse(grepl("^.-.-P$",conf), start.right - 1, read_length)),
        residue_length =  trim_end - trim_start + 1,
        read_newid = paste(read_id, conf,trim_start,trim_end,residue_length,
                read_length,sprintf("Q%.02f", mean_qscore),sep=":"),
        qc_pass = residue_length > minimalResidualLength & mean_qscore >= minimalMeanQscore)

    ## initialise
    fs= FastqStreamer(readsFile,n=readsChunkSize)
    ## read in chunks with FastqStream
    if(saveLog) {
        sink(file=file(processlogFile, open = "wt"),type="message")
    }
    nRead = 0
    nFiltered = 0
    nTrimed = 0
    while (length(fq <- yield(fs))) {
        message("Processing ", length(fq)," Reads starting at read ", nRead+1)

        inputSeq = tibble(read_id = sapply(strsplit(as.character(ShortRead::id(fq))," "),"[",1))
        tabForTrimming = left_join(inputSeq,trimTab,by="read_id")  ## !!same order as fq!!

        ### filtered out reads fail the QC
        if(any(!tabForTrimming$qc_pass)){
            whReadsToDiscard = which(!tabForTrimming$qc_pass)
            outSeq = change_header(fq[whReadsToDiscard],tabForTrimming$read_newid[whReadsToDiscard])

            writeFastq( outSeq, file=outDiscardedReads,
                mode=ifelse(file.exists(outDiscardedReads),"a","w"), compress=FALSE)

            message("Filtered ", length(outSeq)," due to QC")
            nFiltered = nFiltered + length(outSeq)
        }
        ###################### write out processed reads
        if(any(tabForTrimming$qc_pass)) {
            wh = which(tabForTrimming$qc_pass)
            outSeq = change_header(
                narrow(fq[wh],tabForTrimming$trim_start[wh], tabForTrimming$trim_end[wh]),
                tabForTrimming$read_newid[wh])

            writeFastq(outSeq, file=outTrimmedReads,
                mode=ifelse(file.exists(outTrimmedReads),"a","w"), compress=FALSE)
            nTrimed = nTrimed + length(outSeq)
            message("Processed ", length(outSeq)," reads")
        }
        nRead = nRead + length(fq)
    }

    message(
        "Final stats:\n",
        "Processed ", nRead, " in total\n",
        nFiltered," filtered due to QC\t", nFiltered/nRead," perc\n",
        nTrimed," retained after QC\t", nTrimed/nRead," perc\n",
        "trimmed after QC\t",
        sum((trimTab$residue_length != trimTab$read_length)[trimTab$qc_pass])/sum(trimTab$qc_pass),
        " perc"
    )
    close(fs)
    if((nTrimed+nFiltered)!=nRead) { message("we lost reads somewhere, please check\n")}
    system(paste("pigz -p 10", outTrimmedReads))
    system(paste("pigz -p 10", outDiscardedReads))
    message("Done")

}
