## trimming off adapters

library(nanopore)

folder = "."

infolder = file.path(folder, "fastq", "demultiplexed")
outfolder = file.path(folder, "fastq", "trimmed")
if(!file.exists(outfolder)) dir.create(outfolder)

infiles = dir(infolder, "*fastq.gz", full.names=TRUE)

adapter = c("forward"="TTTCTGTTGGTGCTGATATTGC", "reverse"="ACTTGCCTGTCGCTCTATCTTC")

windowSize = 300
readsChunkSize = 1e6
ncpu = 30 ## FIXME this as a option
gapOpening = 1
gapExtension = 3
get_pwa_cutoff = function(M,adapter,FDR){
    require(Biostrings)
    randomStrings = apply(matrix(sample(DNA_BASES, nchar(adapter) * M, replace = TRUE),
                   nrow = M), 1, paste, collapse = "")

    s = pairwiseAlignment(DNAStringSet(randomStrings), DNAString(adapter),
        type="overlap", gapOpening=gapOpening,gapExtension=gapExtension,scoreOnly=TRUE)
    quantile(s,1-FDR)
}
mycutoffs = sapply(adapter, get_pwa_cutoff, M=1e5,FDR=0.005)
if( length(unique(mycutoffs)) ==1 ){
    mycutoff = mycutoffs[1]
} else{
    stop("Changing script for different adapter cutoffs")
}

get_pwa_hit = function(x) {
    ## x is tibble from match_adapter
    ## only work for one adapter
    x %>% filter(score>mycutoff)
}
extract_res_tibble = function(x,suffix="",columnSeletcted){
    ## x is tibble
    if(missing(columnSeletcted))
        columnSeletcted= c("read_id","start","end","score","match_length","pid","adapter_id")
    rv = x[,columnSeletcted]
    colnames(rv) = paste0(colnames(rv),c("",rep(suffix,ncol(rv)-1)))
    return(rv)
}

for(infile in infiles){
    message("Working on ", basename(infile))
    outfile = file.path(outfolder, paste0(sub(".fastq.gz","", basename(infile)),"_identify_adapter.txt.gz"))
    if(file.exists(outfile)) unlink(outfile)

    fs= FastqStreamer(infile,n=readsChunkSize)

    while (length(fqBatch <- yield(fs))) {
        ans = mclapply(split(fqBatch,cut(seq_len(length(fqBatch)),ncpu)),function(fq) {
            qs = get_quality(fq)
            readStats = tibble(
                read_id = sapply(strsplit(as.character(ShortRead::id(fq))," "),"[",1),
                median_qscore = median(qs),
                mean_qscore = calc_mean_qscore(qs),
                median_start_qscore = median(get_quality(extract_seq_end(fq, windowSize,fromEnd=FALSE))),
                median_end_qscore = median(get_quality(extract_seq_end(fq, windowSize,fromEnd=TRUE))),
                read_length = width(fq)
            )

            ## define trimming window
            leftAdapters = match_adapter(extract_seq_end(fq,windowSize=windowSize), adapter,use.names=TRUE,
                    gapOpening=gapOpening, gapExtension=gapExtension, avoidPartialMappingOnEnd=FALSE)
            leftAdapters = get_pwa_hit(leftAdapters)
            leftAdapters$read_id = readStats$read_id[leftAdapters$group]

            rightAdapters = match_adapter(
                reverseComplement(extract_seq_end(fq,windowSize=windowSize,fromEnd=TRUE)), adapter,
                use.names=TRUE,gapOpening=gapOpening, gapExtension=gapExtension,
                avoidPartialMappingOnEnd=FALSE)
            rightAdapters = get_pwa_hit(rightAdapters)
            rightAdapters$read_id = readStats$read_id[rightAdapters$group]
            ## convert back to read coordinates
            readRightAdapterStart = rightAdapters$start
            readRightAdapterEnd = rightAdapters$end
            rightAdapters$start =  readStats$read_length[rightAdapters$group] - readRightAdapterEnd + 1
            rightAdapters$end =  readStats$read_length[rightAdapters$group] - readRightAdapterStart + 1

            ##### polyA
            polyAByLength = find_polyA(fq, n=10, max.mismatch=2L,hitBy="length",use.names=TRUE)
            polyAByLength$read_id = readStats$read_id[polyAByLength$group]

            polyAByDist = find_polyA(fq, n=10, max.mismatch=2L,hitBy="distanceToEnd",use.names=TRUE)
            polyAByDist$read_id = readStats$read_id[polyAByDist$group]

            usefulColsPolyA = c("read_id","start","end","adapter_id","allStart","allEnd","allAdapter_id")
            readStats %>% left_join(extract_res_tibble(leftAdapters,".left"),by="read_id") %>%
                left_join(extract_res_tibble(rightAdapters,".right"), by="read_id") %>%
                left_join(extract_res_tibble(polyAByLength,".polyA",usefulColsPolyA), by="read_id") %>%
                left_join(extract_res_tibble(polyAByDist,".polyA.d",c("read_id","start","end","adapter_id")),
                    by="read_id")
            },mc.cores=ncpu)
        ans= do.call(rbind,ans)

        write_tsv(ans,outfile, append=file.exists(outfile))
    }
    close(fs)
}
