library(nanopore)
library(parallel)
ncpu = 30

make_full_table=function(x){
    myunlist = function(y) unlist(strsplit(y,","))
    isDefined = rowSums(is.na(x[,c("allStart.polyA","allEnd.polyA")]))==0
    subx = x[isDefined,]
    nEach = sapply(strsplit(subx$allStart.polyA,","),length)

    rv = tibble(
        read_id = rep(subx$read_id,nEach),
        read_length = rep(subx$read_length,nEach),
        start = as.integer(myunlist(subx$allStart.polyA)),
        end = as.integer(myunlist(subx$allEnd.polyA)),
        type = myunlist(subx$allAdapter_id.polyA)
    )
    rv
}

folder = "."

infolder = file.path(folder, "fastq", "trimmed")
infiles = dir(infolder, "*_identify_adapter.txt.gz", full.names=TRUE)

windowSize = 300
minimalPolyALength = 10  ## 10 is no filtering
distAllowedToAdapter = c(-5,10)  ## for polyA
chooseByLength = TRUE

#infile = "identify_adapter.txt.gz"
for(infile in infiles) {
    message("Working on ", basename(infile))
    infile = normalizePath(infile)
    outfolder = dirname(infile)

    outfile = file.path(outfolder, paste0(sub("_identify_adapter.txt.gz","", basename(infile)),"_tab_for_trimming.txt.gz"))

    tsv = read_tsv(infile,col_types=list(allStart.polyA="c",allEnd.polyA="c",start.right="d",end.right="d"))

    polyATab = make_full_table(tsv)

    polyATab = polyATab %>% mutate(
        length = end - start + 1,
        dist_to_end = ifelse(type == "polyA", read_length - end, start - 1),
        is_in_end = dist_to_end < windowSize
        ) %>% filter(length>= minimalPolyALength, is_in_end)

    res = left_join(tsv[,c("read_id","median_qscore","mean_qscore", "read_length", "start.left",
        "end.left", "start.right", "end.right" )], polyATab[ , -which(colnames(polyATab) == "read_length")],
        by="read_id")

    ### determining which is the right polyA/T signal
    rvDistAdapter = res %>% mutate(
        dist_to_adapter = ifelse(type=="polyA",start.right - end, start - end.left),
        is_near_adapter = (dist_to_adapter > distAllowedToAdapter[1]) & (dist_to_adapter < distAllowedToAdapter[2]),
        polyA_comp_type = ifelse(
            is.na(type),"no_polyA",ifelse(type=="polyA",
            ifelse(is.na(start.right),"no_right_adapter","polyA_adapter"),
            ifelse(is.na(start.left),"no_left_adapter","polyA_adapter")))
        )
    ### remove false positive polyA/T if we know adapter positions
    ### if adapter position unknown it won't be affected
    selFalsePolyA = which(!rvDistAdapter$is_near_adapter)
    rvDistAdapter$start[selFalsePolyA] = rvDistAdapter$end[selFalsePolyA] =
        rvDistAdapter$length[selFalsePolyA] = as.integer(NA)
    rvDistAdapter$type[selFalsePolyA] = as.character(NA)

    myf = unique(rvDistAdapter$read_id)
    ans = mclapply(split(myf, cut(seq_len(length(myf)),ncpu)), function(thisf) {
        res = rvDistAdapter %>% filter(read_id %in% thisf) %>% group_by(read_id) %>% summarise(
                nhit = sum(is_near_adapter,na.rm=TRUE),
                hit = if( n()==1 ) { 1 }
                      else {
                        if(nhit==1) {
                            which(is_near_adapter) }
                        else if (nhit>1) {
                            wh = which(is_near_adapter)
                            if(chooseByLength) {
                                wh[which.max(length[wh])]
                            }
                            else {
                                wh[which.min(dist_to_end[wh])]
                            }
                        }
                        else { ## in case of no hit
                            wh = which(!is.na(type))
                            if(length(wh)==0) 1
                            else if (length(wh)==1) {wh}
                            else {
                                if(chooseByLength) {
                                    wh[which.max(length[wh])]
                                }
                                else {
                                    wh[which.min(dist_to_end[wh])]
                                }
                            }
                        }
                      },
                mean_qscore=mean_qscore[hit],
                read_length = read_length[hit],
                start.left = start.left[hit],
                end.left = end.left[hit],
                start.right = start.right[hit],
                end.right = end.right[hit],
                start = start[hit],
                end = end[hit],
                type= type[hit],
                length= length[hit],
                n=n())
        },mc.cores=ncpu)
    ans = do.call(rbind,ans)

    ## at this point we have a tab ans with adapter positions and polyA positions
    ## one read can only have one adapter at the start, one adapter at the end and one polyA signal
    ## it's polyT when it's at the start of the read, it's polyA when it's at the end of the read
    ## now we need to obtain left trim position and right trim positions, and read configuration
    ans = ans %>% mutate(
        has_polyA = !is.na(type),
        has_left_adapter = !is.na(start.left),
        has_right_adapter = !is.na(start.right),
        conf = paste(
            ifelse(has_left_adapter, "P", "N"),
            ifelse(has_polyA, ifelse(type=="polyA","A","T"), "N"),
            ifelse(has_right_adapter, "P", "N"),
            sep="-")
        ) %>% select(-one_of("nhit","hit","has_polyA","has_left_adapter","has_right_adapter"))

    ans$conf[which(ans$end.left > ans$start.right)] = "N-N-N"

    write_tsv(ans, outfile)
}
