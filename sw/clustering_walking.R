extract_ss_from_gr = function(tx){
    ## tx is granges object with blocks
    myblocks = blocks(tx)
    nexons = lengths(myblocks)
    idx = cumsum(nexons)

    allstarts = unlist(start(myblocks))
    allends = unlist(end(myblocks))

    myseqnames = rep(as.character(seqnames(tx)), nexons-1)

    qei = GRanges( myseqnames, IRanges(allends[-idx], width=1) )
    qie = GRanges( myseqnames, IRanges(allstarts[-c(1,(idx[-length(idx)]+1))],width=1) )
    qei$qname = qie$qname = rep(tx$name, nexons-1)

    list( qei=qei, qie=qie )
}

get_cluster_optimised = function(x){
    ssSites = extract_ss_from_gr(x)

    nintrons = lengths(blocks(x))-1
    names(nintrons) = x$name

    ovlpsEI = as.matrix(findOverlaps(ssSites$qei, rei,type="within",ignore.strand=myignoreStrand))
    ovlpsIE = as.matrix(findOverlaps(ssSites$qie, rie,type="within",ignore.strand=myignoreStrand))

    stopifnot(all(ssSites$qei$qname == ssSites$qie$qname))

    ovlpsEIans = list(qname = ssSites$qei$qname[ovlpsEI[,1]], refidx = ovlpsEI[,2])
    ovlpsEIansSplt = split(ovlpsEIans$refidx, factor(ovlpsEIans$qname, unique(ssSites$qei$qname)))

    ovlpsIEans = list(qname = ssSites$qie$qname[ovlpsIE[,1]], refidx = ovlpsIE[,2])
    ovlpsIEansSplt = split(ovlpsIEans$refidx, factor(ovlpsIEans$qname, unique(ssSites$qie$qname)))

    isClustered = which((lengths(ovlpsEIansSplt) == nintrons[names(ovlpsEIansSplt)]) &
        (lengths(ovlpsIEansSplt) == nintrons[names(ovlpsIEansSplt)]))

    ans = rep("unclustered", length(x))
    names(ans) = x$name

    if(length(isClustered)){
        ## intron interval FIXME try to use the block class?
        prenames = paste(
            rei$name[ unlist(ovlpsEIansSplt[isClustered]) ],
            rie$name[ unlist(ovlpsIEansSplt[isClustered]) ],
            sep="-")
        myqnames = rep(names(ovlpsEIansSplt[isClustered]),lengths(ovlpsEIansSplt[isClustered]))
        myclnames= tapply(prenames, factor(myqnames,unique(myqnames)),
            paste, collapse=";")

        ans[names(myclnames)] = myclnames
    }
    names(ans)=NULL
    ans
}

library(tidyverse)
library(rtracklayer)

myignoreStrand = TRUE
ncpu=30

args = commandArgs(trailingOnly = TRUE)
#
# ### two inputs both as bed files, ei and ie splice sites
infile = args[1] ## bed file

eiFile = args[2] ## ref ei file
ieFile = args[3] ## ref ie file

rei = import(eiFile)
rie = import(ieFile)

rei$name = paste(seqnames(rei), rei$name,sep="_")
rie$name = paste(seqnames(rie), rie$name,sep="_")

fileBasename= file.path(dirname(infile), strsplit(basename(infile),"\\.")[[1]][1])
outfile = paste0(fileBasename,".clustered.rda")
outfile_clustered = paste0(fileBasename,".clustered.txt.gz")
outfile_singleton = paste0(fileBasename,".singleton.txt.gz")
outfile_unclustered = paste0(fileBasename,".unclustered.txt.gz")

#if(!file.exists(outfile)) {
message("Reading ",basename(fileBasename),"......\n")
myreads = import(infile)

if(anyDuplicated(myreads$name)){
    stop("Reads cannot map to multiple loci. This confuses the clustering")
}

myreads$rcluster="unclustered"

###check singleton
isSingleton = which(elementNROWS(myreads$blocks)==1)
myreads$rcluster[isSingleton] = "singleton"

testIdx = which(elementNROWS(myreads$blocks)>1)

message("Clustering ",basename(fileBasename),"......\n")

myreads$rcluster[testIdx] = unlist(mclapply(split(testIdx, cut(seq_len(length(testIdx)), ncpu)), function(i) {
    x = myreads[i]

    get_cluster_optimised(x)

}, mc.cores=ncpu))

save(myreads,file=outfile)

## export clustered reads as a table
readToClusterMapping = tibble(qname=myreads$name, clname=myreads$rcluster,
    chr=as.character(seqnames(myreads)),
    rstart=as.integer(start(myreads)), rend= as.integer(end(myreads)),
     nexon=lengths(blocks(myreads)) )

isUnclusterd = readToClusterMapping$clname == "unclustered"
isSingleton = readToClusterMapping$clname == "singleton"

write_tsv(readToClusterMapping[(!isUnclusterd) & (!isSingleton), ], outfile_clustered)
write_tsv(readToClusterMapping[isUnclusterd, ], outfile_unclustered)
write_tsv(readToClusterMapping[isSingleton, ], outfile_singleton)
#}



# get_cluster = function(thisReads){
#     unlist(mclapply(1:length(thisReads), function(i){
#         x = thisReads[i]
#         mystarts  = start(x$blocks[[1]])
#         myends = end(x$blocks[[1]])
#         qei = GRanges(seqnames(x),IRanges(start(x) + myends[-length(myends)]  , width=1))
#         qie = GRanges(seqnames(x),IRanges(start(x) + mystarts[-1]-1, width=1))
#
#         ovlpsEI = as.list(findOverlaps(qei, rei,type="within",ignore.strand=myignoreStrand))
#         ovlpsIE = as.list(findOverlaps(qie, rie,type="within",ignore.strand=myignoreStrand))
#         isCluster = all(elementNROWS(ovlpsEI) == 1) & all(elementNROWS(ovlpsIE) == 1)
#         if(isCluster){
#             ## intron interval FIXME try to use the block class?
#              paste(paste(rei$name[unlist(ovlpsEI)],rie$name[unlist(ovlpsIE)],sep="-"),collapse=";")
#          } else {
#             "unclustered"
#         }
#     },mc.cores=ncpu))
# }
