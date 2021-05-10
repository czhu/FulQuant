#### annotation free cluster filtering
library(tidyverse)
library(rtracklayer)
library(parallel)
projectFolder = "."
SCRIPTDIR = file.path(projectFolder, "sw/")
GENOMEDIR = file.path(projectFolder, "genome")
source(file.path(SCRIPTDIR,"clustering_functions.R"))

###################################################################################################
infolder = file.path(projectFolder, "combined")
load(file.path(infolder,"all.clustered.rda"))

outfolder=infolder
if(!file.exists(outfolder)) dir.create(outfolder)

outfile = file.path(outfolder, "cluster_filtered.rda")
message( "Outputting to ", outfile )

#####################################################################################################
##### filter reads
####################################################################################################
myclips = read_clipping_profile_file(file.path(infolder,"clip_profile_at_ends.gz"))
myinserts = read_insert_profile_near_ss_file(file.path(infolder,"insert_profile_near_ss.gz"))

### aggregate myclips and myinserts
myclipsDF = myclips %>% group_by(qname) %>% summarise(max_clip_width=max(width))
myinsertsDF = myinserts %>% group_by(qname) %>% summarise(max_insert=max(width_left + width_right))

ncpu=20 ## not too many cpu to avoid mem overflow on the server
myreads$clname = pvec(myreads$rcluster,standardise_clname,mc.cores=ncpu)

######### readtab contains the information for filtering later
readtab = tibble(
    qname = myreads$name,
    start = start(myreads),
    end = end(myreads),
    strand = as.character(strand(myreads)),
    tx_strand = pvec(qname, extract_txStrand_from_qname, mc.cores=ncpu),
    mappedLength= sum(width(blocks(myreads))),
    qscore = pvec(qname, extract_qscore_from_qname, mc.cores=ncpu),
    isfull = is_fulllength(qname),
    clname = myreads$clname
)

readtab = readtab %>% left_join(myclipsDF,by="qname") %>% left_join(myinsertsDF,by="qname")
readtab$max_clip_width = fill_na_with_value(readtab$max_clip_width)
readtab$max_insert = fill_na_with_value(readtab$max_insert)

##### filtering
minQual = 6
clipCutoff = 30
insertSSSizeCutoff=5

readtab = readtab %>% mutate(
    pass_qual = qscore>minQual,
    pass_clip = max_clip_width < clipCutoff,
    pass_insert= max_insert < insertSSSizeCutoff,
    pass_clusterable = is_clustered(clname),
    isReadPass = pass_qual & pass_clip & pass_insert & pass_clusterable
)

readStatOutfile = file.path(infolder,"read_filter_stat.txt.gz")

sink(file(readStatOutfile,"w"),type="message")
message(mean(readtab$pass_qual)," after qscore filter\n")
message(mean(readtab$pass_clip)," after clip filter\n")
message(mean(readtab$pass_insert)," after insert filter\n")

message(sum(readtab$pass_clip & readtab$pass_clusterable)," after clip filter\n")
message(sum(readtab$pass_insert & readtab$pass_clusterable)," after insert filter\n")
message(sum(readtab$pass_clusterable)," can be clustered\n")

message(mean(readtab$pass_clusterable)," can be clustered\n")

message(mean(readtab$isReadPass)," after all filters\n")
message(sum(readtab$isReadPass)," after all filters\n")
sink(type="message")

### save the essential to
stopifnot(all(myreads$clname==readtab$clname))

mcols(myreads)  = cbind(mcols(myreads), readtab[,setdiff(colnames(readtab),c("qname","start","end","strand","clname"))])
save(myreads, file=file.path(infolder,"myreads.rda"))

readsAfterFilter = myreads[myreads$isReadPass]
message("Working with ",length(readsAfterFilter)," reads\n")

save(readsAfterFilter, file=file.path(infolder,"readsAfterFilter.rda"))

####################################################################################################
### filter at the cluster level
## factor here to make sure it's the same of coercing
##chrname should not contain "-" which is used to concatenate
## proivde readsAfterFilter and outfolder this is a function

spltFac = factor(readsAfterFilter$clname, unique(readsAfterFilter$clname))
readCluster = split(readsAfterFilter, spltFac)
nFullLength = sum(as(split(readsAfterFilter$isfull,spltFac),"LogicalList"))

## we exclude these dubious cases due to very small exon , and hence start>end after clustering to splice sites
myclusterStrandPolyAType = infer_strand_clusteredReads(readCluster, minPercDiff=0.3,method="polyAType")
myclusterStrandSplicingMotif = infer_strand_clusteredReads(readCluster, minPercDiff=0.3,method="splicingMotif")
strandBoth = paste(myclusterStrandPolyAType,myclusterStrandSplicingMotif,sep=":")
readClusterCensus = infer_consensus(readCluster)

strand(readClusterCensus) = myclusterStrandSplicingMotif
readClusterCensus$strand_evidence = ifelse(
    strandBoth %in% c("-:-","+:+"), "both",
        ifelse(strandBoth %in% c("-:*","+:*"),"polyA",
            ifelse(strandBoth %in% c("*:-","*:+"),"motif",
                ifelse( strandBoth %in% c("-:+","+:-"), "conflicting", "none")))
    )

usePolyAStrand = readClusterCensus$strand_evidence %in% c("polyA","conflicting")
strand(readClusterCensus)[usePolyAStrand] = myclusterStrandPolyAType[usePolyAStrand]
names(readClusterCensus) = readClusterCensus$name

### clustered trascrnitp table
stopifnot(identical(names(readCluster),readClusterCensus$name))
stopifnot(identical(names(readCluster),names(nFullLength)))

polyAMaskFile = file.path(GENOMEDIR,"polyA_mask_seed5_windowSize20_cutoff60_recovered.txt.gz")
d_start_polyARich = dist_to_polyARich_region(readClusterCensus,"start",ignore.strand=TRUE,polyAmaskFile)
d_end_polyARich = dist_to_polyARich_region(readClusterCensus,"end",ignore.strand=TRUE,polyAmaskFile)
distPolyARich = dist_to_polyA_rich_region_strand_specific(readClusterCensus, polyAMaskFile)
distPolyARich[is.na(distPolyARich)] = Inf

readClusterCensus$n_full_length = nFullLength
readClusterCensus$d_polyARich = distPolyARich
readClusterCensus$count = lengths(readCluster)
stopifnot(identical(readClusterCensus$name,names(readCluster)))

save(readClusterCensus, file=file.path(outfolder, "readClusterCensus.rda"))

clusterTab = tibble(
    clname = readClusterCensus$name,
    count = readClusterCensus$count,
    n_full_length = readClusterCensus$n_full_length,
    perc_full_length = n_full_length/count,
    chr = as.character(seqnames(readClusterCensus)),
    start = start(readClusterCensus),
    end = end(readClusterCensus),
    strand = as.character(strand(readClusterCensus)),
    d_polyARich=readClusterCensus$d_polyARich
)
message(nrow(clusterTab)," total clusters\n")
####
## CHANGEME
## filtering criteria
clusterFilterTab = tibble(
    pass_count = clusterTab$count >= 3,
    pass_perc_full_length = clusterTab$perc_full_length >= 0.4,
    pass_polyArich = clusterTab$d_polyARich >= 10,
    pass_strand = clusterTab$strand !="*"
)

message(mean(clusterFilterTab$pass_count)," after min count\n")
message(mean(clusterFilterTab$pass_perc_full_length)," after min perc full length\n")
message(mean(clusterFilterTab$pass_polyArich)," after polyA rich region\n")
message(mean(clusterFilterTab$pass_strand)," strand has info\n")

isClusterPass = rowSums(clusterFilterTab)==ncol(clusterFilterTab)

message(mean(isClusterPass)," total ", sum(isClusterPass), " after all filters\n")

## we introduce a new class for filtering tx cluster
names(readClusterCensus) = NULL
mytxCluster = list(
    consensus = readClusterCensus,
    stat = clusterTab,
    filter = clusterFilterTab
)
mytxCluster = find_tu(mytxCluster)

###################################################################################################
##### reference
mytxClusterAdaptive = adaptive_filter(mytxCluster,0.03,ncpu=ncpu)
mytxClusterAdaptive = find_tu(mytxClusterAdaptive)

mytxClusterAdavanced1 = advanced_filter1(mytxClusterAdaptive, minDist=20, ncpu=ncpu,
    smallExonSize=60,deltaDist=20)
mytxClusterAdavanced1 = find_tu(mytxClusterAdavanced1)


mytxClusterAdavanced2 = advanced_filter2(mytxClusterAdavanced1, cutoff =0.25, ncpu=ncpu)
mytxClusterAdavanced2 = find_tu(mytxClusterAdavanced2)

save(mytxCluster, mytxClusterAdaptive, mytxClusterAdavanced1, mytxClusterAdavanced2,
    file=outfile)
