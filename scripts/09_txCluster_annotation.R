## annotate tx clusters using known annotation

library(rtracklayer)
library(tidyverse)
projectFolder = "."
SCRIPTDIR = file.path(projectFolder, "sw")
GENOMEDIR = file.path(projectFolder, "genome")
source(file.path(SCRIPTDIR,"clustering_functions.R"))

infolder = file.path(projectFolder, "combined")

## export reads supporting novel tx
doExtractReads = FALSE

### reference clusters with no mutil hits
load(file.path(GENOMEDIR, "tx.rda"))

### txCluster after filtering
load(file.path(infolder, "cluster_filtered.rda"))

thisTxCluster = mytxClusterAdavanced2
thisTxCluster$stat$flag = calc_bitcode_flag(!thisTxCluster$filter)

outfolder=file.path(infolder,"tx_annot")
if(!file.exists(outfolder)) dir.create(outfolder)

####################################################################################################
thisTxCluster$consensus$clname = thisTxCluster$consensus$name
stopifnot(all(is_clustered(tx$clname)))
stopifnot(all(is_clustered(thisTxCluster$consensus$clname)))
stopifnot(anyDuplicated(tx$clname)==0)

matchToRef = match_to_ref(thisTxCluster$consensus, tx)

thisTxCluster$stat$associated_tx = thisTxCluster$stat$match_type  = as.character(NA)
thisTxCluster$stat$delta_size =  as.integer(NA)

thisTxCluster$stat$associated_tx[matchToRef$query_i] = tx$transcript_id[matchToRef$subject_i]
thisTxCluster$stat$match_type[matchToRef$query_i] = matchToRef$type
thisTxCluster$stat$delta_size[matchToRef$query_i] = matchToRef$delta_size

################################################################################
isToSwapStrand = thisTxCluster$consensus$strand_evidence == "conflicting" & nblocks(thisTxCluster$consensus) > 2
warning(sum(isToSwapStrand), " strand conflicting cases need to be swapped")
strand(thisTxCluster$consensus)[isToSwapStrand] = swap_strand(as.character(strand(thisTxCluster$consensus)[isToSwapStrand]))
thisTxCluster$consensus$strand_evidence[isToSwapStrand] = "fixed"

### for tx matching annotation if the strand disagree with annotation we take the strand from annotation
isMatchToRef = !is.na(thisTxCluster$stat$associated_tx)
strandFromRef = strand(tx)[match(thisTxCluster$stat$associated_tx[isMatchToRef], tx$transcript_id)]
isDiffFromAnnot = strand(thisTxCluster$consensus)[isMatchToRef] != strandFromRef
warning(sum(isDiffFromAnnot), " cases strand disagree with reference" )
strand(thisTxCluster$consensus)[isMatchToRef][isDiffFromAnnot] = strandFromRef[isDiffFromAnnot]
thisTxCluster$consensus$strand_evidence[isMatchToRef][as.logical(isDiffFromAnnot)] = "annot"

################################################################################
### match to genes this can only happen after the strand is fixed
#ovlps = as.list(findOverlaps(thisTxCluster$consensus, mygenes))
matchToGenes = overlap_by_perc(thisTxCluster$consensus, mygenes, ignore.strand=FALSE, simplify=TRUE)

thisTxCluster$stat$associated_gene = as.character(NA)
thisTxCluster$stat$associated_gene[matchToGenes$query_i] = mygenes$gene_id[matchToGenes$subject_i]

convert_name_x_to_y_with_na = function(x,y){
    ans = x
    ans[!is.na(x)] = y[ans[!is.na(x)]]
    fill_na_with_value(ans,"")
}
#### make display names
matchingType = c("exact_match"="EX","fuzzy_match_as"="FS","fuzzy_match_se"="FE")
thisTxCluster$stat$displayName = paste(
        paste0("F",thisTxCluster$stat$flag),
        paste0("C",thisTxCluster$stat$count),
        fill_na_with_value(thisTxCluster$stat$associated_tx,""),
        convert_name_x_to_y_with_na(thisTxCluster$stat$match_type, matchingType),
        fill_na_with_value(thisTxCluster$stat$associated_gene,""),
        thisTxCluster$stat$clname,
        sep=":")

thisTxCluster$clusterType = tibble(
        falsetso = has_flag(thisTxCluster$stat$flag,64),
        tx_sig_diff = has_flag(thisTxCluster$stat$flag,32),
        adaptive = has_flag(thisTxCluster$stat$flag,16),
        #complex_region = has_flag(thisTxCluster$stat$flag,16),
        known_strand = has_flag(thisTxCluster$stat$flag,8),
        polyA_rich =has_flag(thisTxCluster$stat$flag,4) ,
        perc_full_length = has_flag(thisTxCluster$stat$flag,2),
        min_count = has_flag(thisTxCluster$stat$flag,1),
        ## match_tx_pass = !is.na(thisTxCluster$stat$associated_tx) & (thisTxCluster$stat$flag==0),
        ## if it's matching we no longer require to pass all the filter anymore
        match_tx_pass = !is.na(thisTxCluster$stat$associated_tx) & thisTxCluster$stat$flag == 0,
        match_tx_fail = !is.na(thisTxCluster$stat$associated_tx) & (thisTxCluster$stat$flag!=0),
        novel_pass = is.na(thisTxCluster$stat$associated_tx) & (thisTxCluster$stat$flag==0)
        )
thisTxCluster$stat$is_human_final = (thisTxCluster$clusterType$match_tx_pass |
    thisTxCluster$clusterType$match_tx_fail | thisTxCluster$clusterType$novel_pass) &
    !is_control(thisTxCluster$consensus)

thisTxCluster$mycolorScheme = c(
    "#377eb8", ## falsetso
    "#00CDCD", ## cyan3 tx_sig_diff
    "#000000", ## black adaptive
    #"#984ea3",
    "#ff7f00", ## known_strand
    "#a65628", ## polyA_rich
    "#f781bf", ## perc_full_length
    "#999999", ## min_count
    "#4daf4a", ## match_tx_pass
    "#984ea3", ## match_tx_fail
    "#ff7f00" ## novel_pass
)

names(thisTxCluster$mycolorScheme) = colnames(thisTxCluster$clusterType)
plotLegend = TRUE
if(plotLegend) {
    pdf(file.path(outfolder,"mylegend.pdf"),width=5,height=8)
    plot(rep(1,length(thisTxCluster$mycolorScheme)),
        1:length(thisTxCluster$mycolorScheme),col=thisTxCluster$mycolorScheme,pch=19,cex=5,xlim=c(0,4), xaxt='n', yaxt='n',xlab="",ylab="")
    text(rep(2.5,length(thisTxCluster$mycolorScheme)),1:length(thisTxCluster$mycolorScheme),cex=1.2,names(thisTxCluster$mycolorScheme))
    dev.off()
}

thisTxCluster$stat$name= thisTxCluster$stat$displayName
thisTxCluster$stat$score = 128 - thisTxCluster$stat$flag ## 8 bit
thisTxCluster$stat$igv_coord = granges_to_igvCoord(thisTxCluster$consensus)

thisTxCluster$stat$itemRgb = rep(as.character(NA),nrow(thisTxCluster$clusterType))
for(i in 1:ncol(thisTxCluster$clusterType)){
    thisTxCluster$stat$itemRgb[thisTxCluster$clusterType[[i]]] = thisTxCluster$mycolorScheme[i]
}
stopifnot(!anyNA(thisTxCluster$stat$itemRgb))

stopifnot(all(thisTxCluster$consensus$strand_evidence[thisTxCluster$stat$is_human_final] != "none"))
stopifnot(!any(strand(thisTxCluster$consensus)[thisTxCluster$stat$is_human_final] == "*"))

#### collapsing multi machting of tx cluster to reference using fuzzy matching
thisTxCluster$stat$grouping = ifelse(!is.na(thisTxCluster$stat$associated_tx),thisTxCluster$stat$associated_tx,thisTxCluster$stat$clname)

reduce_to_consensus = function(x, ncpu=10){
    # x = thisTxCluster
    annot = as_tibble(cbind(x$stat, x$clusterType[,c("match_tx_pass","match_tx_fail", "novel_pass")]))
    selAnnotColNames = setdiff(colnames(annot), c("chr", "start","end","strand","name","clname"))
    mygr = x$consensus
    mcols(mygr) = cbind(mcols(mygr), annot[,selAnnotColNames])
    ## data for multiple matches to  trace back
    mygr$grand_count = mygr$count
    mygr$n_cl = 1L
    mygr$all_cl = mygr$name

    mygr$name = mygr$displayName
    countPerTx = table(annot$grouping)
    isSingle = annot$grouping %in% names(countPerTx)[countPerTx==1]
    isMulti = annot$grouping %in% names(countPerTx)[countPerTx > 1]

    ansSingle = mygr[isSingle,]
    ansSingle$is_multi = FALSE

    rgList = split(mygr[isMulti], annot$grouping[isMulti])
    mypaste = function(x) paste(x, collapse=",")

    message(length(rgList), " units to process")
    ## impment in a way to sparate the cases no merge and merge
    ## since it's only minority that need merge, the speed gain is significant
    rv = mclapply((1:length(rgList)), function(i) {
        thisRg = rgList[[i]]
        ## in case of multiple tx cluster, take the one with highest expression
        ## if multiple highest expression take the one with flag ==0
        if(length(thisRg)>1){
            isMax = thisRg$count == max(thisRg$count)
            nmaxCount = sum(isMax)
            if ( nmaxCount==1) {
                wh = which.max(thisRg$count)
            } else {
                if (any(thisRg$flag[which(isMax)]==0) ) {
                    wh = which(isMax & thisRg$flag==0)
                    if(length(wh)>1){
                        wh = sample(wh,1)
                    }
                } else {
                    wh = sample(which(isMax),1)
                }
            }
        } else {
            wh = 1
        }
        ans = thisRg[wh]
        ans$grand_count = sum(thisRg$count)
        ans$n_cl = length(thisRg)
        ans$all_cl = mypaste(thisRg$name)
        return(ans)
    }, mc.cores=ncpu)
    ### NOTE this takes for ever, the impementation needs to be changed
    ansMulti = unlist(as(rv,"GRangesList"))
    ansMulti$is_multi = TRUE
    res = c(ansSingle, ansMulti)
    stopifnot( anyDuplicated(res$grouping)==0 )
    names(res)=NULL
    res[order(res)]
}

myclusters = reduce_to_consensus(thisTxCluster,ncpu=30)

export(myclusters[!is.na(myclusters$associated_tx) | myclusters$flag==0],
    file.path(outfolder,"clusters_putative.bed.gz"))

export(
    subset(myclusters, is_human_final),
    file.path(outfolder,"clusters_human.bed.gz")
)

export(myclusters,file.path(outfolder,"clusters_full.bed.gz"))

save(thisTxCluster,myclusters, file=file.path(outfolder,"txCluster.rda"))

########################################################################
### NOTE since the change of code I may still need to adapter here
if(doExtractReads) {
    myfolder = file.path(PROJECTDIR, "longread/clustering")
    load(file.path(myfolder, "readsAfterFilter.rda"))

    thisQnames = readsAfterFilter$name[readsAfterFilter$clname %in% thisTxCluster$stat$clname[thisTxCluster$clusterType$novel_pass]]
    extract_reads_by_qname(thisQnames,
        bamfile=file.path(myfolder,"all.bam"),
        outfile=file.path(outfolder,"reads_cluster_novel.bam"),
        wait=FALSE)

    thisQnames = readsAfterFilter$name[readsAfterFilter$clname %in%
        thisTxCluster$stat$clname[thisTxCluster$clusterType$match_tx_pass | thisTxCluster$clusterType$match_tx_fail]]

    extract_reads_by_qname(thisQnames,
        bamfile=file.path(myfolder,"all.bam"),
        outfile=file.path(outfolder,"reads_cluster_known.bam"),
        wait=FALSE)
}
