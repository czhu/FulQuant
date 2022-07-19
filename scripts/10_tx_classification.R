### given a anntotation annotate exon
library(tidyverse)
library(rtracklayer)
library(parallel)

projectFolder = "."
SCRIPTDIR = file.path(projectFolder, "sw")
GENOMEDIR = file.path(projectFolder, "genome")
source(file.path(SCRIPTDIR,"clustering_functions.R"))

infolder = file.path(projectFolder, "combined/tx_annot")

## reference annotation clusters after cleaning up
load(file.path(GENOMEDIR, "tx.rda"))

## filtering results
load(file.path(infolder, "txCluster.rda"))

## novelTx = sample(novelTx,100) ## debug or development only
annotate_tx = function(thisTx, tx){
    thisTx$align_issue = "none"
    thisTx$is_end_variation = FALSE

    diffAtEitherEnd = is_different_at_either_end(thisTx, tx, simplify=TRUE)
    thisTx$is_end_variation[diffAtEitherEnd$query_i] = TRUE

    smallExonIssue = find_misaligned_small_exon_any(thisTx, tx, ignore.strand = FALSE,
        smallExonSize = 60, deltaDist = 20, simplify = TRUE)
    alignShiftIssue = find_misaligned_ss_shifted_any(thisTx, tx, ss_max_shift_size = 20, simplify = TRUE)
    thisTx$align_issue[alignShiftIssue$query_i] = "align_shift"
    thisTx$align_issue[smallExonIssue$query_i] = "small_exon"

    ## CHANGED 2018-07-10 report where the problematic intron is and how big the difference
    thisTx$delta_size  = as.integer(NA)
    ## for alignment shift there can be multiple hits and we combine by ","
    ## so the returned index is character not integer
    thisTx$ith_intron  = as.character(NA)

    thisTx$delta_size[alignShiftIssue$query_i] = as.integer(alignShiftIssue$delta_size)
    thisTx$delta_size[smallExonIssue$query_i] = as.integer(smallExonIssue$delta_size)

    thisTx$ith_intron[alignShiftIssue$query_i] = as.character(alignShiftIssue$intron_i)
    thisTx$ith_intron[smallExonIssue$query_i] = as.character(smallExonIssue$ith_intron)

    thisTx
}

isNotControl = !is_control(myclusters)
novelTx = subset(myclusters, novel_pass & isNotControl)
novelTx = annotate_tx(novelTx, tx)
novelClusterTab = suppressWarnings( as_tibble(novelTx) )


message(sum(novelClusterTab$is_end_variation), " novel tx due to ends variation")

message(sum(novelClusterTab$align_issue != "none"), " novel tx with align issue")

####################################################################################################
### classify
### comp exon to find out intron retention
### comp intron to find out exon skipping
####################################################################################################

### classify exon intron
knownExonsRanges = exonIV_to_intronIV(exonIV_to_intronIV(tx[nblocks(tx)>2]))
knownExonsRanges = unique(unlist(blocks(knownExonsRanges)))
## introduce new way of naming intervals
knownExonsRanges$name = make_ivname_from_granges(knownExonsRanges, includeStrand=TRUE)
stopifnot(!anyDuplicated(knownExonsRanges$name))

knownIntronsRanges = exonIV_to_intronIV(tx)
knownIntronsRanges = unique(unlist(blocks(knownIntronsRanges)))
## +1L because splice sites defined as last positive exon
## for consistency this is clname!
knownIntronsRanges$name = make_ivname_from_granges(knownIntronsRanges+1L, includeStrand=TRUE)
stopifnot(!anyDuplicated(knownIntronsRanges$name))


splt = unlist(strsplit(knownIntronsRanges$name,"-"))
knownEIJunctions = unique(splt[seq(1,length(splt),by=2)])
knownIEJunctions = unique(splt[seq(2,length(splt),by=2)])


##################################################################
## from our data for data
dataIntronsFull = exonIV_to_intronIV(novelTx)

########
dataIntronRanges = unlist( blocks(dataIntronsFull) )
dataIntronRanges$name =  make_ivname_from_granges( dataIntronRanges+1L, includeStrand=TRUE)

hasMultiExons = is_multiexonic(dataIntronsFull)
dataExonsFull = exonIV_to_intronIV(dataIntronsFull[hasMultiExons])
dataExonRanges = unlist(blocks(dataExonsFull))
dataExonRanges$name = make_ivname_from_granges(dataExonRanges, includeStrand=TRUE)

# myclnameList = strsplit(mytxIsoform$clname, ";")
# myclnameSplt = unlist(strsplit(unlist(myclnameList), "-"))

myclnameSplt = unlist(strsplit(dataIntronRanges$name,"-"))

nElements = nblocks(dataIntronsFull)
indexInTx = unlist(  sapply(nElements, seq_len) )

intronAnnot = tibble(
    clname= rep(novelClusterTab$clname, nElements ),
    name= dataIntronRanges$name,
    index_in_tx = indexInTx,
    pos_type = ifelse(
        rep(nElements==1, nElements) ,"singleton",
        ifelse(indexInTx==1, "first", ifelse( indexInTx == rep(nElements, nElements),"last","middle"))),
    ei = myclnameSplt[seq(1,length(myclnameSplt),by=2)],
    ie = myclnameSplt[seq(2,length(myclnameSplt),by=2)],
    is_ei_novel = !(ei %in% knownEIJunctions),
    is_ie_novel = !(ie %in% knownIEJunctions),
    is_novel = ! (name %in% knownIntronsRanges$name),
    type = ifelse(
        is_ei_novel & is_ie_novel,
            "novel",
            ifelse(
                !is_ei_novel & is_ie_novel, "ie_var", ifelse( is_ei_novel & !is_ie_novel, "ei_var", ifelse(!is_ei_novel & !is_ie_novel & is_novel, "recomb", "known"))
            )
    )
)

myclnameSplt = unlist( strsplit(dataExonRanges$name, "-") )
nElements = lengths( blocks(dataExonsFull) )

indexInTx = unlist( sapply(nElements, seq_len ))

exonAnnot = tibble(
    clname = rep(novelClusterTab$clname[ hasMultiExons ], nElements),
    name= dataExonRanges$name,
    index_in_tx = indexInTx,
    ie = myclnameSplt[seq(1,length(myclnameSplt),by=2)],
    ei = myclnameSplt[seq(2,length(myclnameSplt),by=2)],
    is_ei_novel = !(ei %in% knownEIJunctions),
    is_ie_novel = !(ie %in% knownIEJunctions),
    is_novel = !(name %in% knownExonsRanges$name),
    type = ifelse(
        is_ei_novel & is_ie_novel,
            "novel",
            ifelse(
                !is_ei_novel & is_ie_novel, "ie_var", ifelse( is_ei_novel & !is_ie_novel, "ei_var", ifelse(!is_ei_novel & !is_ie_novel & is_novel, "recomb", "known"))
            )
    )
)

intronAnnot$index = seq_len(nrow(intronAnnot))
exonAnnot$index = seq_len(nrow(exonAnnot))

pos_to_granges = function(posname){
    ## posname is like 18_49843392
    GRanges(
        sub("^(.*?)_(.*?)_(.*?)$","\\1" ,posname),
        IRanges(as.integer(sub("^(.*?)_(.*?)_(.*?)$","\\2" ,posname)),width=1)
    )
}

get_boundary = function(x, type="start"){
    ## x is Granges, return a GRanges
    if(!type %in% c("start","end")) {
        stop("Only start and stop are supported")
    }
    mypos = if(type=="start"){start(x)} else {end(x)}
    ans = GRanges(seqnames(x), IRanges(mypos, width=1), strand=strand(x))
    return(ans)
}


exonAnnot$dist_to_nearest = as.integer(NA)
exonAnnot$ee_name = ""

tmpvar = find_dist_at_one_end(dataExonRanges[exonAnnot$type=="ie_var"],knownExonsRanges,"start")
exonAnnot$dist_to_nearest[exonAnnot$type=="ie_var"][tmpvar$query_i] =  tmpvar$dist_to_nearest
exonAnnot$ee_name[exonAnnot$type=="ie_var"][tmpvar$query_i] =  tmpvar$name

tmpvar = find_dist_at_one_end(dataExonRanges[exonAnnot$type=="ei_var"],knownExonsRanges,"end")
exonAnnot$dist_to_nearest[exonAnnot$type=="ei_var"][tmpvar$query_i] =  tmpvar$dist_to_nearest
exonAnnot$ee_name[exonAnnot$type=="ei_var"][tmpvar$query_i] =  tmpvar$name

intronAnnot$dist_to_nearest = NA

wh = which(intronAnnot$pos_type %in% c("first","singleton") & intronAnnot$type=="ei_var")
thisX = get_boundary(dataIntronRanges[wh], type="start")
thisSubject = pos_to_granges(knownEIJunctions)
distToEI = suppressWarnings(as_tibble(distanceToNearest(
    thisX,
    thisSubject)))

intronAnnot$dist_to_nearest[wh] = distToEI$distance * ifelse(start(thisX[distToEI$queryHits]) > start(thisSubject[distToEI$subjectHits]),-1L,1L)


wh = which(intronAnnot$pos_type %in% c("last","singleton") & intronAnnot$type=="ie_var")
thisX = get_boundary(dataIntronRanges[wh], type="end")
thisSubject = pos_to_granges(knownIEJunctions)
distToIE = suppressWarnings(as_tibble(distanceToNearest(
    thisX,
    thisSubject)))

intronAnnot$dist_to_nearest[wh] = distToIE$distance * ifelse(end(thisX[distToIE$queryHits]) > end(thisSubject[distToIE$subjectHits]),1L,-1L)


findIncludingRanges = function(query, subject, ignore.strand = FALSE){
    ## find all subject that are fully included in query
    ovlps = suppressWarnings(findOverlaps( subject, query, type="within", ignore.strand=ignore.strand))
    ans = tibble(query_i = subjectHits(ovlps), subject_i = queryHits(ovlps), subject_width = width(subject)[subject_i])
    return(ans)
}

## Notice: the results are not simplied i.e. one query could have mutiple hits skipped multiple exons
### exon skipping events

intronsWithExonSkippingTab  = findIncludingRanges(dataIntronRanges, knownExonsRanges)
intronAnnot$skippedExon = c(tapply(
    knownExonsRanges$name[intronsWithExonSkippingTab$subject_i],
    factor(intronsWithExonSkippingTab$query_i,intronAnnot$index),
    paste,collapse=";"))
intronAnnot = intronAnnot %>% mutate(exon_skipped = !is.na(skippedExon), novel_exon_skipped = exon_skipped & type =="recomb")

############## intron retention
exonsWithIntronRetentionTab  = findIncludingRanges(dataExonRanges, knownIntronsRanges)
exonAnnot$retainedIntron = c(tapply(
    knownIntronsRanges$name[exonsWithIntronRetentionTab$subject_i],
    factor(exonsWithIntronRetentionTab$query_i,exonAnnot$index),
    paste,collapse=";"))
exonAnnot = exonAnnot %>% mutate(intron_retained = !is.na(retainedIntron), novel_intron_retained = intron_retained & type == "recomb")

## for debug only
knownIntronsRanges$igvCoord = granges_to_igv_coord(knownIntronsRanges)
knownExonsRanges$igvCoord  = granges_to_igv_coord(knownExonsRanges)


#################
ncpu=30
mypaste = function(x){paste(x, collapse=";")}

## classification according to discussion on April.23.2018
## table including number of
## novel_exon_allnew (both splites are not known)
## novel_exon_extention (only one splite sites is new)
## novel_exon_combnew (both splites are known but the combination is unknown, excluding novel_exon_skipping and novel_inron_retention)
## novel_utr_var (first and/or last exon splite sites change (because we do not know UTR length) changes)
## novel_exon_skipping (both splites sites need to be known)
## novel_intron_retention (both splites sites need to be known))
## novel_comb either 0 or 1 for FALSE, TRUE, if all exon are known, excluding ovel_exon_skipping and utr_variation
## tx type
## if novel_exon >0 -> "NE"
## else if novel_exon_extention >0 -> "ASS" (alternative splice site)
## else if novel_intron_retention -> "IR"
## else if novel_exon_skipping -> "ES"
## else if all exon and intron known -> "NC"
## else -> "others" -> manually check !!!
## novel_exon_combination if known , but all exon and introns are known
## others if none above applies
my_class_code = function(x){
    stopifnot(length(x) == 4)
    ## pos1: num of NE
    ## pos2: num of ASS
    ## pos3: num of IR
    ## pos4: num of ES
    ans = paste0(x,c("N","A","I","S"))
    ans = ifelse(x==0,  "",ans)
    paste(ans, collapse="")
}

is_multi_novel = function(x){
    # same x as in my_class_code
    stopifnot(length(x) == 4)
    sum(x>0) > 1
}

rvList = mclapply(novelClusterTab$clname, function(thisClname){
    thisIntron = intronAnnot %>% filter(clname==thisClname)
    thisExon = exonAnnot %>% filter(clname==thisClname)
    is_exon_ext_exon = thisExon$type %in% c("ei_var","ie_var")
    # is_utr_var = (thisIntron$pos_type %in% c("singleton","first") & thisIntron$type == "ei_var") |
    #     (thisIntron$pos_type %in% c("singleton","last") & thisIntron$type == "ie_var")

    ## FIXME this should not eat up new combination
    #is_utr_var = (thisIntron$pos_type %in% c("singleton","first","last") & thisIntron$type != "known")
    is_utr_var = (thisIntron$pos_type == "singleton" & thisIntron$type != "known") |
                (thisIntron$pos_type == "first" & thisIntron$type == "ei_var") |
                (thisIntron$pos_type == "last" & thisIntron$type == "ie_var")

    #is_utr_var_both = thisIntron$pos_type == "singleton" & thisIntron$type == "novel"
    is_novel_exon_combnew = thisExon$type == "recomb" & !thisExon$novel_intron_retained


    infoTab = tibble(
        clname = thisClname,
        intron_type = mypaste(thisIntron$type),
        exon_type = mypaste(thisExon$type),
        novel_exon_allnew = sum(thisExon$type == "novel"),
        novel_exon_allnew_names = mypaste(thisExon$name[thisExon$type == "novel"]),
        novel_exon_extension = sum(is_exon_ext_exon),
        ## known caveat if the alternative splice site event is a result of a terminal
        ## exon split into two exons, one of which has alternative splice site
        ## we won't be able to get the name of it because terminal exons is not in the knownExons
        ## which is internal exons only
        novel_exon_extension_names = mypaste(thisExon$ee_name[is_exon_ext_exon]),
        novel_utr_var = sum(is_utr_var), #+ sum(is_utr_var_both) * 2,
        novel_exon_combnew = sum(is_novel_exon_combnew),
        novel_exon_combnew_names = mypaste(thisExon$name[is_novel_exon_combnew]),
        intron_retention = sum(thisExon$intron_retained),
        novel_intron_retention = sum(thisExon$novel_intron_retained),
        exon_skipping= sum(thisIntron$exon_skipped),
        novel_exon_skipping = sum(thisIntron$novel_exon_skipped),
        intron_with_exon_skipping = mypaste(thisIntron$name[!is.na(thisIntron$skippedExon)] ),
        exon_with_intron_retention = mypaste(thisExon$name[!is.na(thisExon$retainedIntron)] ),
        #### novel cases
        novel_intron_with_exon_skipping = mypaste(thisIntron$name[!is.na(thisIntron$skippedExon) & thisIntron$type == "recomb"] ),
        novel_exon_with_intron_retention = mypaste(thisExon$name[!is.na(thisExon$retainedIntron) & thisExon$type == "recomb" ]),
        known_retained_intron = mypaste(na.omit(thisExon$retainedIntron[thisExon$type=="known"])),
        novel_retained_intron = mypaste(na.omit(thisExon$retainedIntron[thisExon$type=="recomb"])),
        known_skipped_exon = mypaste(na.omit(thisIntron$skippedExon[thisIntron$type=="known"])),
        novel_skipped_exon = mypaste(na.omit(thisIntron$skippedExon[thisIntron$type=="recomb"])),
        novel_comb = sum(all(thisExon$type == "known") & !novel_exon_skipping & !novel_utr_var),
        novel_comb_strict = sum(all(thisIntron$type == "known") & all(thisExon$type == "known")),
        n_exon_recomb = sum(thisExon$type == "recomb"),
        n_intron_recomb = sum(thisIntron$type == "recomb")) %>%
            ## we simplify here
            mutate(type= if(novel_exon_allnew | novel_exon_combnew) {"NE"}
                else if (novel_exon_extension) {"ASS"}
            else if (novel_intron_retention) {"IR"}
            else if (novel_exon_skipping)  {"ES"}
            else if (novel_comb)  {"NC"}
            else if (novel_utr_var) {"ASS"} ### tx first last exon splice sites variation
            else  {"AMB" },
            all_type = my_class_code(
                c(
                    sum(c(novel_exon_allnew,  novel_exon_combnew)),
                    sum(c(novel_exon_extension, novel_utr_var)),
                    novel_intron_retention,
                    novel_exon_skipping
                )
            ),
            is_multi_novel = is_multi_novel(
                c(
                    sum(c(novel_exon_allnew,  novel_exon_combnew)),
                    sum(c(novel_exon_extension, novel_utr_var)),
                    novel_intron_retention,
                    novel_exon_skipping
                )
            )

        )
    return(infoTab)
}, mc.cores=ncpu)

clusterClass = do.call(rbind, rvList)
## for transcript with multiple novel events
clusterClass$type[clusterClass$is_multi_novel] = "MN"

stopifnot(all(clusterClass$clname == novelClusterTab$clname))

clusterClass = clusterClass %>% left_join(novelClusterTab %>% select(
    clname, associated_tx, associated_gene, igv_coord, align_issue, is_end_variation
), by=c("clname" = "clname"))

save( novelTx, clusterClass, exonAnnot,dataExonRanges, intronAnnot, dataIntronRanges,
    file=file.path(infolder, "novel_cluster_class.rda"))
write_tsv(clusterClass, file.path(infolder, "novel_cluster_class.txt"))
