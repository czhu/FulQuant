##----------------------------------------------------------------------------##
## matching to annotation
match_by_clname = function(query, subject) {
    ## query and subject both have to be multiexonis
    if(is.null(query$clname))
        query$clname = granges_to_clname(query)
    if(is.null(subject$clname))
        subject$clname = granges_to_clname(subject)
    tibble(
        query_i = seq_len(length(query)),
        subject_i = match(query$clname, subject$clname), type="exact_match") %>% filter(!is.na(subject_i))
}

match_to_ref = function(query,subject) {
    defaultAns = tibble(query_i=integer(), subject_i=integer(), type = character(),delta_size=integer(),nhit=integer())
    exactMatch = match_by_clname(query,subject) %>% mutate(delta_size=0, nhit=1)
    indexForFuzyMatch = setdiff(seq_len(length(query)), exactMatch$query_i)
    ans = exactMatch

    if(length(indexForFuzyMatch)){
        message("Fuzzy matching ", length(indexForFuzyMatch), " tx")
        ## NOTE the cuttoff maxDeltaSize, smallExonSize, deltaDist should be the same as used in tx filter
        diffShift= suppressWarnings(find_misaligned_ss_shifted(query[indexForFuzyMatch], subject,
            ignore.strand=FALSE,maxDeltaSize=20, maxShift = 50,simplify=TRUE,  bothEndsShifted=FALSE))
        diffSmallExon = suppressWarnings(find_misaligned_small_exon(query[indexForFuzyMatch], subject,
            ignore.strand = FALSE, smallExonSize = 60, deltaDist=20))

        if(is.null(diffShift)) diffShift = defaultAns
        if(is.null(diffSmallExon)) diffSmallExon = defaultAns

        diffShift$query_i = indexForFuzyMatch[diffShift$query_i]
        diffSmallExon$query_i = indexForFuzyMatch[diffSmallExon$query_i]

        diffShift = diffShift %>% mutate(type="fuzzy_match_as") %>% select(query_i, subject_i, type,delta_size)
        diffSmallExon = diffSmallExon %>% mutate(type="fuzzy_match_se") %>% select(query_i, subject_i, type,delta_size)
        fuzzyMatch = rbind(diffShift, diffSmallExon) %>% group_by(query_i) %>% summarise(
            subject_i = subject_i[which.min(delta_size)],
            type = type[which.min(delta_size)],
            delta_size=delta_size[which.min(delta_size)],
            nhit = sum(delta_size  == min(delta_size)))
        ans = rbind(ans, fuzzyMatch) %>% arrange(query_i)
    }
    return(ans %>% arrange(query_i))
}

##################################################################################################
## helper function for genomewide view plot
granges_to_title = function(x){
    paste(paste0("chr",as.character(seqnames(x)),":"),start(x),"-",end(x))}
##################################################################################################
read_clipping_profile_file = function(infile){
    ans = read_tsv(infile,col_names=FALSE, progress=FALSE,col_types=list(X2="c"))
    colnames(ans) = c("qname","rname","i","width","type")
    ans
}

read_insert_profile_near_ss_file = function(infile){
    ans = read_tsv(infile,col_names=FALSE, progress=FALSE,col_types=list(X2="c"))
    colnames(ans) = c("qname","rname","i","intron_start","intron_end","width_left","width_right","width")
    ans
}

is_clustered = function(x) {
    ## infereed from cluster name
    !(x %in% c("singleton","unclustered"))
}

is_human = function(x){
    as.logical(seqnames(x) %in% c(1:22,"X","Y"))
}

extract_seqname_from_clname = function(clname){
    ## example: 10_1_88935367_15
    sub("^(.*?)_.*","\\1",clname)
}

extract_pos_from_clname = function(x){
    ## example: 10_1_88935367_15 or 10_88935367
    nfields = length(unlist(strsplit(x[1],"_")))
    if(nfields==4){
        ## name format as output by Wave's
        return(sub("^(.*?)_(.*?)_(.*?)_(.*?)$","\\3",x))
    } else if (nfields==2){
        return(sub("^(.*?)_(.*?)$","\\2",x))
    } else {
        stop("Number of fields isn't right!")
    }
}

extract_sampleName_from_qname = function(x){
    sub("^(.*?):.*$","\\1", x)
}

################################################################################
## helper function for exon connecvition analysi
is_overlapping_other_feature = function(x){
    lengths(as.list(findOverlaps(x,x, ignore.strand=FALSE)))>1
}

without_utr = function(x){
    y = blocks(x)
    mystart = unlist(end(y)[as(rep(1, length(x)), "IntegerList")]) + 1L
    myend = unlist(start(y)[as(lengths(y), "IntegerList")]) - 1L

    GRanges(seqnames(x), IRanges(mystart, myend), strand(x))
}

extract_internal_exons = function(x){
    if(!all(is_txClass(x)) | any(nblocks(x)<3) )
        stop("x must all be multi exonic (>3 exons) tx with internal exons")
    exonIV_to_intronIV( exonIV_to_intronIV(x) )
}


### we define a few terms
### clname: cluster name to represent splites sites for a given tx, splites sites is defined as the last/first bit of exon
### ivname: to represent introns
### clname and ivname share the same format chr_pos-chr_pos;chr_pos-chr_pos;chr_pos-chr_pos
### only difference is clname is the first/last bit of exon and can be only used for introns
### clname is a special case of ivname, ivname is more universal
granges_to_ivname = function(x){
    ## clname is standised name
    ## chrIS_77857-chrIS_95198;chrIS_77857_10-chrIS_95198
    if(!is_txClass(x)){stop("x has to be txClass")}
    iv = unlist(blocks(x))
    mystring = paste0(
        as.character(seqnames(iv)), "_", start(iv),"-",
        as.character(seqnames(iv)), "_", end(iv)
    )
    ans = split(mystring,rep(seq_len(length(x)), lengths(blocks(x))))
    ans = sapply(ans, paste, collapse=";")
    names(ans) = names(x)
    return(ans)
}

granges_to_clname = function(x){
    if ( !all(is_multiexonic(x)) ) { stop("x has to be multiexonic") }
    myintrons = exonIV_to_intronIV(x)
    ### +1 bp because the define ei and ei boudanries as the last bp of exon
    myiv = asBED(relist( unlist(blocks(myintrons)) +1, blocks(myintrons)))
    ans = granges_to_ivname(myiv)
    ans
}

granges_to_igvCoord = function(x){
    ## x is granges object
    ## return igv coordnates for quick copy and past
    paste0(as.character(seqnames(x)), ":",start(x),"-",end(x))
}
granges_to_igv_coord = granges_to_igvCoord

igvCoord_to_granges = function(x){
    ## x is "chr13:113,045,818-113,842,282"
    ## remove comma and space in case of any
    ans = gsub(",", "" , x)
    ans = gsub(" ", "" , ans)
    ## deal with long or short dash
    myregexp = "^(.+?):(.+)[-−](.+)$"
    GRanges(
        sub("chr","",sub(myregexp,"\\1",ans)),
        IRanges(as.integer(sub(myregexp,"\\2",ans)), as.integer(sub(myregexp,"\\3",ans)))
    )
}

ivname_to_granges = function(x){
    ## example chrIS_1800_77857_10-chrIS_850_95198_20
    ## this is vectoriesd version, the output is ranges for introns
    ## does NOT include nay bit of exon, ie.not incluing the junctions
    splt = strsplit(x,";")
    ## exon include UTR
    nintrons = lengths(splt) ## cluster name represents introns hence +1 for the number of exons

    splt2 = unlist(strsplit(unlist(splt),"-"))

    ssSitesChr = extract_seqname_from_clname(x)
    ssSites = as.integer(extract_pos_from_clname(splt2))

    ## we fill the starts and ends into the vector
    allstarts = ssSites[seq(1,length(ssSites),by=2)]
    allends = ssSites[seq(2,length(ssSites),by=2)]

    rv = GRanges(rep(ssSitesChr,nintrons),
        IRanges(
        allstarts,
        allends))

    rv = asBED(relist(rv, splt))
    names(rv) = x
    rv
}

## this should be the same as ivnames_to_granges
## works for ivname with strand
## FIXME minus strand not working
ivnames_to_granges2 = function(x){
        ## x is "2_74483410_+-2_74490244_+;2_74483410_+-2_74490244_+"
        mysplt = unlist(strsplit(x,";"))
        fmt = "^(.+)_(.+)_(.)-(.+)_(.+)_(.)$"
        ans = GRanges(
            as.character(sub(fmt,"\\1",mysplt)),
            IRanges(as.integer(sub(fmt,"\\2",mysplt)), as.integer(sub(fmt,"\\5",mysplt))),
            sub("=","-",as.character(sub(fmt,"\\3",mysplt))) )
        return(ans)
}

clname_to_granges = function(clname){
    ## this does not include UTR start and end!!
    ## for that you need read cluster please refer to
    ## infer_consensus function
    myranges = ivname_to_granges(clname)

    ans = asBED(relist( unlist(blocks(myranges)) - 1, blocks(myranges) ) )
    names(ans) = clname
    ans
}

ivname_to_granges_w_strand = function(clname){
    ## example chrIS_1800_+-chrIS_9000_+
    ## this is vectoriesd version, the output is ranges for introns
    ## does NOT include nay bit of exon, ie.not incluing the junctions
    splt = strsplit(clname,";")
    ## exon include UTR
    nintrons = lengths(splt) ## cluster name represents introns hence +1 for the number of exons

    spf = rep(1:length(splt), nintrons)
    splt2 = unlist(strsplit(unlist(splt),"-"))

    ssSitesChr = sub("^(.*?)_.*","\\1",clname)
    ssSites = as.integer(sub("^(.*?)_(.*?)_(.*?)$","\\2",splt2))

    ## we fill the starts and ends into the vector
    allstarts = ssSites[seq(1,length(ssSites),by=2)]
    allends = ssSites[seq(2,length(ssSites),by=2)]

    rv = GRanges(rep(ssSitesChr,nintrons),
        IRanges(
        ## !!!no +1/-1 here !!!!
        allstarts,
        allends))

    rv = asBED(split(rv,spf))
    rv$name = clname
    rv
}

has_left_end = function(x){
    ## x is qname
    grepl("^.*?:P-.-.:.*$",x) | grepl("^.*?:.-T-.:.*$",x)
}

has_right_end = function(x){
    grepl("^.*?:.-.-P:.*$",x) | grepl("^.*?:.-A-.:.*$",x)
}


infer_consensus = function(x) {
    ## x is GRangesList with cluster name as the element name
    stopifnot(is(x,"GRangesList"))

    ## FIXME is this the most accurate estimate?
    ## NOTE this is the bits how we infer the UTR start/end
    ## we estimate the start and end sites by taking median of all reads
    ## start and end are inferred using full length tx if available
    hasLeftEnd = relist(has_left_end(unlist(x)$name), x)
    hasRightEnd = relist(has_right_end(unlist(x)$name), x)

    ### for those without any adapter in the cluster we use them all
    hasLeftEnd[!any(hasLeftEnd)] = !hasLeftEnd[!any(hasLeftEnd)]
    hasRightEnd[!any(hasRightEnd)] = !hasRightEnd[!any(hasRightEnd)]

    mystart = round(median(start(x)[hasLeftEnd] ))
    myend  = round(median(end(x)[hasRightEnd]))

    myclusterName = names(x)
    ## BEGIN FIXME the part should be wrapped up as a function to replace clname_to_granges
    ## the problem: we need to esimates of the start and end of the trranscript
    ## this can only be done when we infer the consensus

    splt = strsplit(myclusterName,";")

    ## exon include UTR
    nexons = lengths(splt)+1 ## cluster name represents introns hence +1 for the number of exons

    spf = rep(1:length(splt),nexons)
    splt2 = unlist(strsplit(unlist(splt),"-"))

    ssSitesChr = extract_seqname_from_clname(myclusterName)
    ssSites = as.integer(extract_pos_from_clname(splt2))

    allstarts = allends = vector(mode="integer",length=sum(nexons))
    idx = cumsum(nexons)

    ## we fill the starts and ends into the vector
    allstarts[c(1,(idx+1)[-length(idx)])] = mystart
    allstarts[-c(1,(idx+1)[-length(idx)])] = ssSites[seq(2,length(ssSites),by=2)]

    allends[idx] = myend
    allends[-idx] = ssSites[seq(1,length(ssSites),by=2)]
    ## END

    ###  occcasionally we have end>start
    ###  this is artefact caused by exon smaller than the clustering winwdow, i.e exon less than 10 or 8 bp
    ###  quick fix we take the median of the start and end -> one bp exon magic
    ###  we sign coensesus_type "dubious" to these cases and should be excluded from the downstream analysis
    isArteFact = allstarts > allends
    if(any(isArteFact)){
        warning("check\t",paste(unique(spf[isArteFact]),collapse=" "))
        allstarts[isArteFact] = allends[isArteFact] =
            round(rowMeans(cbind(allstarts,allends)[isArteFact,,drop=FALSE]))
    }

    rv = GRanges(rep(ssSitesChr,lengths(splt)+1),
        IRanges(
        allstarts,
        allends))
    hasArteFact = any(as(split(isArteFact, spf), "LogicalList"))
    rv = asBED(split(rv,spf))
    rv$name = myclusterName
    rv$consensus_type = ifelse(hasArteFact,"dubious","putative")
    rv
}

overlap_by_perc = function(query,subject,ignore.strand=FALSE, simplify=TRUE){
    ## both query and subjct are grangesList, need to consider exon
    ## choose which covers the the query most, not considering mutual ovlp
    if(is(query,"GRangesList") & is(subject,"GRangesList")){
        mywidth = function(x) sum(width(x))
        myintersect = GenomicRanges::intersect
    } else if (is(query,"GRanges") & is(subject,"GRanges")){
        mywidth = function(x) width(x)
        myintersect = GenomicRanges::pintersect
    } else {
        stop("Both query and subject can only be either GRanges or GRangesList")
    }

    ovlps = suppressWarnings(findOverlaps(query, subject, ignore.strand=ignore.strand))

    myintsct = suppressWarnings(myintersect(
        query[queryHits(ovlps)],
        subject[subjectHits(ovlps)],ignore.strand=ignore.strand))

    ## mutual overlap
    res = tibble(
        query_i = queryHits(ovlps),
        subject_i = subjectHits(ovlps),
        ovlp_width = mywidth(myintsct),
        query_width = mywidth(query[query_i]),
        subject_width = mywidth(subject[subject_i]),
        per_ovlp_query = ovlp_width / query_width,
        per_ovlp_subject =  ovlp_width / subject_width
        )

    if(simplify){
        ans = res %>% group_by(query_i) %>% summarise(
            #subject_i = tx_i[which.max(per_ovlp_qseq)]
            n_mutual = n(),
            ##wh = if(sum(is_hit)==1){which(is_hit)} else {which.max(0.5*per_ovlp_qseq + per_ovlp_tx)},
            ### NOTE: if only one hit, take the one maximises per_ovlp_query
            ## if multiple hits with the same per_ovlp_query, take the one that maximises per_ovlp_subject
            wh = if(n_mutual==1){which.max(per_ovlp_query)} else {
                is_max_size = per_ovlp_query == max(per_ovlp_query)
                if(sum(is_max_size)==1){
                    which.max(per_ovlp_query)
                } else {
                    which(is_max_size)[which.max(per_ovlp_subject[which(is_max_size)])]
                }
                },
            subject_i = subject_i[wh],
            ovlp_width = ovlp_width[wh],
            query_width = query_width[wh],
            subject_width = subject_width[wh],
            per_ovlp_query = per_ovlp_query[wh],
            per_ovlp_subject =  per_ovlp_subject[wh]
        )
    }

    ans
}

extract_qscore_from_qname =  function(x){
    ## x is th qname such as eac94a59-de38-4eed-9da0-12e316f91f92:P-N-P:60:459:400:498:Q6.96
    as.numeric(sub("^.*:Q(.*)$","\\1",x))
}
extract_qscore = extract_qscore_from_qname

extract_txStrand_from_qname = function(x){
        ## x is th qname such as eac94a59-de38-4eed-9da0-12e316f91f92:P-N-P:60:459:400:498:Q6.96
    sub("^(.+?):(.+?):(.+?):.*$","\\3",x)
}


# this version based on gene boundary not exon based much more primitivs
# overlap_by_perc = function(query,subject,perc=0.5,ignore.strand=TRUE){
#     ## both query and subjct are granges, need to consider for spliced data
#     ## perc for query
#     ovlps = findOverlaps(query, subject, ignore.strand=ignore.strand)
#
#     myintsct = GenomicRanges::pintersect(
#         query[queryHits(ovlps)],
#         subject[subjectHits(ovlps)],ignore.strand=ignore.strand)
#
#     ## mutual overlap
#     res = tibble(
#         qseq_i = queryHits(ovlps),
#         tx_i = subjectHits(ovlps),
#         ovlp_size = width(myintsct),
#         qseq_size = width(query[queryHits(ovlps)]),
#         tx_size = width(subject[subjectHits(ovlps)]),
#         per_ovlp_qseq = ovlp_size / qseq_size,
#         per_ovlp_tx =  ovlp_size/ tx_size,
#         is_valid = per_ovlp_qseq > perc
#         )
#     res = res[res$is_valid,]
#
#     ans = res %>% group_by(qseq_i) %>% summarise(
#         #subject_i = tx_i[which.max(per_ovlp_qseq)]
#
#         n_mutual = n(),
#         ##wh = if(sum(is_hit)==1){which(is_hit)} else {which.max(0.5*per_ovlp_qseq + per_ovlp_tx)},
#         ### FIXME
#         wh = if(n_mutual==1){which.max(per_ovlp_qseq)} else {
#             is_max_size = per_ovlp_qseq == max(per_ovlp_qseq)
#             if(sum(is_max_size)==1){
#                 which.max(per_ovlp_qseq)
#             } else {
#                 which(is_max_size)[which.max(per_ovlp_tx[which(is_max_size)])]
#             }
#             },
#         subject_i = tx_i[wh]
#     )
#
#     ans
# }

fill_na_with_value = function(x,value=0){
    ## change NA to a defined value
    ifelse(is.na(x), value, x)
}

is_fulllength = function(x){
    ## x qname
    grepl("P-A-.",x) | grepl(".-T-P",x)
}

### copied from the plotting function
swap_strand = function(x){
    stopifnot(all(x %in% c("+","-")))
    ifelse(x=="+","-","+")
}

extract_polyAType_from_qname = function(x) {
    ## x is a standard qname like
    ## d98e37df-d929-4391-8a89-55248982993e:P-N-P:124:2648:2525:2760:Q11.41
    sub("^.*?:.-(.)-.:.*$","\\1",x)
}
extract_polyA_type = extract_polyAType_from_qname

# infer_strand_stat = function(polyAType,mappedStrand) {
#     ## tabluate infer strand
#     inferredStrand = ifelse(
#         polyAType=="N","*", ifelse(polyAType=="A",mappedStrand,swap_strand(mappedStrand))
#     )
#     table(factor(inferredStrand,c("+","-","*")))
# }
#
# infer_strand_pl = function(x, minPercDiff=0.3, ncpu=1){
#     ## GRangesList with each element  with name, include transcripts belong to the same cluster
#     ## pl stands for paralleled
#     ans = mclapply( 1:length(x), function(i) {
#         gr = x[[i]]
#         infer_strand_stat(extract_polyA_type(gr$name), as.character(strand(gr)))
#     }, mc.cores=ncpu)
#     ans = do.call(rbind, ans)
#
#     mydiffPerc = abs((ans[,"+"] - ans[,"-"]))/rowSums(ans[,c("+","-")])
#     rv = ifelse(rowSums(ans[,c("+","-")])==0,"*", ifelse(mydiffPerc>=minPercDiff, ifelse(ans[,"+"]>ans[,"-"],"+","-"),"*"))
#     return(rv)
# }

infer_strand_clusteredReads = function(x, minPercDiff=0.3,method="polyAType"){
    ## GRangesList with each element  with name, include transcripts belong to the same cluster
    ## method 1: polyAType use the polyA tail to determine txStrand
    ## method 2: splicingMotif: use the txStrand inferred by minimap2 tag ts using splicing motif
    mappedStrand = suppressWarnings(as.character(unlist(strand(x))))

    txStrand = switch(method,
        polyAType = {
        mypolyAType = suppressWarnings(extract_polyAType_from_qname(unlist(x)$name))
        ifelse(mypolyAType=="N","*", ifelse(mypolyAType=="A","+","-"))
        }, splicingMotif = {
            extract_txStrand_from_qname(unlist(x)$name)
        },{stop("method must either be polyAtype or splicingMotif!")}
    )

    inferredStrand = ifelse(
        txStrand=="*","*", ifelse(txStrand=="+",mappedStrand,swap_strand(mappedStrand))
    )
    inferredStrand = relist(factor(inferredStrand,c("+","-","*")), x)
    ans = table(inferredStrand)
    rownames(ans) = NULL

    mydiffPerc = abs((ans[,"+"] - ans[,"-"]))/rowSums(ans[,c("+","-"),drop=FALSE])
    rv = ifelse(rowSums(ans[,c("+","-"),drop=FALSE])==0,"*", ifelse(mydiffPerc>=minPercDiff, ifelse(ans[,"+"]>ans[,"-"],"+","-"),"*"))
    return(rv)
}

infer_strand= function(polyAType,mappedStrand,minDiffTwoStrandPerc = 0.4) {
    ## both vector
    ## polyAType: "A","T","N"
    ## mappedStrand: +,-
    ## inference can only be made on knownPolyAType,
    ## minDiffTwoStrandPerc: minimum difference between the percentage of two strand
    ## for cluster with 0.2 minus strand and 0.6 plus we consider a plus strand
    ## however if we have 0.2 minus and 0.3 we consider it to ambiguous so "*" instead of plus
    if(all(polyAType=="N")) {
        return("*")
    } else {
        wh = polyAType != "N"
        mystrd = ifelse(
            polyAType[wh]=="A", mappedStrand[wh], swap_strand(mappedStrand[wh]))
        tb = table(factor(mystrd,c("+","-")))

        myplusPerc = tb["+"]/length(mystrd)
        myminusPerc = tb["-"]/length(mystrd)

        ## cannot resolve if one strand is not more than the other swap_strand
        ## in ideal world, all tx should have only one
        if(abs(myplusPerc - myminusPerc)<minDiffTwoStrandPerc){
            return("*")
        } else {
            if(myplusPerc>0.5) {
                return("+")
            } else {
                return("-")
            }
        }
    }
}

infer_strand_cluster = function(x){
    ## x GRanges object with name, include transcripts belong to the same cluster
    infer_strand(extract_polyA_type(x$name),as.character(strand(x)))
}

read_polyA_file = function(x){
    read_tsv(polyAMaskFile,col_types=list(chr="c",start="i",end="i",perc="d",type="c"),progress=FALSE)
}

dist_to_polyARich_region = function(x,type="start",ignore.strand=TRUE,
    polyAMaskFile="/g/steinmetz/czhu/project/pore/data/genome/GRCh38.91/polyA_site_filter/polyA_mask_seed5_windowSize20_cutoff60_recovered.txt.gz"
    ){
    pos_to_gr = function(pos){ GRanges(seqnames(x),IRanges(pos,pos)) }
    polyATab = read_polyA_file(polyAMaskFile)
    mypolyGR = GRanges(polyATab$chr,IRanges(polyATab$start,polyATab$end))

    if(type=="start") {
        posgr = pos_to_gr(start(x))
    } else{
        posgr = pos_to_gr(end(x))
    }

    ans = distanceToNearest(posgr, mypolyGR,ignore.strand=ignore.strand)
    mcols(ans)$distance
}

dist_to_polyA_rich_region_strand_specific = function(
    x,
    polyAMaskFile="/g/steinmetz/czhu/project/pore/data/genome/GRCh38.91/polyA_site_filter/polyA_mask_seed5_windowSize20_cutoff60_recovered.txt.gz"
    ){
    ## x is tx class
    polyATab = read_polyA_file(polyAMaskFile)

    mypolyGR = GRanges( polyATab$chr, IRanges(polyATab$start,polyATab$end),
        ifelse(polyATab$type=="A","+","-") )
    poi = threeend(x)
    ans = rep(as.integer(NA),length=length(poi))

    distanceToNearestMC = function(x,y,ignore.strand=FALSE,ncpu=2,...){
        x$index = 1:length(x)

        rv = mclapply(split(x, cut(seq_len(length(x)),ncpu)),function(z){
            rv = distanceToNearest(z, y,ignore.strand=ignore.strand)
            tibble(
                queryHits=z$index[queryHits(rv)], subjectHits=subjectHits(rv),
                distance=mcols(rv)$distance )
        },mc.cores=ncpu)
        do.call(rbind,rv)
    }

    thisDist = distanceToNearestMC(poi, mypolyGR,ignore.strand = FALSE, ncpu=30)
    #ans = distanceToNearest(poi, mypolyGR,ignore.strand=FALSE)
    #ans = mcols(ans)$distance
    ans[thisDist$queryHits] = thisDist$distance

    if(any(strand(x) == "*")) {
        wh = which(strand(x) == "*")
        distLeft = distanceToNearestMC(
            GRanges( seqnames(x)[wh], IRanges(start(x)[wh],width=1) ),
            mypolyGR,ignore.strand = TRUE,ncpu=30)

        distRight = distanceToNearestMC(
            GRanges(seqnames(x)[wh], IRanges(end(x)[wh],width=1) ),
            mypolyGR, ignore.strand = TRUE, ncpu=30)

        mymat = matrix(rep(as.integer(NA) , 2*length(wh)),ncol=2)

        mymat[distLeft$queryHits,1] = distLeft$distance
        mymat[distRight$queryHits,2] = distRight$distance

        ans[wh] = matrixStats::rowMins(mymat,na.rm=TRUE)
    }
    return(ans)
}


extract_reads_by_qname = function(qnames,bamfile,outfile,outfile2,doRun=TRUE,wait=TRUE){
    stopifnot(file.exists(bamfile)) ### path to the bam file not valid
    ### outfile for included reads
    ### outfile2 for excluded reads if defined

    myfile = tempfile()
    writeLines(qnames,myfile)
    #picard FilterSamReads I=$infile O=reads_w_problem.bam READ_LIST_FILE=problem_read_names.txt FILTER=includeReadList
    ans = paste("picard FilterSamReads",
        paste0("I=", bamfile),
        paste0("O=", outfile),
        paste0("READ_LIST_FILE=", myfile),
        "FILTER=includeReadList VALIDATION_STRINGENCY=LENIENT",
        paste("&& samtools index",outfile))
    if(doRun){
        system(ans,wait=wait)
    }

    if(!missing(outfile2)){
        ## output good reads
        ##picard FilterSamReads I=$infile O=reads_wo_problem.bam READ_LIST_FILE=problem_read_names.txt FILTER=excludeReadList
        cmd2 = paste("picard FilterSamReads",
            paste0("I=", bamfile),
            paste0("O=", outfile2),
            paste0("READ_LIST_FILE=", myfile),
            "FILTER=excludeReadList VALIDATION_STRINGENCY=LENIENT",
            paste("&& samtools index",outfile2))
        ans = c(ans,cmd2)
        if(doRun) system(cmd2,wait=wait)
    }
    return(ans)
}

# myquery = GRangesList(GRanges("X",IRanges(c(1,20),c(50,70))))
# mysubject = GRangesList(GRanges("X",IRanges(c(1,20),c(50,70))))

# dist_between_gr = function(query,subject,ignore.strand=TRUE, simplify=TRUE){
#     require(GenomicRanges)
#     ## both query and subjct are grangesList, need to consider exon
#     ## choose which covers the the query most, not considering mutual ovlp
#     stopifnot(is(query,"GRangesList") & is(subject,"GRangesList"))
#     stopifnot(length(query) == length(subject))
#
#     ovlps = findOverlaps(query, subject, ignore.strand=ignore.strand)
#
#     if(length(ovlps)) {
#         myintsct = GenomicRanges::intersect(
#             query[queryHits(ovlps)],
#             subject[subjectHits(ovlps)],ignore.strand=ignore.strand)
#
#         myunion = GenomicRanges::union(
#                 query[queryHits(ovlps)],
#                 subject[subjectHits(ovlps)],ignore.strand=ignore.strand)
#
#         ## mutual overlap
#         res = tibble(
#             qseq_i = queryHits(ovlps),
#             tx_i = subjectHits(ovlps),
#             ovlp_size = sum(width(myintsct)),
#             union_size = sum(width(myunion)),
#             qseq_size = sum(width(query[queryHits(ovlps)])),
#             tx_size = sum(width(subject[subjectHits(ovlps)])),
#             diff_size = union_size - ovlp_size)
#         if(simplify){
#             ans = res %>% group_by(qseq_i) %>% summarise(
#                 #subject_i = tx_i[which.max(per_ovlp_qseq)]
#
#                 n = n(),
#                 ##wh = if(sum(is_hit)==1){which(is_hit)} else {which.max(0.5*per_ovlp_qseq + per_ovlp_tx)},
#                 wh = if(n==1){1} else {
#                     is_min_diff_size = diff_size == min(diff_size)
#                     if(sum(is_min_diff_size)==1){
#                         which.min(diff_size)
#                     } else {
#                         which(is_min_diff_size)[which.max(ovlp_size[which(is_min_diff_size)])]
#                         }
#                     },
#                 subject_i = tx_i[wh],
#                 diff_size = diff_size[wh]
#             )
#             return(ans$diff_size)
#         } else {
#             return(res$diff_size)
#         }
#     } else { return(Inf)}
# }

dist_between_gr = function(query,subject,ignore.strand=TRUE){
    require(GenomicRanges)
    ## both query and subjct are grangesList, need to consider exon
    ## choose which covers the the query most, not considering mutual ovlp
    stopifnot(is(query,"GRangesList") & is(subject,"GRangesList"))
    stopifnot(length(query) == length(subject))

    #isOvlp = overlapsAny(query,subject, ignore.strand=ignore.strand)
    #ans = rep(as.integer(NA), length(query))

    myintsct = GenomicRanges::intersect(
        query,
        subject,ignore.strand=ignore.strand)

    myunion = GenomicRanges::union(
            query,
            subject,ignore.strand=ignore.strand)

        ## mutual overlap
    ans = sum(width(myunion)) - sum(width(myintsct))
    ## no overlaps
    ans[sum(width(myintsct))==0] = Inf
    return(ans)
}


to_grangs_w_blocks = function(x){
    ## x is GRanges belonging to the same Tx
    asBED(GRangesList(x))
}


exonIV_to_intronIV = function(x){
    ## x is Granges with blokcs to exons
    ## this is a vectoriesd implementation
    ## this should give the same results as clname_to_iv
    if(!all(is_multiexonic(x)))
        stop("exonIV_to_intron_IV only works with multi exonic tx")
    tmpVal = blocks(x)
    nexons = lengths(tmpVal)
    allstarts = unlist(start(tmpVal))
    idxShift = cumsum(nexons)
    allstarts = allstarts[-c(1, (idxShift+1)[-length(idxShift)] )]
    allends = unlist(end(tmpVal))
    allends = allends[ -idxShift]
    myfac = rep(1:length(x), nexons-1)
    ans = GRanges(
        rep(as.character(seqnames(x)), nexons-1),
        IRanges(start=allends+1, end=allstarts-1),
        strand=rep(as.character(strand(x)), nexons-1)
    )
    ans = asBED(split(ans, factor(myfac,1:length(x))))
    if( !is.null(names(blocks(x))) ){
        names(ans) = names(blocks(x))
    }
    return(ans)
}
## for backward compatibility
exonIV_to_intron_IV = exonIV_to_intronIV

dist_between_tx_intron_size = function(tx1, tx2, ignore.strand=TRUE) {
    dist_between_gr(blocks(exonIV_to_intron_IV(tx1)),blocks(exonIV_to_intron_IV(tx2)),ignore.strand = ignore.strand)
}

dist_between_tx_exon_size = function(tx1, tx2, ignore.strand=TRUE) {
    dist_between_gr(blocks(tx1),blocks(tx2), ignore.strand = ignore.strand)
}

is_different_at_five_end  = function(shortTx,longTx){
    ## query and subject are are GRAnges object with blocks, having strand!!!
    ## representing exons
    if(as.character(strand(shortTx)) != as.character(strand(longTx)) ){
        return(FALSE)
        #stop("Both short tx and long tx must be on the same strand\n")
    }
    if(width(shortTx)>width(longTx)){
        warning("Short tx should be shorter than long tx!!\n")
        return(FALSE)
    }
    mystrand = as.character(strand(shortTx))

    ## x and y are both GRanges object with blocks, only representing introns
    is_fully_within = function(x,y) {
        ## is x within y
        ans = all(blocks(x) %in% blocks(y)) & (nblocks(x) < nblocks(y))
                ## rule out intron retention
        if(any(ans)){
            sel = which(ans)
            isIntronRetained = any(diff(which( suppressWarnings(blocks(y)[sel] %in% blocks(x)[sel])))>1)
            ans[sel][isIntronRetained] = FALSE
        }
        return(ans)
    }
    shortTxIntrons = exonIV_to_intron_IV(shortTx)
    longTxIntrons = exonIV_to_intron_IV(longTx)

    if( is_fully_within(shortTxIntrons,longTxIntrons) ){
        if(mystrand=="+"){
            return( (start(shortTxIntrons)>start(longTxIntrons)) & ( end(shortTxIntrons)==end(longTxIntrons) ))
        } else {
            return( (end(shortTxIntrons)<end(longTxIntrons)) & (start(shortTxIntrons) == start(longTxIntrons)) )
        }
    } else {
        return(FALSE)
    }
    FALSE
}

extract_filter_flag_from_clname = function(x){
    ## x is cluster name for displaying in IGV
     as.integer(sub("^F(\\d+):.*","\\1",x))
}

standardise_clname = function(x) {
    ## remove 4840 and 4 from 7_4840_127588566_4
    ## so it only contains  chr_pos, the cluster name looks then like
    ## chr_pos-chr_pos;chr_pos-chr_pos
    ## 4840 is index (count high to) needed for the clustering pipeleine from Wave
    ## 4 is the count number
    ## not comparable between the runs
    ans = gsub("(.*?)_.*?_(.*?)_.*?-(.*?)_.*?_(.*?)_.*?;","\\1_\\2-\\3_\\4;",paste0(x,";"))
    sub( ";$", "" , ans )
}

calc_bitcode_flag = function(x){
    ## turn a matrix of filter into bit code and summarise them
    ## similar to samtools bit flag
    bitcodePerCol = 2^(1:ncol(x)-1)
    mat = t(x)
    ## mat[is.na(mat)]=0
    colSums(mat*bitcodePerCol,na.rm=TRUE)
}

which2keep = function(x){
    ## x is QC matrix or data frame with TRUE/FALSE or NA
    ## returns row index passing all qc filters (TRUE means pass)
    which(rowSums(x,na.rm=TRUE) == ncol(x))
}

adaptive_filter = function(x, cutoff = 0.01, ncpu=1){
    ## x is a cluster filter object with slots consensus, filter
    ## this adds an additional splot pass_adaptive or overrites the existing
    ## in the filter matrix
    ## not the best code , lots of room for improvement
    idxKeep = which2keep(x$filter)
    rv = rep(as.logical(NA), nrow(x$filter))
    rv[idxKeep] = FALSE

    ovlps = findOverlaps(x$consensus[idxKeep], x$tu,type="within")

    mymat = x$stat
    mymat$index = 1:nrow(mymat)
    mymat = mymat[idxKeep,]

    ## make sure each entry in mymat correspondds to the one tu
    stopifnot(nrow(mymat) == length(ovlps))

    ans = mclapply(split(mymat, subjectHits(ovlps)), function(thisTu){
        thisTu$index[which(thisTu$count >= (max(thisTu$count)*cutoff))]
    },mc.cores=ncpu)
    ans = unlist(ans)
    rv[ans] = TRUE
    x$filter$pass_adaptive = rv
    x
}
#
# advanced_filter1 = function(x, minDist=50, ncpu=1){
#     ## x is a cluster filter object with slots consensus, filter
#     ## this adds an additional splot pass_adaptive or overrites the existing
#     ## in the filter matrix
#     ## this filter remove transcript that is very close to other transcripts in the same transricpt unit
#     ## if the distance is shorter by given cutoff minDist. At this stage, we do not have a relative count
#     ## threshold.
#     ## returns cluster filer object with addtional column pass_sigdiff in the filter matrix
#     idxKeep = which2keep(x$filter)
#     rv = rep(as.logical(NA), nrow(x$filter))
#     rv[idxKeep] = FALSE
#
#     mymat = x$stat
#     mymat$index = 1:nrow(mymat)
#     mymat = mymat[idxKeep,]
#
#     ovlps = findOverlaps(x$consensus[idxKeep], x$tu, type="within")
#     message(length(subjectHits(ovlps))," total clusters to filter")
#
#     stopifnot(nrow(mymat) == length(subjectHits(ovlps)))
#
#     dat = split(mymat, subjectHits(ovlps))
#     ans = mclapply( (1:length(dat)), function(ith) {
#         thisTu=dat[[ith]]
#         message( "process ", ith," th Cluster" )
#
#         y = thisTu[order(thisTu$count,decreasing=TRUE),]
#         if(nrow(y)>1){
#             shouldKeep = rep(TRUE, nrow(y))
#             for(i in 2:nrow(y)){
#                 thistx = x$consensus[y$index[i]]
#                 distToOtherTx = sapply(
#                     (1:(i-1))[shouldKeep[1:(i-1)]], function(j){
#                     reftx = x$consensus[y$index[j]]
#                     dist_between_tx_intron_size(thistx, reftx)}
#                 )
#                 if(any(distToOtherTx<minDist)){
#                     shouldKeep[i]=FALSE
#                 }
#             }
#             return(y$index[shouldKeep])
#         } else {
#             return(y$index)
#         }
#     }, mc.cores=ncpu)
#     ans = unlist(ans)
#     rv[ans]=TRUE
#     x$filter$pass_sigdiff= rv
#     x
# }

# advanced_filter1 = function(x, minDist=20, ncpu=1,smallExonSize=60,deltaDist=30){
#     ## x is a cluster filter object with slots consensus, filter
#     ## this adds an additional splot pass_adaptive or overrites the existing
#     ## in the filter matrix
#     ## this filter remove transcript that is very close to other transcripts in the same transricpt unit
#     ## if the distance is shorter by given cutoff minDist. At this stage, we do not have a relative count
#     ## threshold.
#     ## returns cluster filer object with addtional column pass_sigdiff in the filter matrix
#     ## there're two types of artefacts we are dealing with
#     ## share the all exons but only different in 1 or 2 splice sites by a small amount
#     ## exon is missing because it's too small this results in adjacent splice sites shifted
#     ## we merge on the fly
#     idxKeep = which2keep(x$filter)
#     isPassed = rep(as.logical(NA), nrow(x$filter))
#     isPassed[idxKeep] = FALSE
#
#     mymat = x$stat
#     mymat$index = 1:nrow(mymat)
#
#     mymat = mymat[idxKeep,]
#
#     ovlps = findOverlaps(x$consensus[idxKeep], x$tu, type="within")
#     message(length(unique(subjectHits(ovlps))) ," total clusters to filter")
#
#     mycons = x$consensus
#     mycons$index = 1:length(mycons)
#
#     stopifnot(nrow(mymat) == length(ovlps))
#
#     dat = split(mymat, subjectHits(ovlps))
#     idx = (1:length(dat))
#     ans = mclapply( idx, function(ith) {
#         thisTu=dat[[ith]]
#         message( "process ", ith," th Cluster" )
#
#         y = thisTu[order(thisTu$count,decreasing=TRUE),]
#         rv = tibble(index=y$index, index_reassigned = y$index, keep=rep(TRUE, nrow(y)),type=rep("None",nrow(y)))
#
#         if(nrow(y)>1){
#             #shouldKeep = rep(TRUE, nrow(y))
#             for (i in 2:nrow(y)) {
#                 thistx = mycons[y$index[i]]
#                 reftx = mycons[y$index[(1:(i-1))[rv$keep[1:(i-1)]]]]
#
#                 #diffShift = find_misaligned_ss_shifted_v2(thistx,reftx, ignore.strand=FALSE, maxDist=minDist)
#                 diffShift = find_misaligned_ss_shifted(thistx,reftx, ignore.strand=FALSE, maxDeltaSize=minDist,maxShift = 50, simplify=TRUE,
#                     bothEndsShifted=FALSE)
#                 diffSmallExon = find_misaligned_small_exon(thistx, reftx,ignore.strand=FALSE, smallExonSize = smallExonSize, deltaDist=deltaDist)
#
#                 if(is.null(diffShift) & is.null(diffSmallExon)) next
#                 hasShiftIssue = hasSmallExonIssue = FALSE
#
#                 ## NOTE after running the code another time I realised that we do not need to
#                 ## consider the case of self match here, because for tx cluster filter purpose
#                 ## this case would not happen
#                 if(!is.null(diffShift)){
#                     diffShift = diffShift %>% filter(delta_size>0)
#                     hasShiftIssue = nrow(diffShift)>0
#                 }
#
#                 if(!is.null(diffSmallExon)){
#                     diffSmallExon = diffSmallExon %>% filter(delta_size>0)
#                     hasSmallExonIssue = nrow(diffSmallExon)>0
#                 }
#
#                 #if(!is.null(diffShift) | !is.null(diffSmallExon)) rv$keep[i]=FALSE
#
#                 ## decide which is the closest to assign to
#                 if( hasShiftIssue & hasSmallExonIssue ) {
#                 ### should be rare case
#                     minDiffShift = min(diffShift$delta_size)
#                     minDiffSmallExon = min(diffSmallExon$delta_size)
#                     if(minDiffShift<=minDiffSmallExon) {
#                         rv$index_reassigned[i] = reftx$index[diffShift$subject_i[which.min(diffShift$delta_size)]]
#                         rv$type[i] = "ss_shift"
#                     } else {
#                         rv$index_reassigned[i] = reftx$index[diffSmallExon$subject_i[which.min(diffSmallExon$delta_size)]]
#                         rv$type[i] = "small_exon"
#                     }
#                     rv$keep[i]=FALSE
#                 } else if ( hasShiftIssue ){
#                     rv$index_reassigned[i] = reftx$index[diffShift$subject_i[which.min(diffShift$delta_size)]]
#                     rv$type[i] = "ss_shift"
#                     rv$keep[i]=FALSE
#                 } else if ( hasSmallExonIssue ){
#                     rv$index_reassigned[i] = reftx$index[diffSmallExon$subject_i[which.min(diffSmallExon$delta_size)]]
#                     rv$type[i] = "small_exon"
#                     rv$keep[i]=FALSE
#                 } else next
#             }
#         }
#         return(rv)
#     }, mc.cores=ncpu)
#     ans = do.call(rbind,ans)
#     isPassed[ans$index[ans$keep]] = TRUE
#     x$filter$pass_sigdiff = isPassed
#     x$tx_reassigned = ans[,c("index","index_reassigned","type")]
#     return(x)
# }

advanced_filter1 = function(x, minDist=20, ncpu=1,smallExonSize=60,deltaDist=30){
    ## x is a cluster filter object with slots consensus, filter
    ## this adds an additional splot pass_adaptive or overrites the existing
    ## in the filter matrix
    ## this filter remove transcript that is very close to other transcripts in the same transricpt unit
    ## if the distance is shorter by given cutoff minDist. At this stage, we do not have a relative count
    ## threshold.
    ## returns cluster filer object with addtional column pass_sigdiff in the filter matrix
    ## there're two types of artefacts we are dealing with
    ## share the all exons but only different in 1 or 2 splice sites by a small amount
    ## exon is missing because it's too small this results in adjacent splice sites shifted
    ## we merge on the fly
    idxKeep = which2keep(x$filter)
    isPassed = rep(as.logical(NA), nrow(x$filter))
    isPassed[idxKeep] = FALSE

    mymat = x$stat
    mycons = x$consensus
    mymat$index = mycons$index = 1:nrow(mymat)

    mymat = mymat[idxKeep,]
    mycons = mycons[idxKeep]

    ovlps = findOverlaps(mycons, x$tu, type="within")
    message(length(unique(subjectHits(ovlps))) ," total clusters to filter")

    stopifnot(nrow(mymat) == length(ovlps))

    dat = split(mymat, subjectHits(ovlps))
    idx = (1:length(dat))

    ## FIXME consider either diffShift diffSmallExon is null
    diffShift = find_misaligned_ss_shifted(mycons,mycons, ignore.strand=FALSE, maxDeltaSize=minDist,
        maxShift = 50, simplify=FALSE, bothEndsShifted=FALSE)  %>%
        filter(query_i!=subject_i) %>% mutate(query_i=mycons$index[query_i], subject_i = mycons$index[subject_i])
    diffSmallExon = find_misaligned_small_exon(mycons, mycons,ignore.strand=FALSE,
        smallExonSize = smallExonSize, deltaDist=deltaDist, simplify=FALSE) %>%
        filter(query_i!=subject_i) %>% mutate(query_i=mycons$index[query_i], subject_i = mycons$index[subject_i])

    ans = mclapply( idx, function(ith) {
        message( "process ", ith," th Cluster" )

        y = dat[[ith]] %>% arrange(desc(count)) %>% mutate(keep=TRUE)

        if(nrow(y)>1){
            #shouldKeep = rep(TRUE, nrow(y))
            for (i in 2:nrow(y)) {
                thistx = y$index[i]
                reftx = y$index[(1:(i-1))[y$keep[1:(i-1)]]]

                hasShiftIssue = nrow(diffShift %>% filter(query_i==thistx, subject_i %in% reftx))>0
                hasSmallExonIssue = nrow(diffSmallExon %>% filter(query_i==thistx, subject_i %in% reftx))>0
                if( hasShiftIssue  | hasSmallExonIssue ){
                    y$keep[i]=FALSE
                }
            }
        }
        return(y)
    }, mc.cores=ncpu)
    ans = do.call(rbind,ans)
    isPassed[ans$index] = ans$keep
    x$filter$pass_sigdiff = isPassed
    return(x)
}

nblocks = function(x){
    lengths(blocks(x))
}

is_multiexonic = function(x) {
    if(!is_txClass(x))
        stop("only works with txClass")
    nblocks(x) > 1
}

is_txClass = function(x){
    is(x, "GRanges") & (!is.null(blocks(x)))
}

is_exon_all_overlap = function(query,subject) {
    ### no strand consideration here
    ### do all exons in each query and subject pairwise overlap
    if( !all( nblocks(query) == nblocks(subject)) )
        stop("both query and subject should have same number of exons")

    mycoords = tibble(
        start_query = unlist(start(blocks(query))),
        end_query = unlist(end(blocks(query))),
        start_subject= unlist(start(blocks(subject))),
        end_subject = unlist(end(blocks(subject))))
    ans = (mycoords$end_query > mycoords$start_subject) & (mycoords$start_query < mycoords$end_subject)
    ans = all(relist(ans,blocks(query)))
    return(ans)
}


tx_width = function(x){
    if(!is_txClass(x)){
        stop("x must be a txClass object")
    }
    sum(width(blocks(x)))
}

represent_minus_strand = function(x){
    ## replace "-" with "=" to represent minus strand
    sub("-","=",x)
}

make_ivname_from_granges = function(x, includeStrand=FALSE){
    ## x is granges

    if(includeStrand) {
        mysuffix = paste0("_",represent_minus_strand(as.character(strand(x))))
    } else {
        mysuffix = ""
    }

    ans = paste0(
        as.character(seqnames(x)), "_", start(x), mysuffix,
        "-",
        as.character(seqnames(x)), "_", end(x), mysuffix)
    return(ans)
}

## development of find_misaligned_ss_shifted_v2 function
##ref         ----      ----      ----
##true pos    -----      ---      ----  mostly likely due to an alignment issue
##false pos   ----     -----       ---  both splice sites are novel and and should be retained
##ignore      -----      ---             ----
# library(rtracklayer)
#source("/Users/czhu/Dropbox/Developer/Eclipse/stuffs/nanopore/dev_tx_mapping/reads_clustering/clustering_functions.R")
# myref = asBED( GRangesList(GRanges("x",IRanges(start=c(1,20,40), end=c(10,30,50)))) )
# ## first one is , second is a false positive hit
# mydata = asBED( GRangesList(
#     ## the same
#     GRanges("x",IRanges(start=c(1,20,40), end=c(10,30,50))),
#     ## a true positive hit, shift of
#     GRanges("x",IRanges(start=c(1,24,40), end=c(14,30,50))),
#     GRanges("x",IRanges(start=c(1,16,44), end=c(10,30,50))),
#     ## negative case
#     GRanges("x",IRanges(start=c(1,24,60), end=c(14,30,70)))
# ))
# plot_feature_vpr(c(myref,mydata),featureCols=c(rep("firebrick",length(myref)),rep("steelblue",length(mydata))))
# tx_width(myref) == tx_width(mydata)
# ## function 1 compare to a reference annotation
# find_misaligned_ss_shifted_v2(mydata,myref,maxDeltaSize=1)

find_misaligned_ss_shifted_v2 = function(query, subject, ignore.strand=FALSE, maxDeltaSize=20,maxShift = 50,
    bothEndsShifted=FALSE, simplify=TRUE) {
    ## query is input
    ## subject is usualy reference to compare against
    ## anything that differ less than maxDist is considered the same (due to alignment error)
    ## maxShift distance to known splice site any allowed
    ## maxDiffSize: diff for one given intron
    ## simplify: reduce to one query to one subject mapping, if FALSE, one query could map multiple
    ## bothEndsShifted: TRUE require both end are at least shifted, FALSE one end shift is fine
    ## subject
    if(!is_txClass(query) | !(is_txClass(subject))) {
        stop("Both query and subject have to be txClass")}

    if(!all(is_multiexonic(query)) | !all(is_multiexonic(subject))) {
    ## FIXME singleton can be safely ignored without stopping the function
        stop("Both query and subject have to be multiexonic")}

    ## all exon overlap
    queryIntrons = exonIV_to_intronIV(query)
    subjectIntrons= exonIV_to_intronIV(subject)

    ovlps = findOverlaps(blocks(queryIntrons), blocks(subjectIntrons))
    if(length(ovlps) == 0) return(NULL)

    ## same number of exons
    hasSameNumExons = nblocks(queryIntrons[queryHits(ovlps)]) == nblocks(subjectIntrons[subjectHits(ovlps)])
    if(!any(hasSameNumExons)) return(NULL)
    ovlps = ovlps[hasSameNumExons]

    mycoords = tibble(
        start_query = unlist(start(blocks(queryIntrons[queryHits(ovlps)]))),
        end_query = unlist(end(blocks(queryIntrons[queryHits(ovlps)]))),
        start_subject= unlist(start(blocks(subjectIntrons[subjectHits(ovlps)]))),
        end_subject = unlist(end(blocks(subjectIntrons[subjectHits(ovlps)]))))

    sk = blocks(queryIntrons[queryHits(ovlps)])

    isExonsAllOverlap = (mycoords$end_query > mycoords$start_subject) & (mycoords$start_query < mycoords$end_subject)
    isExonsAllOverlap = all(relist(isExonsAllOverlap, sk))

    diff_size_start = relist(mycoords$start_query - mycoords$start_subject,sk)
    diff_size_end = relist(mycoords$end_query - mycoords$end_subject,sk)

    delta_size = abs(diff_size_start - diff_size_end)
    passMaxShift = all(abs(diff_size_start) <= maxShift) & all(abs(diff_size_end) <= maxShift)

    ## delta size < cutoff and shift
    belowCutoff = delta_size <= maxDeltaSize

    if(bothEndsShifted){
        hasBothEndsShifted = abs(diff_size_start)>0 & abs(diff_size_end)>0
        belowCutoff =  belowCutoff & hasBothEndsShifted
    }

    ans = tibble(
        query_i = queryHits(ovlps), subject_i = subjectHits(ovlps),
        delta_size =  sum(delta_size),
        n_diff_both_end = sum(belowCutoff & (abs(diff_size_start) >0  & abs(diff_size_end) >0)),
        n_diff_either_end = sum(belowCutoff & (abs(diff_size_start) >0  | abs(diff_size_end) >0)) - n_diff_both_end
    )[passMaxShift & all(belowCutoff),]

    if(simplify) {
        ans = ans %>% group_by(query_i) %>% summarise(
            subject_i = subject_i[which.min(delta_size)],
            n = n(),
            nhit =sum(delta_size==min(delta_size)),
            delta_size=delta_size[which.min(delta_size)],
            n_diff_both_end = n_diff_both_end[which.min(delta_size)],
            n_diff_either_end = n_diff_either_end[which.min(delta_size)])
    }

    return(ans)
}

## function 2 not tx annotation exists but only splice sites
## FIXME this needs to be implemented

## old implementation without considering if alignment issues belong to adjacent exons
find_misaligned_ss_shifted_v1 = function(query, subject, ignore.strand=FALSE, maxDist=20) {
    ## query is input
    ## subject is usualy reference to compare against
    ## anything that differ less than maxDist is considered the same (due to alignment error)
    if(!is_txClass(query) | !(is_txClass(subject))) {
        stop("Both query and subject have to be txClass")}

    ## all exon overlap
    mydist = distanceToNearestTx(query, subject, ignore.strand = ignore.strand, simplify=TRUE)
    if(is.null(mydist) |  all(!is.finite(mydist$delta_size)) | all(mydist$delta_size>maxDist)) return(NULL)

    ## same number of exons
    hasSameNumExons = lengths(blocks(subject[mydist$subject_i])) == lengths(blocks(query[mydist$query_i]))
    ans = mydist[hasSameNumExons,]
    if(nrow(ans)==0) return(NULL)
    ## exon all overlap
    isExonsAllOverlap = is_exon_all_overlap(query[ans$query_i] , subject[ans$subject_i])
    if(all(!isExonsAllOverlap)) return(NULL)

    ans = ans[isExonsAllOverlap, ] %>% filter(delta_size <= maxDist)
    return(ans)
}
find_misaligned_ss_shifted = find_misaligned_ss_shifted_v2

## test dataset
# mydata = asBED(GRangesList(
#     c(GRanges("x",IRanges(5,12)), GRanges("x",IRanges(40,50))), ## false pos
#     c(GRanges("x",IRanges(5,12)), GRanges("x",IRanges(20,25)), GRanges("x",IRanges(40,50))) ## true hit
# ))
#
# myref = asBED(GRangesList(
#     c(GRanges("x",IRanges(5,10)), GRanges("x",IRanges(20,25)), GRanges("x",IRanges(40,50)))
# ))
## alternative implentation consider using introns as the unit
find_misaligned_ss_shifted_any_using_ss = function(query, subject, ss_max_shift_size=20, simplify=TRUE){
    ## NOTE: I cannot control if both splice sites come from the same intron in reference
    ## due to the way the computation is done with splice sites sparately
    ## problematch cases, both splice sites are shifted,
    ## however one splice site happens to overlap a known refereence splice sites
    ## the nearest used will consider the wrong splice sites as the reference point
    ## because the distances is the shortest, this combination (intron) won't occur in
    ## reference introns, we will miss these cases
    ## find tx with any splice sites deviate from the known
    ## ss_max_shift_size max distance to a known splice site
    ## used here to filter for true novel splice sites

    ## we compare with splice sites extracted from reference
    ## dist_ei defines as ei_read - ei_ref
    ## dist_ie difine as ie_ref - ie_read
    ##  ei + sign
    ##  ref  ----
    ##  read -------
    ##  ei - sign
    ##  ref  ----
    ##  read --
    ##  ie  +
    ##  ref     ----
    ##  read  ------
    ##  ref     ----   - sign
    ##  read     --
    ## the sign indicates if the extension is longer or shorter than the reference in the respective
    ## oriention
    if( !is_txClass(query) | !is_txClass(subject) ) {
        stop("Both query and subject have to be txClass")
    }

    if( !all(is_multiexonic(query)) |  !all(is_multiexonic(subject))){
        stop("Both query and subject have to be multiexonic")
    }

    queryIntrons = exonIV_to_intronIV(query)
    subjectIntrons = exonIV_to_intronIV(subject)

    ## -1 EI is the last base of exon
    myextract_ss = function(x, fix){
        shift(resize( unlist(blocks(x)), 1L, fix=fix, ignore.strand=TRUE),ifelse(fix=="start",-1L,1L) )
    }

    queryEI = myextract_ss(queryIntrons, "start")
    queryIE = myextract_ss(queryIntrons, "end")
    queryEI$tx_i = queryIE$tx_i = rep(seq_len(length(queryIntrons)), nblocks(queryIntrons))
    queryEI$ith_intron = queryIE$ith_intron = unlist(sapply(nblocks(queryIntrons), seq_len))


    subjectEI = unique(myextract_ss(subjectIntrons, "start"))
    subjectIE = unique(myextract_ss(subjectIntrons, "end"))

    mydistanceToNearest = function(x,y){
        ## distanceToNearest gives a distance 0 for case pos 1 versus 2
        ## this should be distance of 1 but since it's 0 it doesn't differentite the following case
        ## pos 1 versus pos 1
        ## pos 1 versus pos 2
        if( any(width(x)>1) | any(width(y)) >1 ){
            stop("both x and y can only contain range with width of 1")
        }
        suppressWarnings( as_tibble( nearest(x, y, select="all") ) ) %>%
            mutate( distance=abs(start(x)[queryHits] - start(y)[subjectHits]) ) %>%
            #filter(distance>0) %>%
            group_by(queryHits) %>% summarise(
                subjectHits = subjectHits[which.min(distance)],
                distance = distance[which.min(distance)])
    }

    distEI= mydistanceToNearest(queryEI, subjectEI) %>%
        filter(distance < ss_max_shift_size) %>%
        mutate( signed_dist= start(queryEI)[queryHits] - start(subjectEI)[subjectHits] )

    distIE= mydistanceToNearest(queryIE, subjectIE) %>%
        filter(distance < ss_max_shift_size) %>%
        mutate( signed_dist = start(subjectIE)[subjectHits] - start(queryIE)[queryHits])

    ans = full_join(distEI,distIE,by=c("queryHits"="queryHits"), suffix = c(".ei", ".ie")) %>%
        ## only if both distance can be calculated
        filter(!is.na(subjectHits.ei), !is.na(subjectHits.ie)) %>%
        mutate(query_i=queryEI$tx_i[queryHits],intron_i=queryEI$ith_intron[queryHits],
            dist_ei=fill_na_with_value(signed_dist.ei,0L),
            dist_ie=fill_na_with_value(signed_dist.ie,0L))  %>%
            filter(dist_ei!=0 | dist_ie!=0, ) %>%
            mutate(
                delta_size=abs(dist_ei+dist_ie),
                ivname = paste0(
                    seqnames(subjectEI)[subjectHits.ei],"_",
                    start(subjectEI)[subjectHits.ei],"-",
                     seqnames(subjectIE)[subjectHits.ie],"_",
                     start(subjectIE)[subjectHits.ie]) ) %>%
            arrange(query_i) #%>% select(query_i,intron_i,dist_ei,dist_ie)

    if(simplify){
        mypaste = function(x) paste(x,collapse=",")
        ans = ans %>% group_by(query_i) %>% summarise (
            intron_i = mypaste(intron_i),
            n_diff = n(),
            dist_ei = mypaste(dist_ei),
            dist_ie = mypaste(dist_ie),
            delta_size = sum(delta_size)
        )
    }

    return(ans)
}


find_misaligned_ss_shifted_any_using_intron =  function(query, subject, ss_max_shift_size=20, simplify=TRUE){
    ## alternative implentation to find_misaligned_ss_shifted_any_using_ss
    ## this should address the problen above: control for if the ss combination
    ## did occur in reference as introns

    if( !is_txClass(query) | !is_txClass(subject) ) {
        stop("Both query and subject have to be txClass")
    }

    if( !all(is_multiexonic(query)) |  !all(is_multiexonic(subject))){
        stop("Both query and subject have to be multiexonic")
    }

    queryIntrons = exonIV_to_intronIV(query)
    subjectIntrons = exonIV_to_intronIV(subject)

    queryIntronsSep = unlist( blocks(queryIntrons) )
    queryIntronsSep$tx_i = rep( seq_len(length(queryIntrons)), nblocks(queryIntrons) )
    queryIntronsSep$ith_intron = unlist( sapply(nblocks(queryIntrons), seq_len) )

    subjectIntronsSep = unique(unlist( blocks(subjectIntrons) ))

    ans = suppressWarnings( as_tibble(findOverlaps( queryIntronsSep, subjectIntronsSep )) ) %>%
        mutate(
            dist_ei = start(queryIntronsSep)[queryHits] - start(subjectIntronsSep)[subjectHits] ,
            dist_ie = end(subjectIntronsSep)[subjectHits] - end(queryIntronsSep)[queryHits],
            query_i=queryIntronsSep$tx_i[queryHits],
            intron_i=queryIntronsSep$ith_intron[queryHits],
            delta_size = abs(dist_ei + dist_ie)) %>% group_by(queryHits) %>%
            dplyr::slice(which.min(delta_size)) %>%
            filter(dist_ei!=0 | dist_ie!=0, abs(dist_ei) < ss_max_shift_size,
                abs(dist_ie) < ss_max_shift_size) %>%
            arrange(queryHits)

    if(anyDuplicated(ans$queryHits))
        stop( "one intron can be only compared to one reference intron" )

    if(simplify){
        mypaste = function(x) paste(x,collapse=",")
        ans = ans %>% group_by(query_i) %>% summarise (
            intron_i = mypaste(intron_i),
            n_diff = n(),
            dist_ei = mypaste(dist_ei),
            dist_ie = mypaste(dist_ie),
            delta_size = sum(delta_size)
        )
    }

    return(ans)
}

find_misaligned_ss_shifted_any = find_misaligned_ss_shifted_any_using_intron

##ref         ----   --    ----      ----
##true pos    -----       -----      ----  mostly likely due to an alignment issue
##false pos   ----         -----    -----  both splice sites are novel and and should be retained

# library(rtracklayer)
#source("/Users/czhu/Dropbox/Developer/Eclipse/stuffs/nanopore/dev_tx_mapping/reads_clustering/clustering_functions.R")
# myref = asBED( GRangesList(GRanges("x",IRanges(start=c(1,20,30,50), end=c(10,25,40,60)))) )
# # ## first one is , second is a false positive hit
# mydata = asBED( GRangesList(
#     ## the same
#     GRanges("x",IRanges(start=c(1,20,30,50), end=c(10,25,40,60))),
#     ## a true positive hit, shift of
#     GRanges("x",IRanges(start=c(1,28,50), end=c(12,40,60))),
#     GRanges("x",IRanges(start=c(1,30,48), end=c(10,42,60)))
# ))
#
# plot_feature_vpr(c(myref,mydata),featureCols=c(rep("firebrick",length(myref)),rep("steelblue",length(mydata))))
# tx_width(myref)
# tx_width(mydata)
# find_misaligned_small_exon(mydata,myref,smallExonSize=7,deltaDist=10)

### old implementation false positive, we don't check if the rest is exact match
find_misaligned_small_exon_v1 = function(query,subject, ignore.strand = FALSE, smallExonSize = 60, deltaDist=20){
    ## deltaDist is per exon difference, there could be multiple exon difference, although this is rare case
    ### subject should be reference, to compare against
    ## only allow one small exon difference
    if(!is_txClass(query) | !(is_txClass(subject))) {
        stop("Both query and subject have to be txClass")}

    query$index = seq_len(length(query))
    subject$index = seq_len(length(subject))

    if( all(nblocks(subject)<2) ) return(NULL)

    subject = subject[nblocks(subject)>2]
    internalExons = exonIV_to_intronIV(exonIV_to_intronIV(subject))

    hasSmallExon = width(blocks(internalExons)) <= smallExonSize
    subjectIndexWithSmallExon = which(any(hasSmallExon))

    if(length(subjectIndexWithSmallExon) ==0) return(NULL)

    subject = subject[subjectIndexWithSmallExon]
    hasSmallExon = hasSmallExon[subjectIndexWithSmallExon]

    queryIntrons = exonIV_to_intronIV(query)
    subjectIntrons = exonIV_to_intronIV(subject)

    ovlps = suppressWarnings(findOverlaps(query, subject, ignore.strand=ignore.strand))
    if(!length(ovlps)) return(NULL)
    isOnlyOneExonDiff = (nblocks(subject[subjectHits(ovlps)]) - nblocks(query[queryHits(ovlps)]))==1
    if(all(!isOnlyOneExonDiff)) return(NULL)
    ovlps = ovlps[isOnlyOneExonDiff]

    mydist = dist_between_gr(blocks(queryIntrons[queryHits(ovlps)]),blocks(subjectIntrons[subjectHits(ovlps)]),
        ignore.strand = ignore.strand )

    res = tibble(
        query_i = queryHits(ovlps), subject_i = subjectHits(ovlps),
        delta_size = mydist) %>%
            group_by(query_i) %>% summarise(
                subject_i = subject_i[which.min(delta_size)],
                n = n(),
                nhit = sum(delta_size == min(delta_size)),
                delta_size = min(delta_size),
            ) %>% filter(!is.infinite(delta_size)) %>% mutate(delta_size = as.integer(delta_size))
    if(!nrow(res)) return(NULL)

    res$n_exons_diff = nblocks(subject[res$subject_i]) - nblocks(query[res$query_i])
    ans = res %>% filter(n_exons_diff>0)
    if(!nrow(ans)) return(NULL)

    diffQueryToRef  = GenomicRanges::setdiff(blocks( subjectIntrons[ans$subject_i] ), blocks( queryIntrons[ans$query_i]) )
    diffRefToQuery  = GenomicRanges::setdiff(blocks( queryIntrons[ans$query_i] ), blocks( subjectIntrons[ans$subject_i]) )
    ans$len_missing_exon = sum(width(diffRefToQuery))
    ans$len_ss_added  = sum(width(diffQueryToRef))

    isDiffOnlySmallExon = all(width(diffRefToQuery) <= smallExonSize) & (lengths(diffRefToQuery) == ans$n_exons_diff)
    isDiffSmallerThanDeltaDist = abs(sum(width(diffQueryToRef)) - sum(width(diffRefToQuery))) <= (deltaDist * ans$n_exons_diff)

    if( !any(isDiffOnlySmallExon & isDiffSmallerThanDeltaDist)) return(NULL)
    ans = ans[isDiffOnlySmallExon & isDiffSmallerThanDeltaDist,] %>% filter(n_exons_diff==1)

    ## new way of verifying
    # validationDat = find_misaligned_small_exon_any(query[ans$query_i], subject[ans$subject_i],
    #     ignore.strand =ignore.strand,smallExonSize = smallExonSize, deltaDist = deltaDist, simplify=FALSE)
    # ans = ans[validationDat$query_i, ]
    # ans$delta_size = validationDat$delta_size

    ## NOTE: this piece of code is here inefficient and I am not sure if it works 100%
    ## we verify that the missiong small exon is indeed aligned to adjacent exons in the data
    # smallExonIndex = which(hasSmallExon)
    # validationDat = lapply( seq_len(nrow(ans)), function(i) {
    #     thisSmallExonIndex = smallExonIndex[[ ans$subject_i[[i]] ]]
    #     thisSubjectExons = unlist(blocks( subject[ans$subject_i[i]] ))
    #     thisQueryIntrons = unlist(blocks( queryIntrons[ans$query_i[i]] ))
    #     thisQueryExons = unlist(blocks(query[ans$query_i[i]]))
    #
    #     tmpOvlps = suppressWarnings(findOverlaps( thisSubjectExons[thisSmallExonIndex],thisQueryIntrons,
    #         ignore.strand=ignore.strand,type="within" ))
    #     if(length(tmpOvlps)==0) {
    #         return( tibble(is_hit=FALSE, delta_size=Inf) )
    #     } else {
    #         tmpvar = tibble(
    #             query_i = queryHits(tmpOvlps),
    #             subject_i = subjectHits(tmpOvlps),
    #             query_left = end(thisQueryExons)[subject_i],
    #             query_right = start(thisQueryExons)[subject_i+1],
    #             subject_left = end(thisSubjectExons)[thisSmallExonIndex[query_i]-1],
    #             subject_right = start(thisSubjectExons)[thisSmallExonIndex[query_i]+1]
    #         ) %>% mutate(d1=query_left - subject_left,d2=subject_right-query_right,
    #             d=width(thisSubjectExons[thisSmallExonIndex[queryHits(tmpOvlps)]]),
    #             delta_size=abs(d-d1-d2),
    #             is_hit = delta_size<=deltaDist & (d1 !=0 | d2 !=0) )
    #         return(tibble(is_hit=all(tmpvar$is_hit), delta_size=sum(tmpvar$delta_size)))
    #     }
    # })
    # validationDat = do.call(rbind, validationDat)
    # ans$delta_size = validationDat$delta_size
    # ans = ans[validationDat$is_hit,]
    ans = ans %>% mutate(query_i=query$index[query_i],subject_i=subject$index[subject_i])
    return(ans)
}

find_misaligned_small_exon = function(query,subject, ignore.strand = FALSE, smallExonSize = 60, deltaDist=20, simplify=TRUE){
    ## deltaDist is per exon difference, there could be multiple exon difference, although this is rare case
    ### subject should be reference, to compare against
    if(!is_txClass(query) | !(is_txClass(subject))) {
        stop("Both query and subject have to be txClass")}

    ## NOTE convert list to vector, generalise this idea
    myshift = function(x,lengthShift){
        ## x is a list, lengthShift is shift for for element of x
        myls = rep(c(0,cumsum(lengthShift[-length(lengthShift)])), lengths(x))
        sign(unlist(x)) *(abs(unlist(x)) + myls)
    }
    myunlist = function(x){
        ans = unlist(blocks(x))
        ans$index = rep(seq_len(length(x)),lengths(blocks(x)))
        ans
    }

    ans = find_misaligned_small_exon_any(query, subject,
            ignore.strand =ignore.strand,smallExonSize = smallExonSize, deltaDist = deltaDist, simplify=FALSE)

    if(nrow(ans)==0) return(NULL)
    ans$n_exons_diff = nblocks(subject[ans$subject_i]) - nblocks(query[ans$query_i])
    ans = ans %>% filter(n_exons_diff>0)
    if(nrow(ans)==0) return(NULL)

    ans$index = seq_len(nrow(ans))

    queryIntronIndex = as(
        split(-ans$ith_intron,ans$index) , "IntegerList")
    subjectNeighboringIntronIndex = as(split(
        -as.vector( t(cbind(ans$ith_exon - 1 , ans$ith_exon ))),
        rep(ans$index, each=2) ),"IntegerList")

    queryIntrons = exonIV_to_intronIV(query)
    subjectIntrons = exonIV_to_intronIV(subject)


    queryVar1 = myunlist(queryIntrons[ans$query_i])
    queryVar2 = queryVar1[myshift(queryIntronIndex,nblocks(queryIntrons[ans$query_i]))]

    subjectVar1 = myunlist(subjectIntrons[ans$subject_i])
    subjectVar2 = subjectVar1[myshift(subjectNeighboringIntronIndex,nblocks(subjectIntrons[ans$subject_i]))]

    queryVar = split( queryVar2,factor(queryVar2$index,seq_len(max(queryVar1$index)) ) )
    subjectVar = split( subjectVar2,factor(subjectVar2$index,seq_len(max(subjectVar1$index)) ) )
    wh = lengths(queryVar) == lengths(subjectVar)

    ### special case
    isHitSpecialCase = lengths(queryVar) == 0 & lengths(subjectVar) ==0
    ansSpecial = ans[isHitSpecialCase,]

    sel = wh & !isHitSpecialCase
    if(any(sel)){
        ansNormalCase = ans[sel, ] %>% mutate(
            queryClname = granges_to_ivname(asBED(queryVar[sel])),
            subjectClname =  granges_to_ivname(asBED(subjectVar[sel]))) %>%
            filter(queryClname==subjectClname) %>% select(-queryClname, -subjectClname)

        res = rbind(ansNormalCase, ansSpecial)
    } else{
        res = ansSpecial
    }
    res = res %>% arrange(query_i) %>% select(-index)

    if(nrow(res)==0) return(NULL)

    if(simplify){
        res = res %>%
            group_by(query_i) %>% summarise(
                subject_i = subject_i[which.min(delta_size)],
                ith_intron = ith_intron[which.min(delta_size)],
                ith_exon = ith_exon[which.min(delta_size)],
                exon_size = exon_size[which.min(delta_size)],
                ss_len = ss_len[which.min(delta_size)],
                delta_size = min(delta_size),
                n_exons_diff = n_exons_diff[which.min(delta_size)]
            )
    }

    return(res)
}

find_misaligned_small_exon_any = function(query,subject, ignore.strand = FALSE, smallExonSize = 60,
    deltaDist=20, simplify=TRUE){
    ## deltaDist is per exon difference, there could be multiple exon difference, although this is rare case
    ### subject should be reference, to compare against
    ## only allow one small exon difference
    if(!is_txClass(query) | !(is_txClass(subject))) {
        stop("Both query and subject have to be txClass")}
    defaultAns = tibble(query_i=integer(), subject_i = integer(),ith_intron = integer(),
        ith_exon=integer(), exon_size=integer(), ss_len=integer(), delta_size= integer())

    myfunc = function(x){
        ans = unlist(blocks(x))
        ### ith tx
        ans$tx_i = rep(seq_len(length(x)), lengths(blocks(x)))
        ### ith exon/intron of tx
        ans$pos_i = unlist(lapply(lengths(blocks(x)), seq_len))
        ans$id = paste(ans$tx_i, ans$pos_i,sep="_")
        names(ans) = ans$id
        return(ans)
    }

    ## for reducing the data later
    query$index = seq_len(length(query))
    subject$index = seq_len(length(subject))

    if( all(nblocks(subject)<3) ) return(defaultAns)

    subject = subject[nblocks(subject)>2]
    internalExons = exonIV_to_intronIV(exonIV_to_intronIV(subject))

    hasSmallExon = width(blocks(internalExons)) <= smallExonSize
    subjectIndexWithSmallExon = which(any(hasSmallExon))

    if(length(subjectIndexWithSmallExon) ==0) return(defaultAns)

    subject = subject[subjectIndexWithSmallExon]
    hasSmallExon = hasSmallExon[subjectIndexWithSmallExon]
    internalExons = internalExons[subjectIndexWithSmallExon]

    queryIntrons = exonIV_to_intronIV(query)
    subjectIntrons = exonIV_to_intronIV(subject)

    subjectIntronsList = myfunc(subjectIntrons)
    queryIntronsList = myfunc(queryIntrons)
    subjectInternalExonsList = myfunc(internalExons)
    smallExonList = subjectInternalExonsList[width(subjectInternalExonsList)<=smallExonSize]

    ### NOTE: is it necessary type within
    ovlps = suppressWarnings(findOverlaps(smallExonList, queryIntronsList,
        ignore.strand=ignore.strand, type="within"))

    if(!length(ovlps)) return(defaultAns)

    id_intron_prec = paste(smallExonList$tx_i[queryHits(ovlps)], smallExonList$pos_i[queryHits(ovlps)], sep = "_")
    id_intron_succ = paste(smallExonList$tx_i[queryHits(ovlps)], smallExonList$pos_i[queryHits(ovlps)] +1, sep="_")

    ans = tibble(
        query_i = queryIntronsList$tx_i[subjectHits(ovlps)], subject_i = smallExonList$tx_i[queryHits(ovlps)] ,
        ith_intron = queryIntronsList$pos_i[subjectHits(ovlps)] ,
        ith_exon = smallExonList$pos_i[queryHits(ovlps)] + 1 ,
        exon_size = width(smallExonList[paste(smallExonList$tx_i[queryHits(ovlps)] , ith_exon-1, sep="_")]),
        ### both should be positive
        dist_left = start(queryIntronsList[subjectHits(ovlps)])  - start(subjectIntronsList[id_intron_prec]),
        dist_right = end(subjectIntronsList[id_intron_succ]) - end(queryIntronsList[subjectHits(ovlps)]),
        ss_len = dist_left + dist_right ,
        delta_size = abs(ss_len - exon_size)
    ) %>% filter(delta_size<=deltaDist, ss_len>0)

    if(nrow(ans)==0) return(defaultAns)

    if(simplify)
        ans = ans %>%
        ###NOTE is simplify a good idea here? we need to do it to avoid multiple mapping
            group_by(query_i) %>% summarise(
                subject_i = subject_i[which.min(delta_size)],
                ith_intron = ith_intron[which.min(delta_size)],
                ith_exon = ith_exon[which.min(delta_size)],
                exon_size = exon_size[which.min(delta_size)],
                ss_len = ss_len[which.min(delta_size)],
                delta_size = min(delta_size)
            )
    ## convert back the index
    ans$query_i = query$index[ans$query_i]
    ans$subject_i = subject$index[ans$subject_i]

    return(ans)
}

## for ie_var and ei_var, we would like to know how long extension is
find_dist_at_one_end = function(query, subject, diffEndType="start"){
    ## endType start for distance at start, end for distance at end
    ## this indicates the other end most be the same
    if( !diffEndType %in% c("start","end") ) {
        stop("endType must be either start or end")
    }

    ovlps = suppressWarnings(findOverlaps(query,subject,type = ifelse(diffEndType=="start","end","start")))
    if(length(ovlps)==0) return(NULL)

    distToNearest = if(diffEndType=="start"){
        start( subject[ subjectHits(ovlps) ] ) - start(query[queryHits(ovlps)])
    }  else {
        end(query[queryHits(ovlps)]) - end( subject[ subjectHits(ovlps) ] )
    }

    mymakename = function(chr,strd,pos1,pos2){
        if(!all(pos2>=pos1)) stop("pos2 have to be greater or equal to pos1")
        paste0(chr,"_",pos1,"_",strd,"-",chr,"_",pos2,"_",strd)
    }

    ans = tibble(
        query_i = queryHits(ovlps),
        subject_i = subjectHits(ovlps),
        dist_to_nearest = distToNearest,
        name = if(diffEndType=="start"){
            mystart = ifelse(dist_to_nearest>0, start(query)[query_i] - dist_to_nearest,start(query)[query_i])
            myend = ifelse(dist_to_nearest>0, start(query)[query_i] ,start(query)[query_i] - dist_to_nearest)
            mymakename(
                as.character(seqnames(query)[query_i]),
                represent_minus_strand(as.character(strand(query)[query_i])),
                mystart,myend
            )
        } else {
            mystart = ifelse(dist_to_nearest>0, start(query)[query_i] ,start(query)[query_i] + dist_to_nearest)
            myend = ifelse(dist_to_nearest>0, start(query)[query_i] + dist_to_nearest,start(query)[query_i])
            mymakename(
                as.character(seqnames(query)[query_i]),
                represent_minus_strand(as.character(strand(query)[query_i])),
                mystart,myend
            )
        }
    )

    ans %>% group_by(query_i) %>% summarise(
        subject_i=subject_i[which.min(abs(dist_to_nearest))],
        dist_to_nearest = dist_to_nearest[which.min(abs(dist_to_nearest))],
        name=name[which.min(abs(dist_to_nearest))] )

    return(ans)
}


### code to plot all cases report find_misaligned_small_exon_any
# mysanitise = function(x){
#     mcols(x)[, setdiff(colnames(mcols(x)),"blocks")]  = NULL
#     return(x)
# }
# pdf("alignment_small_exon.pdf")
# for(i in seq_len(nrow(ans)) ) {
#     grid.newpage()
#     myref = mysanitise(subject[ans$subject_i[i]])
#     mydata = mysanitise(query[ans$query_i[i]])
#     plot_feature_vpr(suppressWarnings(c(myref, mydata )),featureCols=c(rep("firebrick",length(myref)),rep("steelblue",length(mydata))))
# }
# dev.off()
########################################################################
## FIXME both find_misaligned_small_exon and find_misaligned_ss_shifted we need to control
## if the difference occur at the two ends of the same introns, this is currently not the case
## resulting in possible false matching
# smallExonAtEndsCutoff = 30
# txWithSmallExonAtEnds = txMultiExonic[width(extract_nthExon(txMultiExonic,1))<smallExonAtEndsCutoff | width(extract_nthExon(txMultiExonic,nblocks(txMultiExonic)))<smallExonAtEndsCutoff]
# length(txWithSmallExonAtEnds)/
#
# x= asBED( GRangesList(GRanges("x",IRanges(start=c(25,40), end=c(30,50)))) )
# y= asBED( GRangesList(GRanges("x",IRanges(start=c(1,20,40), end=c(10,30,50))) ))
#
#
# x = GRanges("x",IRanges(start=c(1,20,40), end=c(10,30,50)))
# x1= asBED( GRangesList(GRanges("x",IRanges(start=c(1,25,42), end=c(15,30,50))), x) )
# y1= asBED( GRangesList( x,x ))


# x=asBED(GRangesList(GRanges("x",IRanges(c(1,20,30),c(10,25,40)),strand="+")))
# y=asBED(GRangesList(GRanges("x",IRanges(c(5,28),c(12,40)),strand="+")))
#
# find_misaligned_small_exon(y,x)
#find_misaligned_small_exon(thisTxCluster$consensus,tx)
#find_misaligned_ss_shifted(thisTxCluster$consensus, tx, ignore.strand=FALSE,maxDist=15)


### for 5 TSO event
## slow implementation but robust
advanced_filter2_v1 = function(x,cutoff=0.2,ncpu=1){
    ## x is a cluster filter object with slots consensus, filter, tu
    ## this adds an additional splot pass_falsetso or overrites the existing
    ## in the filter matrix
    ## this filter remove transcript that is likely to arise due to early TSO events
    ## in the same transricpt unit
    ## if the distance is shorter by given cutoff minDist. At this stage, we do not have a relative count
    ## threshold.
    ## returns cluster filer object with addtional column pass_falsetso in the filter matrix
    ## currently I only compare thisTx to all tx that are longer in 5' ends
    ## closest tx in distance, that are different at 5 ends only

    idxKeep = which2keep(x$filter)

    x$filter$pass_falsetso = as.logical(NA)
    x$filter$pass_falsetso[idxKeep] = FALSE

    mymat = x$stat
    mymat$index = 1:nrow(mymat)
    mymat = mymat[idxKeep,]

    ovlps = findOverlaps(x$consensus[idxKeep], x$tu, type="within")
    message(length(unique(subjectHits(ovlps)))," total clusters to filter")

    stopifnot(nrow(mymat) == length(subjectHits(ovlps)))

    dat = split(mymat, subjectHits(ovlps))
    ans = mclapply((1:length(dat)), function(ith){
            thisTu=dat[[ith]]
            message( "process ", ith," th Cluster" )

            y = thisTu[order(thisTu$count,decreasing=TRUE),]
            if(nrow(y)>1){
                shouldKeep = rep(TRUE, nrow(y))
                for(i in 2:nrow(y)){
                    #message(i)
                    thistx = x$consensus[y$index[i]]
                    isDiffAtFiveSite = sapply((1:(i-1))[shouldKeep[1:(i-1)]], function(j){
                        reftx = x$consensus[y$index[j]]
                        suppressWarnings(is_different_at_five_end(thistx,reftx))
                    })
                    if(any(isDiffAtFiveSite)){
                        ### picking the clostest one in terms of distance
                        # diffFiveSite = sapply((1:(i-1))[shouldKeep[1:(i-1)]][isDiffAtFiveSite], function(j){
                        #     reftx = x$consensus[y$index[j]]
                        #     dist_between_tx_intron_size(thistx, reftx)
                        # })
                        #reftxIdx = y$index[1:(i-1)][(shouldKeep[1:(i-1)])][isDiffAtFiveSite][which.min(diffFiveSite)]
                        reftxIdx = y$index[1:(i-1)][(shouldKeep[1:(i-1)])][isDiffAtFiveSite]
                        #if(y$count[i] < x$stat$count[reftxIdx] *cutoff){
                        #if(y$count[i] < max(x$stat$count[reftxIdx]) *cutoff){
                        if(y$count[i] < sum(x$stat$count[reftxIdx]) *cutoff){
                            shouldKeep[i] = FALSE
                        }
                    }
                }
                return(y$index[shouldKeep])
            } else { return(y$index) }
            },
         mc.cores=ncpu)
    ans = unlist(ans)
    x$filter$pass_falsetso[ans]= TRUE
    return(x)
}

## fast implement but there are some edge casess need to be looked at
advanced_filter2_v2 = function(x,cutoff=0.2,ncpu=1){
    ## x is a cluster filter object with slots consensus, filter, tu
    ## this adds an additional splot pass_falsetso or overrites the existing
    ## in the filter matrix
    ## this filter remove transcript that is likely to arise due to early TSO events
    ## in the same transricpt unit
    ## if the distance is shorter by given cutoff minDist. At this stage, we do not have a relative count
    ## threshold.
    ## returns cluster filer object with addtional column pass_falsetso in the filter matrix
    ## currently I only compare thisTx to all tx that are longer in 5' ends
    ## closest tx in distance, that are different at 5 ends only
    x$filter$pass_falsetso = NULL ## make sure we initialise
    idxKeep = which2keep(x$filter)

    x$filter$pass_falsetso = as.logical(NA)
    #x$filter$pass_falsetso[idxKeep] = FALSE
    x$filter$pass_falsetso[idxKeep] = TRUE

    mymat = x$stat
    mycons = x$consensus
    mymat$index = mycons$index = 1:nrow(mymat)

    mymat = mymat[idxKeep,]
    mycons = mycons[idxKeep]

    ovlps = findOverlaps(mycons, x$tu, type="within")
    message(length(unique(subjectHits(ovlps)))," total clusters to filter")

    stopifnot(nrow(mymat) == length(subjectHits(ovlps)))

    diffAtEitherEnd = is_different_at_either_end(mycons, mycons, ncpu=ncpu) %>%
        mutate(strand=as.character(strand(mycons))[query_i]) %>%
        filter(query_i!=subject_i, strand!="*") %>%
        mutate(is_five_var= (strand =="+" & type=="left") | (strand =="-" & type=="right")) %>%
        filter(is_five_var) %>% mutate(query_i = mycons$index[query_i], subject_i = mycons$index[subject_i])

    dat = split(mymat, subjectHits(ovlps))

    ## 5' truncation product less than cutoff * all transcripts higher/equal this gets filtered out
    diffAtEitherEndWithCount = diffAtEitherEnd %>% mutate(
        count_query=mymat$count[match(query_i,mymat$index)],
        count_subject=mymat$count[match(subject_i,mymat$index)]) %>% group_by(query_i) %>% summarise(
        count_query=unique(count_query), count_query_sum= sum(count_subject[count_subject>=count_query]))
    tofilter = diffAtEitherEndWithCount %>% filter(count_query < count_query_sum * cutoff)

    # ans = mclapply((1:length(dat)), function(ith){
    #         message( "process ", ith," th Cluster" )
    #
    #         y = dat[[ith]] %>% arrange(desc(count))
    #         if(nrow(y)>1){
    #             keep = rep(TRUE, nrow(y))
    #             for(i in 2:nrow(y)){
    #                 #message(i)
    #                 thistx = y$index[i]
    #                 reftx = y$index[(1:(i-1))[keep[1:(i-1)]]]
    #
    #                 isDiffAtFiveSite = sapply(reftx, function(j){
    #                     nrow(diffAtEitherEnd %>% filter(query_i == thistx, subject_i==j))>0
    #                 })
    #
    #                 if(any(isDiffAtFiveSite)){
    #                     ### picking the clostest one in terms of distance
    #                     # diffFiveSite = sapply((1:(i-1))[shouldKeep[1:(i-1)]][isDiffAtFiveSite], function(j){
    #                     #     reftx = x$consensus[y$index[j]]
    #                     #     dist_between_tx_intron_size(thistx, reftx)
    #                     # })
    #                     #reftxIdx = y$index[1:(i-1)][(shouldKeep[1:(i-1)])][isDiffAtFiveSite][which.min(diffFiveSite)]
    #                     reftxIdx = y$index[1:(i-1)][(keep[1:(i-1)])][isDiffAtFiveSite]
    #                     #if(y$count[i] < x$stat$count[reftxIdx] *cutoff){
    #                     #if(y$count[i] < max(x$stat$count[reftxIdx]) *cutoff){
    #                     if(y$count[i] < sum(x$stat$count[reftxIdx]) *cutoff){
    #                         keep[i] = FALSE
    #                     }
    #                 }
    #             }
    #             return(y$index[keep])
    #         } else { return(y$index) }
    #         },
    #      mc.cores=ncpu)
    # ans = unlist(ans)
    # x$filter$pass_falsetso[ans]= TRUE
    x$filter$pass_falsetso[tofilter$query_i]= FALSE
    return(x)
}

advanced_filter2 = advanced_filter2_v2

find_tu = function(x){
    ## find transcript unit for tx clusters
    x$tu = GenomicRanges::reduce(
        x$consensus[which2keep(x$filter)], ignore.strand=FALSE)
    return(x)
}

threeend = function(x){
    ## x is GRanges
    GRanges(
        seqnames(x),
        IRanges(ifelse(strand(x)=="+", end(x), start(x)),width=1), strand(x)
    )
}

has_flag = function(x,flag){
    as.logical(bitwAnd(x,flag))
}

is_control= function(x){
    ## x is txClass
    chr = as.character(seqnames(x))
    grepl("ERCC",chr) | chr == "chrIS" | chr == "SIRV"
}

extract_nth_list = function(x, i){
    ## x is list , i is index in each element
    ## return i th element of each entry in list(-like) x
    ne = lengths(x)
    if(any(i > ne)) stop("index out of range for the list x")
    indexShift = cumsum(ne)
    unlist(x)[c(0, indexShift[-length(indexShift)]) + i]
}

extract_nthExon = function(x,nthExon){
    ## x is tx class, i.e GRanges with blocks
    ## nthExon is n th exon irrelevant of strand
    extract_nth_list(blocks(x),nthExon)
}

threeutr = function(x){
    ## x is tx class, i.e GRanges with blocks
    if(!all(strand(x) %in% c("+","-"))) stop("All strand must be either + or -.")
    nthExon = ifelse(strand(x)=="+",lengths(blocks(x)), 1)
    extract_nthExon(x,nthExon)
}

distanceToNearestTx = function(query, subject, ignore.strand=FALSE, simplify=TRUE ) {
    ## both query and subjct are GRanges, with column blocks for exons
    ## query can be my read, subject is the reference
    ## the delta size is based on intron space!!!
    stopifnot(is(query,"GRanges") & is(subject,"GRanges"))

    #ovlps = suppressMessages(findOverlaps(exonIV_to_intron_IV(query), exonIV_to_intron_IV(subject), ignore.strand=ignore.strand))
    ovlps = suppressMessages(findOverlaps(query, subject, ignore.strand=ignore.strand))
    ### nothing overlaps
    if(length(ovlps)==0){ return(NULL) }

    mydist = dist_between_tx_intron_size( query[queryHits(ovlps)], subject[subjectHits(ovlps)] , ignore.strand=ignore.strand)

    ## mutual overlap
    res = tibble(
        query_i = queryHits(ovlps),
        subject_i = subjectHits(ovlps),
        delta_size = mydist
        )
    if(simplify){
        ans = res %>% group_by(query_i) %>% summarise(
            subject_i = subject_i[which.min(delta_size)],
            n = n(),
            nhit = sum(delta_size == min(delta_size)),
            delta_size = min(delta_size),
        )
    } else {
        ans = res
    }
    #ans[is.finite(ans$delta_size),]
    ans = ans %>% filter(!is.infinite(delta_size)) %>% mutate(delta_size = as.integer(delta_size) )
    return(ans)
}

## working with txFilterClass
subset_txFilterClass = function(x, index){
    list(
        consensus=x$consensus[index],
        stat= x$stat[index,],
        "filter"= x$filter[index,],
        "clname" = x$clname[index],
        "associated_tx" = x$associated_tx[index],
        "associated_refannot" = x$associated_refannot[index,]
    )
}

validate_txFilterClass = function(x) {
    length(x$consensus) == nrow(x$filter)
}

merge_txFilterClass = function(x, mapping){
    ### given a remapping matrix, we merge entries in x
    ### this is mainly to deal with alignment artefacts
    if( all(mapping$keep) ) return(x)
    if(length(unique(mapping$index)) != length(x$consensus) ) stop("Mapping should be as long as x")
    myfac = factor( mapping$index_reassigned, sort(unique(mapping$index_reassigned)) )
    indexKeep = as.integer(levels(myfac))
    ans = list(
        consensus = x$consensus[indexKeep],
        stat  = tibble(
            clname = x$stat$clname[indexKeep],
            count = tapply(x$stat$count,myfac, sum),
            n_full_length  = tapply(x$stat$n_full_length,myfac, sum),
            perc_full_length = n_full_length/count,
            ### FIXME do we really need the following columns?
            chr = x$stat$chr[indexKeep],
            start = x$stat$start[indexKeep],
            end = x$stat$end[indexKeep],
            strand = x$stat$strand[indexKeep],
            d_polyARich = x$stat$d_polyARich[indexKeep]
        ),
        filter = x$filter[indexKeep,],
        clname=x$clname[indexKeep],
        associated_tx = x$associated_tx[indexKeep],
        associated_refannot = x$associated_refannot[indexKeep,],
        mapping= mapping,
        is_remapped = TRUE,
        filter_org = x$filter,
        consensus_org = x$consensus,
        stat_org = x$stat
    )
    return(ans)
}

simplify_diff_at_either_end = function(x){
    ## x is tibble from is_different_at_either_end
    ## return value:
    x %>% group_by(query_i) %>% summarise(n=n(), nhit= sum(nexon_diff==min(nexon_diff)),
    ## if multiple hits with the same minimum exon num diff, we choose the one with smallest size difference
    wh = if(nhit==1){which.min(nexon_diff)} else {sel=which(nexon_diff==min(nexon_diff)); sel[which.min(size_diff[sel])] },
    subject_i = subject_i[wh],nexon_diff=nexon_diff[wh],type=type[wh],size_diff=size_diff[wh] )
}

quick_comp = function(x,y, ncpu=20){
   rv = mclapply(split( seq_len(length(x)), cut(seq_len(length(x)),ncpu)) , function(i){
        x[i] %in% y[i]
    }, mc.cores=ncpu)
    Reduce(c, rv)
}

is_different_at_either_end = function(query, subject, ignore.strand=FALSE, simplify=FALSE, ncpu=20){
    ## return a tible with index in query subject, number of different exons, type of difference,
    ## total size of exons that are difffer
    if(!is_txClass(query) | !is_txClass(subject))
        stop("Both query and subject should be txClass")
    if(!all(is_multiexonic(query)) | !all(is_multiexonic(subject)))
        stop("Both query and subject should be multiexonic")

    queryIntrons = exonIV_to_intronIV(query)
    subjectIntrons = exonIV_to_intronIV(subject)

    ovlps = suppressWarnings(findOverlaps(blocks(queryIntrons), blocks(subjectIntrons),type="within", ignore.strand=ignore.strand))
    if(length(ovlps)==0) return(NULL)

    isShorter = nblocks(queryIntrons[queryHits(ovlps)]) < nblocks(subjectIntrons[subjectHits(ovlps)])
    if(!any(isShorter)) return(NULL)
    ovlps = ovlps[isShorter]

    isDifferAtEnds = all(suppressWarnings(
        quick_comp(
            blocks(queryIntrons[queryHits(ovlps)]),
            blocks(subjectIntrons[subjectHits(ovlps)]), ncpu=ncpu) ) )
    if(!any(isDifferAtEnds)) return(NULL)
    ovlps = ovlps[isDifferAtEnds]

    ### infer which end or both ends and the
    subjectIndexSharedIntrons =  which(suppressWarnings(
        quick_comp(blocks(subjectIntrons[subjectHits(ovlps)]), blocks(queryIntrons[queryHits(ovlps)]), ncpu=ncpu))
    )
    ### there are cases the middle introns are skipped, this is not considered 5' 3 variation
    ## FIXME a new function taking advantage of this to detect intron retention?
    isIntronRetained =  any(diff(subjectIndexSharedIntrons)>1)

    subjectAllIntronIndex = as(mclapply(nblocks(subjectIntrons[subjectHits(ovlps)]), seq_len, mc.cores=ncpu),"CompressedIntegerList")
    spIndex = split(seq_len(length(subjectAllIntronIndex)), cut(seq_len(length(subjectAllIntronIndex)), ncpu))
    diffIndex = mclapply(spIndex, function(i) {
        setdiff(subjectAllIntronIndex[i], subjectIndexSharedIntrons[i])
    } ,mc.cores=ncpu)
    diffIndex = Reduce(c, diffIndex)
    ## downstrem subsetting does not work with names list for some unknown reason
    names(diffIndex) = NULL

    isLeftEndDiffer = any(diffIndex < min(subjectIndexSharedIntrons))
    isRightEndDiffer = any(diffIndex > max(subjectIndexSharedIntrons))
    exonHitWidth = width(blocks(subject[subjectHits(ovlps)]))

    ans = tibble(
        query_i = queryHits(ovlps),
        subject_i = subjectHits(ovlps),
        nexon_diff = lengths(diffIndex),
        type = ifelse( isLeftEndDiffer & isRightEndDiffer, "both", ifelse(isLeftEndDiffer, "left", ifelse(isRightEndDiffer,"right","none"))) ,  ## left, right, both
        size_diff = sum(exonHitWidth[diffIndex[diffIndex<min(subjectIndexSharedIntrons)]]) +
            ## for right end, plus 1 to convert from intron index to exon index
            sum(exonHitWidth[diffIndex[diffIndex>max(subjectIndexSharedIntrons)]+1] ))

    ans = ans[!isIntronRetained,]
    if(simplify) ans = simplify_diff_at_either_end(ans)
    return(ans)
}


# export_to_bed = function(x, outfile="test.bed.gz"){
#     ## x is a tx clustering object
#     ## with filter class
#     require(rtracklayer)
#     filterBitcodes =calc_bitcode_filter(!x$filter)
#     mycol = ifelse(
#         ## color for pass all filters
#         !filterBitcodes,"red",
#         ### highlight the ones that failed polyA rich and TSO artefact filter
#         ifelse(!x$filter$pass_polyArich,"blue",
#             ifelse(!x$filter$pass_falsetso,"green","black"))
#     cs = x$consensus
#     cs$name = paste(
#         paste0("F",filterBitcodes),
#         paste0("C",x$stat$count), names(cs),
#         sep=":")
#     cs$itemRgb = mycol
#     export(cs[x$keep],outfile)
# }


# x1 = to_grangs_w_blocks(clname_to_granges(xx))
# x2 = to_grangs_w_blocks(clname_to_granges(xx))
#
#
# x=GRanges("x",IRanges(c(1,20,30),c(10,25,35)))
# y=GRanges("x",IRanges(c(20,30),c(25,35)))
# library(rtracklayer)
# x=asBED(GRangesList(GRanges("x",IRanges(c(1,10,30,40),c(5,20,35,45)),strand="+")))
# y=asBED(GRangesList(GRanges("x",IRanges(c(10,30),c(20,35)),strand="+")))
# x2=asBED(GRangesList(GRanges("x",IRanges(c(3,10,30),c(5,20,35)),strand="+")))
# is_different_at_five_end(y,x)
