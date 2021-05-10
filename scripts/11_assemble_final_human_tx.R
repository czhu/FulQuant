library(tidyverse)
library(rtracklayer)
projectFolder = "."
SCRIPTDIR = file.path(projectFolder, "sw")
GENOMEDIR = file.path(projectFolder, "genome")
source(file.path(SCRIPTDIR, "clustering_functions.R"))

## change here
labPrefix = "My"
labSuffix = "Lab"

infolder = file.path(projectFolder, "combined/tx_annot")

load( file.path(infolder, "novel_cluster_class.rda") )
load( file.path(infolder, "txCluster.rda") )
load( file.path(GENOMEDIR, "tx.rda"))

myclusters = subset(myclusters, is_human_final)
myclusters$type = clusterClass$type[match(myclusters$clname, clusterClass$clname)]
myclusters$type[is.na(myclusters$type)] = "UNKNOWN"
myclusters$type[myclusters$match_tx_pass | myclusters$match_tx_fail] = "REF"

novelClNameWithEndVariation = clusterClass$clname[clusterClass$is_end_variation]
novelClNameWithAlignIssues = clusterClass$clname[ clusterClass$align_issue!="none" ]

myclusters$is_end_variation = FALSE
myclusters$is_end_variation[myclusters$clname %in% novelClNameWithEndVariation] = TRUE

myclusters$has_align_issue = FALSE
myclusters$has_align_issue[myclusters$clname %in% novelClNameWithAlignIssues] = TRUE

COLOR_SCHEME_HUMAN_FINAL_TX = c(
    "base" = "black",
    "match_ref_filter_pass" = "#4daf4a", ## match_tx_pass
    "match_ref_filter_fail" = "#73D171", ## match_tx_fail light green
    "novel_align_issue" =  "#a65628",
    "novel_ir" = "#e41a1c" ,
    "novel_es" = "#e41a1c" ,
    "novel_nc" = "#e41a1c" ,
    "novel_ne" = "#e41a1c" ,
    "novel_ev" = "#e41a1c",
    "novel_multi" = "#e41a1c"
)

mt = match(myclusters$clname,clusterClass$clname)
myclusters$align_issue = clusterClass$align_issue[mt]

## Benjamin's data has higher count
finalTx = subset(myclusters, !is_end_variation & (!has_align_issue | (novel_pass & count >= 10)))

## new color code
finalTx$itemRgb = COLOR_SCHEME_HUMAN_FINAL_TX["base"]
finalTx$itemRgb[finalTx$match_tx_fail] = COLOR_SCHEME_HUMAN_FINAL_TX["match_ref_filter_fail"]
finalTx$itemRgb[finalTx$clname %in%  clusterClass$clname[ clusterClass$type == "IR" ] ]  = COLOR_SCHEME_HUMAN_FINAL_TX["novel_ir"]
finalTx$itemRgb[finalTx$clname %in%  clusterClass$clname[ clusterClass$type == "ES" ] ]  = COLOR_SCHEME_HUMAN_FINAL_TX["novel_es"]
finalTx$itemRgb[finalTx$clname %in%  clusterClass$clname[ clusterClass$type == "NC" ] ]  = COLOR_SCHEME_HUMAN_FINAL_TX["novel_nc"]
finalTx$itemRgb[finalTx$clname %in%  clusterClass$clname[ clusterClass$type == "EV" ] ]  = COLOR_SCHEME_HUMAN_FINAL_TX["novel_ev"]
finalTx$itemRgb[finalTx$clname %in%  clusterClass$clname[ clusterClass$type %in% c("NE","ASS") ] ]  = COLOR_SCHEME_HUMAN_FINAL_TX["novel_ne"]
finalTx$itemRgb[finalTx$clname %in%  clusterClass$clname[ clusterClass$type == "MN" ] ]  = COLOR_SCHEME_HUMAN_FINAL_TX["novel_multi"]

stopifnot(all(finalTx$itemRgb %in% COLOR_SCHEME_HUMAN_FINAL_TX))

## NOTE remove novel tx with only 1 intron
hasNoGeneAssigned = is.na( finalTx$associated_gene )
message( sum( hasNoGeneAssigned ) ," ncRNA to assign gene name")

ncRNATx = finalTx[ hasNoGeneAssigned ]
ncRNAloci = GenomicRanges::reduce( ncRNATx, ignore.strand=FALSE)
## XXX strand must be defined, we have a case with undefined strand
ncRNAloci = subset(ncRNAloci, strand(ncRNAloci) != "*")
tx2locus = findOverlaps(ncRNATx, ncRNAloci, ignore.strand=FALSE, type = "within")
stopifnot( all(lengths(as.list(tx2locus))==1) )

## assign gene tx name etc
finalTx$gene_name = mygenes$gene_name[match(finalTx$associated_gene, mygenes$gene_id)]
## gene_source to source GENCODE
finalTx$gene_source = as.character(mygenes$source[match(finalTx$associated_gene, mygenes$gene_id)])
finalTx$gene_type = mygenes$gene_type[match(finalTx$associated_gene, mygenes$gene_id)]
finalTx$transcript_name = tx$transcript_name[match(finalTx$associated_tx, tx$transcript_id)]
## transcript_source to source GENCODE
finalTx$transcript_source = as.character(tx$source[match(finalTx$associated_tx, tx$transcript_id)])

allowedGeneBiotype = c("protein_coding","lincRNA")
ovlps = findOverlaps( ncRNAloci, mygenes[mygenes$gene_type %in% allowedGeneBiotype] )

stopifnot( length(ovlps) ==0 )
## define antise
ovlps2 = findOverlaps(invertStrand(ncRNAloci), mygenes[mygenes$gene_type %in% allowedGeneBiotype])
locusType = ifelse(lengths( as.list(ovlps2) ) > 0, "antisense" , "intergenic" )
locusName = sprintf("NCRNA%05d", 1:length(ncRNAloci))

ncRNAType = as_tibble(tx2locus) %>% mutate(
    clname = ncRNATx$clname[queryHits],
    locus_name = locusName[subjectHits],
    locus_type = locusType[subjectHits]
)

ncRNAType = ncRNAType %>% mutate(count = ncRNATx$count[queryHits]) %>%
    group_by(subjectHits) %>% mutate(rank=rank( -count, tie="first" ),
    tx_name = paste(locus_name, rank, sep="-"),
    locus_range = granges_to_igvCoord(ncRNAloci)[subjectHits])

stopifnot( identical(finalTx$clname[hasNoGeneAssigned], ncRNAType$clname) )

finalTx$gene_name[hasNoGeneAssigned] = ncRNAType$locus_name
finalTx$gene_source[hasNoGeneAssigned] = labPrefix
finalTx$gene_type[hasNoGeneAssigned] = paste(labPrefix,ncRNAType$locus_type, sep="_")
finalTx$transcript_name[hasNoGeneAssigned] = ncRNAType$tx_name
finalTx$transcript_source[hasNoGeneAssigned] = labPrefix

## assign tx name
hasNoTxName = is.na(finalTx$transcript_name)
stopifnot( !anyNA(finalTx$associated_gene[hasNoTxName]))

txCountOrderInKnownGenes = as_tibble(finalTx[hasNoTxName,c("count","associated_gene", "clname")]) %>%
    group_by(associated_gene) %>% mutate(rank=rank(-count,tie="first"))

stopifnot( identical(txCountOrderInKnownGenes$clname, finalTx$clname[hasNoTxName]) )

finalTx$transcript_name[hasNoTxName] = paste0( finalTx$gene_name[hasNoTxName],
    "-",labSuffix,txCountOrderInKnownGenes$rank )
finalTx$transcript_source[hasNoTxName] = labPrefix

stopifnot( !anyNA(finalTx$transcript_name) )
stopifnot( !anyNA(finalTx$transcript_source) )
stopifnot( !anyNA(finalTx$gene_name) )
stopifnot( !anyNA(finalTx$gene_source) )
stopifnot( !anyNA(finalTx$gene_type))

## define read through
mygenesAllowed = mygenes[mygenes$gene_type %in% allowedGeneBiotype]
mygenesReduced = GenomicRanges::reduce(mygenesAllowed)
ovlps3 = findOverlaps(finalTx, mygenesReduced)
hasReadThrough = lengths( as.list(ovlps3) ) > 1

readThrough = ifelse(hasReadThrough,
    sapply( as.list(findOverlaps(finalTx[which(hasReadThrough)], mygenesAllowed)),
        function(x) paste(mygenesAllowed$gene_name[x], collapse=";")),
    as.character(NA))

finalTx$read_through = readThrough


finalTx$name = finalTx$transcript_name

export(myclusters, file.path(infolder, "tx_human_full.bed.gz") )
export(finalTx, file.path(infolder, "tx_human_final.bed.gz") )
txAnnot = mcols(finalTx)[, c( "name", "is_human_final", "gene_name", "transcript_name", "count", "type", "displayName","score", "igv_coord")]
write.table(txAnnot, file.path(infolder, "tx_human_final_annot.txt"), quote=FALSE, sep="\t",row.names=FALSE)

save(myclusters,finalTx, file=file.path(infolder, "tx_human.rda") )
