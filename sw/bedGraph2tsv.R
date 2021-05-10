## convert the bed file coverage for clustering with Wave's java software
## cannot work with bedGraph!!! because bedGraph's range is not supported in such conversion
library(tidyverse)
#library(rtracklayer)

args = commandArgs(trailingOnly = TRUE)
infile = args[1] ## reads file
outfile = paste0(tools::file_path_sans_ext(infile),".reformatted.tsv")

hasHeader = grepl("^chrom",readLines(infile,1))

if(hasHeader){
    x = read_tsv(infile,col_names=hasHeader, col_types=cols(chrom="c"))
} else {
    x = read_tsv(infile,col_names=hasHeader, col_types=cols(X1="c"))
}

tb = tibble(chr=x[[1]],strand="+",
    ## plus1 bcause it's bed format
    pos=x[[2]]+1L,count=x[[3]]) %>% arrange(chr,desc(count))
tb$count = as.integer(tb$count)

write_tsv(tb,outfile, col_names=FALSE)
