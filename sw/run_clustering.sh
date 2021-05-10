# @Author: Chenchen Zhu <czhu>
# @Date:   30-10-2019
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 2021-02-26

#! /bin/bash
## master script to run the whole clustering process
## bam file as input

infile=$1
ncpu=20
SCRIPTDIR=sw/
outfile=$(dirname $infile)/$(basename $infile .bam)

## extract the splice sites
bash $SCRIPTDIR/extract_read_feature.sh $infile

## reformatt
Rscript $SCRIPTDIR/bedGraph2tsv.R "$outfile".ie.cov.bed.gz &
Rscript $SCRIPTDIR/bedGraph2tsv.R "$outfile".ei.cov.bed.gz &

wait

## clustering for two types splice sites
bash $SCRIPTDIR/ss_clustering.sh "$outfile".ei.cov.bed.reformatted.tsv &
bash $SCRIPTDIR/ss_clustering.sh "$outfile".ie.cov.bed.reformatted.tsv &

wait

bgzip -f -i -@ $ncpu "$outfile".ei.cov.bed.reformatted.clustered.bed &
bgzip -f -i -@ $ncpu "$outfile".ie.cov.bed.reformatted.clustered.bed &

wait

## get the clusters
Rscript $SCRIPTDIR/clustering_walking.R "$outfile".bed.gz "$outfile".ei.cov.bed.reformatted.clustered.bed.gz "$outfile".ie.cov.bed.reformatted.clustered.bed.gz


#find . -name "*bed" | xargs -n 1 -P 5 bgzip -i -@ 2
