# @Author: Chenchen Zhu <czhu>
# @Date:   30-10-2019
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 2021-02-26


#! /bin/bash
# bam file as input
infile=$1
outfile=$(dirname $infile)/$(basename $infile .bam)
ncpu=30
genomefile=genome/hg38_SIRV.fa.genomeFile
scriptFolder=sw
## exclude secondary alignment and supplementary alignment
## 256 is secondary and 2048 is supplementary we don't take either
samtools view -@ $ncpu -u -F 2304 $infile | bedtools bamtobed -bed12 -cigar > "$outfile".bed
perl $scriptFolder/extract_splice_sites.pl "$outfile".bed

## we don't do sort, takes too much time and resource
## should not be needed
cat "$outfile".ie.bed  | bedtools genomecov -bg -i - -g $genomefile > "$outfile".ie.bedGraph &
cat "$outfile".ei.bed  | bedtools genomecov -bg -i - -g $genomefile > "$outfile".ei.bedGraph &

## sorting should not be necessary but just in case it's needed
## the output isn't real bed file
cat "$outfile".ie.bed  | bedtools genomecov -dz -i - -g $genomefile > "$outfile".ie.cov.bed &
cat "$outfile".ei.bed  | bedtools genomecov -dz -i - -g $genomefile > "$outfile".ei.cov.bed &

# ## the bedGraph files are for visualisation purposes
# cat "$outfile".ie.bed  | sort -k 1,1 | bedtools genomecov -bg -i - -g $genomefile > "$outfile".ie.bedGraph &
# cat "$outfile".ei.bed  | sort -k 1,1 | bedtools genomecov -bg -i - -g $genomefile > "$outfile".ei.bedGraph &
#
# ## sorting should not be necessary but just in case it's needed
# ## the output isn't real bed file
# cat "$outfile".ie.bed  | sort -k 1,1 | bedtools genomecov -dz -i - -g $genomefile > "$outfile".ie.cov.bed &
# cat "$outfile".ei.bed  | sort -k 1,1 | bedtools genomecov -dz -i - -g $genomefile > "$outfile".ei.cov.bed &

wait
#cat "$outfile".left.bed  | sort -k 1,1 | bedtools genomecov -bg -i - -g $genomefile > "$outfile".left.cov.bed
#cat "$outfile".right.bed  | sort -k 1,1 | bedtools genomecov -bg -i - -g $genomefile > "$outfile".right.cov.bed

#find . -name "*bedGraph" | xargs -n 1 -P 5 bgzip --force -i -@ 2
bgzip -f -i -@ $ncpu "$outfile".ie.bedGraph
bgzip -f -i -@ $ncpu "$outfile".ei.bedGraph

bgzip -f -i -@ $ncpu "$outfile".ie.bed
bgzip -f -i -@ $ncpu "$outfile".ei.bed
bgzip -f -i -@ $ncpu "$outfile".bed
bgzip -f -i -@ $ncpu "$outfile".ie.cov.bed
bgzip -f -i -@ $ncpu "$outfile".ei.cov.bed

wait
#bgzip -i -@ 5 "$outfile".bed "$outfile".ie.bed "$outfile".ei.bed "$outfile".ie.cov.bed "$outfile".ei.cov.bed
#find . -name "*bed" | xargs -n 1 -P 5 bgzip --force -i -@ 2
#bgzip -i "$outfile".left.cov.bed
#bgzip -i "$outfile".right.cov.bed
