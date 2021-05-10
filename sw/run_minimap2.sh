# @Author: Chenchen Zhu <czhu>
# @Date:   30-10-2019
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 2021-02-26

#!/bin/bash
## index
#faFile=GRCh37_sequin.fa.gz
index=hg38_SIRV.fa.gz.idx

## one time index
# minimap2 -t 30 -x splice -d GRCh37_sequin.fa.idx GRCh37_sequin.fa.gz

infile=$1
if [[ "$infile" != /* ]]; then
    infile=$PWD/$infile
fi
outfile=$2
if [ -z "$outfile" ]
  then
      outfolder=$(dirname $(dirname $infile))/alignment_minimap2
      fileBasename=$(basename $(dirname $(dirname $infile)) | cut -d. -f1)
      outfile=$outfolder/$fileBasename.bam
fi
mkdir -p $(dirname $outfile)

ncpu=$3
if [ -z "$ncpu" ]
  then
      ncpu=10
fi

flank=$4
if [ -z "$flank" ]
  then
      flank=yes
fi

BN=$(basename $outfile .bam)
RG="@RG\tID:"$BN"\tSM:"$BN


echo "input $infile out $outfile with read group $RG using $ncpu cpu"
minimap2 -K500m -t $ncpu --secondary=no -a -x splice $index $infile -R $RG --splice-flank=$flank | \
    samtools sort -@ $ncpu > $outfile && samtools index -@ $ncpu $outfile
echo "Done"
## N is secondary allowed
## splice  preset  (-k15 -w5 --splice -g2000 -G200k -A1 -B2 -O2,32 -E1,0 -z200 -ub --cost-non-gt-ag 5)
## do we  but cost on canonical splicing sites??
## gmap
