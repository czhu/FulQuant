# @Author: Chenchen Zhu <czhu>
# @Date:   30-10-2019
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 2021-02-26


#! /bin/bash
# bam file as input
INFILE=$1
if [ "${INFILE:0:1}" != "/" ] ; then
    INFILE=$PWD/$INFILE
fi
OUTFILE=`dirname $INFILE`/`basename $INFILE tsv`clustered.bed
CURRENTDIR=$PWD

CHRLENFILE=genome/hg38_SIRV.fa.genomeFile

printf "input is $INFILE \noutput is $OUTFILE \n"

cd sw/TScluster
java -Xmx100G -cp .:jars/* TIFSeq.TScluster F1=$INFILE F2=$CHRLENFILE outF=$OUTFILE HFWINSIZE=5 DISTHRES=8 PTHRESHOLD=3 SUMHOLD=5
cd $CURRENTDIR
#javasol TIFSeq.TScluster F1=$INFILE F2=$CHRLENFILE outF=$OUTFILE HFWINSIZE=5 DISTHRES=11 PTHRESHOLD=2 SUMHOLD=3
#java -cp /g/steinmetz/czhu/project/pore/sw/TScluster/:jars/* TIFSeq.TScluster
# INFILE=/g/steinmetz/czhu/project/pore/analysis/20170731_wildtype/basecalling/test/testdir/reads_align.ei.reformatted.tsv
# INFILE=/g/steinmetz/czhu/project/pore/analysis/20170731_wildtype/basecalling/test/testdir/reads_align.ie.reformatted.tsv
#
# OUTFILE=`dirname $INFILE`/`basename $INFILE tsv`clustered.bed
# CHRLENFILE=/g/steinmetz/czhu/project/pore/data/genome/GRCh38.82/GRCh38_ERCC_sequin.fa.genomeFile
# cd /g/steinmetz/czhu/project/pore/sw/TScluster
# java -cp .:./src/:jars/* TIFSeq.TScluster F1=$INFILE F2=$CHRLENFILE outF=$OUTFILE HFWINSIZE=5 DISTHRES=10 PTHRESHOLD=2 SUMHOLD=3
