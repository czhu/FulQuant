# @Author: Chenchen Zhu <czhu>
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 2021-02-26
ncpu=30
PROG=sw/run_minimap2.sh
folder="."
indir=$folder/fastq_trimmed/
outdir=$folder/alignments/
mkdir -p $outdir

for infile in $(find $indir -name "*fastq.gz" -not -name "*dicarded*"); do
    $PROG $infile $outdir/$(basename $infile .fastq.gz).bam $ncpu
done
