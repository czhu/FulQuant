# @Author: Chenchen Zhu <czhu>
# @Date:   29-10-2019
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 2021-02-26


## merge bam files for this dataset both long and read reads
folder="."
indir=$folder/alignments_no_supp/
outdir=$folder/combined/
mkdir -p $outdir
outfile=$outdir/ont_combined_no_supp.bam
samtools merge -u -@ 20 - $(ls $indir/*bam) | samtools sort -m 4G -@ 20 > $outfile && samtools index -@ 20 $outfile

ln -s $outfile $outdir/all.bam
ln -s $outfile.bai $outdir/all.bam.bai
