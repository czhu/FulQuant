# @Author: Chenchen Zhu <czhu>
# @Date:   31-10-2019
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 2021-02-26
SCRIPTDIR="."
folder="."
indir=$folder/alignments
outdir1=$folder/alignments_no_supp/
outdir2=$folder/alignments_supp/
outdir3=$folder/alignments_unspliced
mkdir -p $outdir1
mkdir -p $outdir2
mkdir -p $outdir3
ncpu=20
### filter reads with supplementary
find $indir -name "*bam" | xargs -n 1 -P $ncpu -I% bash -c 'python sw/filter_reads_for_clustering.py % alignments_no_supp/$(basename %) alignments_supp/$(basename %) alignments_unspliced/$(basename %)'
