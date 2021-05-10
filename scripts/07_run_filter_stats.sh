# @Author: Chenchen Zhu <czhu>
# @Date:   29-10-2019
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 2021-02-26


## run various tables for filtering later
scriptFolder=sw
folder="."
infile=$folder/combined/all.bam

infolder=$(dirname $infile)
cd $infolder
ncpu=15

samtools view -@ $ncpu -u $infile | python $scriptFolder/count_insert_size_near_splice_site.py \
    | pigz -c -p $ncpu > insert_profile_near_ss.gz &

samtools view -@ $ncpu -u $infile | python $scriptFolder/print_clipping_size.py \
    | pigz -c -p $ncpu > clip_profile_at_ends.gz &

## run polyA mask scan
