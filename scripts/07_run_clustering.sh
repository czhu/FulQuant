# @Author: Chenchen Zhu <czhu>
# @Date:   29-10-2019
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 2021-02-26


## first merge relevant bam files into one
SCRIPTDIR=sw/
folder="."

cd $folder/combined

bash $SCRIPTDIR/run_clustering.sh all.bam
