# @Author: Chenchen Zhu <czhu>
# @Date:   30-10-2019
# @Email:  czhu@me.com
# @Last modified by:   czhu
# @Last modified time: 01-11-2019


#!/usr/bin/python
## filter reads with supplementary alignment into two files
## this is better than simple samtools view -F 2048
## samtools view only filters alignment but still retain the primanyly alignment
## i.e ie read has both primary and supplementary alignment, the primary alignment
## remains after samtools filter
## with this filter the reads are remove if any portion is a supplementary alignment
## this filter essentially identifies fusion alignment
## this script also adds RG to the qname
## 2019-05-09 add a few tags to qname for downstream analysis
## tags are defined in the manual of minimap2
##  tp |  A   | Type of aln: P/primary, S/secondary and I,i/inversion |
##  NM |  i   | Total number of mismatches and gaps in the alignment  |
##  ts |  A   | Transcript strand (splice mode only)                  |
"""
python filter.py - file_with_supp.bam file_without_supp.bam file_unspliced.bam
"""

import pysam;
import sys;
import os;
import re;

def is_spliced(thisRead):
    if(re.search('(\d+)N',read.cigarstring)):
        return True
    else:
        return False

infile = sys.argv[1]
outfile = sys.argv[2] ## good reads
outfile2 = sys.argv[3] ## reads with supplements
outfile3 = sys.argv[4] ## unspliced reads
#outfile = os.path.abspath(sys.argv[2])
#outfile2 = os.path.splitext(outfile)+"_withSupp.bam"

bamFP = pysam.AlignmentFile(infile, "rb")
otfile = pysam.AlignmentFile(outfile, "wb", template=bamFP)
otfile2 = pysam.AlignmentFile(outfile2, "wb", template=bamFP)
otfile3 = pysam.AlignmentFile(outfile3, "wb", template=bamFP)

for read in bamFP:
    if( not( read.is_unmapped ) ):
        if( read.has_tag("ts") ):
            ts = read.get_tag("ts")
        else:
            ts = "*"
        read.qname = "%s:%s:%s:%s:%s" % (read.get_tag("RG"),read.get_tag("tp"),
                                ts, read.get_tag("NM"), read.qname)

        if( read.has_tag("SA") ):
            otfile2.write(read)
        else:
            if( is_spliced(read)):
                # spliced reads
                otfile.write(read)
            else:
                # unspliced reads
                otfile3.write(read)
