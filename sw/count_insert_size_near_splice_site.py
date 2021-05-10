#!/usr/bin/python
# https://www.biostars.org/p/76119/
import sys;
import pysam;
import HTSeq;
import re;

#bamFile    = sys.argv[1];
bamFP = pysam.Samfile("-", "rb")
#bamFP = pysam.Samfile("./chrIS_subreads_no.bam", "rb")
halfWindowSize = 15;
## FIXME provide this as a option
minInsertToReport = 3

def extract_from_cigarOpList(cigarOpList, cigarOp="I"):
    myCigarOp = []
    for thisCigarOp in cigarOpList:
        if(thisCigarOp.type == cigarOp):
            myCigarOp.append(thisCigarOp)
    return myCigarOp

def calc_gi_ovlp_count(gi,featList):
    ovlpFeatLen = []
    for thisFeat in featList:
        if(gi.contains(thisFeat.ref_iv)):
            ovlpFeatLen.append(thisFeat.size)
    return sum(ovlpFeatLen)

for read in bamFP:
    ## read = bamFP.next()
    if( not( read.is_unmapped ) ):   #if it's mapped
        #print read.cigarstring
        if(re.search('(\d+)N',read.cigarstring)):  # only spliced reads
            thisCigar = HTSeq.parse_cigar( read.cigarstring, read.reference_start, read.reference_name)

            myIntrons = extract_from_cigarOpList(thisCigar,cigarOp="N")
            myInserts = extract_from_cigarOpList(thisCigar,cigarOp="I")

            if(len(myInserts)) :
                i = 0
                for thisIntron in myIntrons:
                    i += 1
                    thisLeft=HTSeq.GenomicInterval(
                        chrom=read.reference_name,
                        start=thisIntron.ref_iv.start-halfWindowSize,
                        end=thisIntron.ref_iv.start)

                    thisRight=HTSeq.GenomicInterval(
                        chrom=read.reference_name,
                        start=thisIntron.ref_iv.end,
                        end=thisIntron.ref_iv.end + halfWindowSize)

                    insertSizeLeft = calc_gi_ovlp_count(thisLeft,myInserts)
                    insertSizeRight = calc_gi_ovlp_count(thisRight,myInserts)

                    totalSize = insertSizeLeft + insertSizeRight

                    if( totalSize >=minInsertToReport) :
                        print read.query_name,"\t",read.reference_name,"\t", i,"\t",\
                            thisIntron.ref_iv.start,"\t",thisIntron.ref_iv.end,"\t", \
                            insertSizeLeft,"\t", insertSizeRight, "\t", totalSize
