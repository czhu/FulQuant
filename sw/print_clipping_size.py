#!/usr/bin/python
# https://www.biostars.org/p/76119/
import sys;
import HTSeq;
import pysam;
import re;


#thisCigar = HTSeq.parse_cigar( read.cigarstring, read.reference_start, read.reference_name)
#thisCigar = HTSeq.parse_cigar( cigar, 100, "1s")
## provide this as a option
minClipSizeToReport = 10;

def extract_from_cigarOpList(cigarOpList, cigarOp="I"):
    myCigarOp = []
    for thisCigarOp in cigarOpList:
        if(thisCigarOp.type == cigarOp):
            myCigarOp.append(thisCigarOp)
    return myCigarOp

def is_spliced(thisRead):
    if(re.search('(\d+)N',read.cigarstring)):
        return True
    else:
        return False

def extract_exons(thisRead):
    if(is_spliced(thisRead)):
        thisCigar = HTSeq.parse_cigar( thisRead.cigarstring, thisRead.reference_start, thisRead.reference_name)
        myIntrons = extract_from_cigarOpList(thisCigar,cigarOp="N")
        exonStart = [thisRead.reference_start] + [x.ref_iv.end for x in myIntrons]
        exonEnd = [x.ref_iv.start for x in myIntrons] + [thisRead.reference_end]
    else:
        exonStart = thisRead.reference_start
        exonEnd  = thisRead.reference_end
    rv = []
    for idx in range(0,len(exonStart)):
        rv.append(HTSeq.GenomicInterval(thisRead.reference_name, exonStart[idx], exonEnd[idx]))
    return rv

def extract_introns(thisRead):
    if(is_spliced(thisRead)):
        thisCigar = HTSeq.parse_cigar( thisRead.cigarstring, thisRead.reference_start,
                                      thisRead.reference_name)
        myIntrons = extract_from_cigarOpList(thisCigar,cigarOp="N")
        rv = [x.ref_iv for x in myIntrons]
    else:
        rv = []
    return rv

def extract_clipping(thisRead):
    thisCigar = HTSeq.parse_cigar( thisRead.cigarstring, thisRead.reference_start,
                                  thisRead.reference_name)
    softClippings = extract_from_cigarOpList(thisCigar,cigarOp="S")
    hardClippings = extract_from_cigarOpList(thisCigar,cigarOp="H")
    clippings = softClippings + hardClippings
    return clippings

def read_to_iv(thisRead):
    iv = HTSeq.GenomicInterval(thisRead.reference_name, thisRead.reference_start,
                          thisRead.reference_end)
    return iv

bamFP = pysam.Samfile("-", "rb")
#bamFP = pysam.Samfile("./chrIS_subreads_no.bam", "rb")

for read in bamFP:
    ## read = bamFP.next()
    if( not( read.is_unmapped ) ):   #if it's mapped
        #print read.cigarstring
        cp = extract_clipping(read)

        if(len(cp)) :
            i = 0
            for thisCp in cp:
                i += 1
                if(thisCp.size>=minClipSizeToReport):
                    print read.query_name,"\t",read.reference_name,"\t", i,"\t",\
                        thisCp.size,"\t",thisCp.type
