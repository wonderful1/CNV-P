import numpy as np
import pandas as pd
import math

def Get_Depth(bamfile,chrom,pos_l,pos_r,RG_bas_dict,sample):
    cova=RG_bas_dict[sample].split("\t")[-1]
    if pos_l == pos_r:
        ave_depth=0
        gc_cont=0
    else:
        depth = bamfile.count_coverage(chrom,pos_l,pos_r)
        np_depth=np.array(list(depth))
        sum_depth = np_depth.sum(axis=0)
        ave_depth=sum_depth.sum()/(pos_r-pos_l+1)
        if ave_depth > 0: ave_depth=math.log2(ave_depth/float(cova))
        if np_depth.sum()==0 : gc_cont=0
        else:gc_cont=(np_depth[1].sum()+np_depth[2].sum())/np_depth.sum()
    return ave_depth,gc_cont
#Get_Depth(bamfile,"19",60006,60100)

def Get_Feature(bamfile,chrom,pos_l,pos_r,RG_bas_dict,sample):
    '''1)proper_pair  2)mapping_quality  3)insert size  4)split read '''
    SR=0
    QUA=0
    PAIR=0
    PEM=0
    cova=RG_bas_dict[sample].split("\t")[-1]
    for read in bamfile.fetch(chrom, pos_l, pos_r):
        if read.cigarstring != None:
            split_read=(1 if(read.get_cigar_stats()[0][4]>= 2) else 0)
            mean,stdev,co=RG_bas_dict[sample].split("\t")
            mean=float(mean)
            stdev=float(stdev)
            ab_ins_size=(1 if ((abs(read.next_reference_start-read.get_reference_positions()[1]+read.query_alignment_length)-mean) >= 3*stdev) else 0)
            pair_is=(0 if(read.is_proper_pair) else 1)
            pe_is=ab_ins_size
            qua_is=(0 if(read.mapping_quality> 10) else 1)
            #split read
            SR+=split_read
            #quality
            QUA+=qua_is
            #paired
            PAIR+=pair_is
            #pair end
            PEM+=pe_is
    SR=SR/float(cova)
    QUA=QUA/float(cova)
    PAIR=PAIR/float(cova)
    PEM=PEM/float(cova)
    return SR,QUA,PAIR,PEM
