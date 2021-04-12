import numpy as np
import sys
import argparse
import Read_File as rf
import Get_Feature as gf
import os
import warnings
warnings.filterwarnings("ignore")

#Defined the USAGE
ap = argparse.ArgumentParser(description='CNV-P: A machine-learning framework for filtering copy number variations')

ap.add_argument(
    '-b',
    '--bam',
    type=str,
    required=True,
    metavar='bamfile',
    help='bam file that used to extract features'
)

ap.add_argument(
    '-bed',
    '--CNV_bed',
    type=str,
    required=True,
    metavar='bedfile',
    help='bed file that provide the CNV calls'
)

ap.add_argument(
    '-bas',
    '--basfile',
    type=str,
    required=True,
    metavar='basfile',
    help='file that provide mean insert-size and sequencing depth'
)

ap.add_argument(
    '-sam',
    '--Samplename',
    type=str,
    default='CNV-P',
    metavar='Samplename',
    help='Samplename that use as prefix of result'
)

ap.add_argument(
    '-o',
    '--output',
    type=str,
    default='./',
    metavar='outdir',
    help='the directory into which output files will go'
)
args = ap.parse_args()

#extract traning features
def Feature_extract(CNV_bed,sam_file,RG_bas_dict,outdir,sample):
     f = open(outdir+ sample + ".feature.txt",'w')

     bed_len = len(CNV_bed)
     print("bed_len %d"%bed_len)
     for i in range(bed_len):
          print("i = %d"%i)
          depth,gc_cont=gf.Get_Depth(sam_file,CNV_bed[i][0],int(CNV_bed[i][1]),int(CNV_bed[i][2]),RG_bas_dict,sample)
          depth_l,gc_cont_l=gf.Get_Depth(sam_file,CNV_bed[i][0],int(CNV_bed[i][1])-1000,int(CNV_bed[i][1]),RG_bas_dict,sample)
          depth_r,gc_cont_r=gf.Get_Depth(sam_file,CNV_bed[i][0],int(CNV_bed[i][2]),int(CNV_bed[i][2])+1000,RG_bas_dict,sample)
          sr_l,qua_l,pair_l,pem_l=gf.Get_Feature(sam_file,CNV_bed[i][0],int(CNV_bed[i][1])-500,int(CNV_bed[i][1])+500,RG_bas_dict,sample)
          sr_r,qua_r,pair_r,pem_r=gf.Get_Feature(sam_file,CNV_bed[i][0],int(CNV_bed[i][2])-500,int(CNV_bed[i][2])+500,RG_bas_dict,sample)
          out="\t".join(CNV_bed[i])
          fea="\t".join(map(str,[depth,gc_cont,depth_l,gc_cont_l,depth_r,gc_cont_r,sr_l,qua_l,pair_l,pem_l,sr_r,qua_r,pair_r,pem_r]))
          f.write(out+"\t"+fea+"\n")

if __name__ == "__main__":
    print("running: check the parameter...")
    if not os.path.exists(args.bam):
        print("ERROR: " + args.bam + " not exists, please check the path")
        exit(1)
    if not os.path.exists(args.CNV_bed):
        print("ERROR: " + args.CNV_bed + " not exists, please check the path")
        exit(1)
    if not os.path.exists(args.basfile):
        print("ERROR: " + args.basfile + " not exists, please check the path")
        exit(1)
    if not os.path.exists(args.output):
        print("ERROR: " + args.output + " not exists, please check the path")
        exit(1)
    print("Check finished, start running...")
    print("Feature extraction...")
    CNV_bed, outbed = rf.Read_bed(args.CNV_bed)
    bamfile = rf.Read_bam(args.bam)
    RG_bas = rf.Read_bas(args.basfile)
    outdir = args.output
    Feature_extract(CNV_bed, bamfile, RG_bas, outdir, args.Samplename)
    featfile = outdir + '/' + args.Samplename + ".feature.txt"