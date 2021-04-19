from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.preprocessing import label_binarize
from sklearn.externals import joblib
import matplotlib.pyplot as plt
from itertools import cycle
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
    '-m',
    '--model',
    type=str,
    default="RF",
    metavar='model',
    help='which model you want to used'
)

ap.add_argument(
    '-b',
    '--bam',
    type=str,
    required=True,
    metavar='bamfile',
    help='bam file that used to extract features'
)

ap.add_argument(
    '-s',
    '--soft',
    type=str,
    required=True,
    metavar='CNVcaller',
    help='which CNV caller you used'
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
    '-n',
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
    help='output directory'
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
    #if args.model not in ["RF", "GBC", "SVM"]:
    #    print("ERROR: model(-m) should be one of RF, GBC and SVM. For other models, you need to build your own")
    #    exit(1)
    #if args.soft not in ["Lumpy", "Manta", "Pindel", "Delly" , "breakdancer"]:
    #    print("ERROR: CNVcaller(-s) should be one of Lumpy, Manta, Pindel, Delly and breakdancer. For other CNVcaller, see CNV-P_training_main.py to build a new model")
    #    exit(1)
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
    print("Step1: feature extraction...")
    CNV_bed, outbed = rf.Read_bed(args.CNV_bed)
    bamfile = rf.Read_bam(args.bam)
    RG_bas = rf.Read_bas(args.basfile)
    outdir = args.output
    Feature_extract(CNV_bed, bamfile, RG_bas, outdir, args.Samplename)
    featfile = outdir + '/' + args.Samplename + ".feature.txt"
    x, info = rf.Read_feat(featfile)

    print("Step2: Loading mode and predicting.... ")
    mod = sys.path[0] + "/../model/" + args.soft + "." + args.model + ".train_model.m"
    if not os.path.exists(mod):
        print("ERROR: " + model + " not exists, please check the parameters: " + args.model + "(-m) " +args.soft + "(-s)" )
        exit(1)
    print(mod)
    yscore_dict={}
    Classifier = joblib.load(mod)
    y1_pre_score=Classifier.predict_proba(x)
    yscore_dict[mod]=y1_pre_score
    y_pre=Classifier.predict(x)
    yp=np.array(y_pre)
    ys=np.array(y1_pre_score[:,1])
    print(yp.shape,ys.shape)
    #y:true #yp:predict #ys:score
    pre_res=np.vstack((yp,ys)).T
#    out1=np.vstack(np.array(y_pre),np.array(y1_pre_score[:,1]))
#    out=np.vstack(out0,out1)
    outf= np.hstack((outbed,pre_res))
    header= np.array(["ChrID","start","end","length","CNV_type","class","probability_score"])
    outfinal= np.vstack((header,outf))

    print("Step3: Save the prediction results... ")
    np.savetxt(outdir +"/" + args.Samplename + ".pre.prop.txt",outfinal,fmt="%s",delimiter="\t")
    print("All steps done: the results show in " + args.output)