import pysam as ps
import os
import numpy as np

def Read_bed(bed):
     if not os.path.isfile(bed):
          print('Bedfile does not exist at location provided, please check (tried "{}")'.format(bed))
          exit(2)
     print("Loading bed file: {} ...".format(bed))
     bed_f = open(bed)
     CNV_bed=[]
     for line in bed_f:
          line=line.strip("\n")
          line_list=line.split("\t")
          CNV_bed.append(line_list)
     outbed = np.array(CNV_bed)
     bed_f.close()
     return CNV_bed, outbed

def Read_bas(bas):
     if not os.path.isfile(bas):
          print('Bedfile does not exist at location provided, please check (tried "{}")'.format(bas))
          exit(2)
     print("Loading bas file: {} ...".format(bas))
     bas_f = open(bas)
     RG_bas={}
     for line in bas_f:
          line=line.strip("\n")
          line_list=line.split("\t")
          RG_bas[line_list[6]]=line_list[17]+"\t"+line_list[18]+"\t"+line_list[-1]
     bas_f.close()
     return RG_bas

def Read_bam(bam):
     print("Loading bam file: {} ...".format(bam))
     bamfile = ps.AlignmentFile(bam, 'rb')
     try:
          bamfile.check_index()
          print('Index file found - {}'.format(bam + '.bai'))
     except ValueError:
          print('Index file not found (should be named "{}.bai" and located in the same directory)'.format(bam))
     return bamfile


def Read_feat(feat_f):
     feat_f = open(feat_f, "r")
     data_list1 = []
     col4_li = []
     for line0 in feat_f:
          line0 = line0.strip("\n")
          num = line0.split("\t")
          if num[0] == "space": feat_labels = num[3:20];continue
          col4 = num[0:4]
          num = num[3:20]
          data_list1.append(num)
          col4_li.append(col4)
     x1 = np.array(data_list1)
     info = np.array(col4_li)
     print(x1.shape)
     print(info.shape)
     return x1, info

def Read_label(lab_f):
     lab_f= open(lab_f, "r")
     ylist = []
     for line1 in lab_f:
          line1 = line1.strip("\n")
          # line1=line1.split("\t")
          # ylist=line1
          ylist.append(line1)
     # ylist.pop()
     y1 = np.array(ylist)
     print(y1.shape)
     return y1
