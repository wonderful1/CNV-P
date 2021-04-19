# CNV-P: A machine-learning framework for filtering copy number variations
CNV-PG is a novel and post–processing approach for CNV filtering.  

## Prerequisites:
[python3](https://www.python.org/)  
[sklearn](https://pypi.org/project/sklearn/)  
[matplotlib](https://pypi.org/project/matplotlib/)  
[pysam](https://pypi.org/project/pysam/)  
[pandas](https://pypi.org/project/pandas/)  
[numpy](https://pypi.org/project/numpy/)  

## Install by conda
conda install python=3.7.0  
conda install -c anaconda scikit-learn=0.21  
conda install -c conda-forge matplotlib  
conda install -c bioconda pysam  
conda install -c anaconda pandas  
conda install -c anaconda numpy  

## Getting started
### 1. CNV Predicting
run the "python  script/CNV-P_predict_main.py -h" to see the USAGE;  
```
usage: CNV-P_predict_main.py [-h] [-m model] -b bamfile -s CNVcaller -bed
                             BEDfile -bas basfile [-sam Samplename]
                             [-o outdir]
optional arguments:
  -h, --help            show this help message and exit
  -m model, --model model
                        which model you want to used
  -b bamfile, --bam bamfile
                        file provides features
  -s CNVcaller, --soft CNVcaller
                        which CNV caller you used
  -bed bedfile, --CNV_bed bedfile
                        the format of input file
  -bas basfile, --basfile basfile
                        file that provide mean insert-size and sequencing
                        depth
  -n Samplename, --Samplename Samplename
                        Samplename that use as prefix of result
  -o outdir, --output outdir
                         output directory
```
#### 1.1 input parameters
**model**:  should be one of RF (Random Forest), GBC (Gradient Boosting classifier) and SVM (Support Vector Machine)  
**CNVcaller**: Lumpy, Manta, Pindel, Delly and breakdancer is currently supported,  for Other software needs to be pre-trained(see 2.training for other CNV callers)  
**bamfile**: BAM file should generated by a read aligner that supports partial read alignments, such as BWA-MEM  
**bedfile**: This file should be 5 Columns: chromsome, start, end, length of CNV, type of CNV (DUP:1,DEL:0) 
for example (test_data/HG002.Lumpy.fil.mer.bed):  
```
chr19	350768	351961	1194	1
chr19	434243	434587	345	1
chr19	566222	569347	3126	0
chr19	878739	879857	1119	1
chr19	1182660	1183097	438	0
chr19	1572816	1573149	334	0
chr19	2033040	2033182	143	0
chr19	2713161	2714159	999	0
```
**basfile**: this file should be 4 columns: Samplename, median value of insert size, standard deviation of insert size, coverage
for example (test_data/HG002.bam.bas):
```
Samplename	median_insert_size	insert_size_median_sd	coverage
HG002	568.177944	163.819637	35.41
```

#### 1.2 output
**samplename.feature.txt**: Extracted feature matrix.  
**samplename.pre.prop.txt**: The prediction result and probability score. Including 7 columns:
```
ChrID: Chromosome (e.g. chr3, chrY)
start: Start coordinate on the chromosome 
end: End coordinate on the chromosome
length: length of CNV
CNV_type: type of CNV (DUP:1,DEL:0)
class: predicting results (true CNV：1 ,false CNV: 0)
probability_score: Probability of this CNV to be true
```


#### 1.3 running example
```
python  script/CNV-P_predict_main.py  -m RF -b Test_data/HG002.test.bam -s Lumpy -bed Test_data/HG002.Lumpy.fil.mer.bed -bas Test_data/HG002.bam.bas -sam HG002 -o Test_data/out/
```
  
  
### 2. training for other CNV callers
For  training a model for other CNV callers, use 'CNV-P_featureExtract_main.py' to perform features extraction:  
```
python script/CNV-P_featureExtract_main.py -b test-data/HG002.test.bam -bed test-data/HG002.Lumpy.fil.mer.bed -bas test-data/HG002.bam.bas -sam HG002 -o test-data/out/
```
then，run the "script/CNV-P_training_main.py" to train a model  
run " python script/CNV-P_training_main.py -h " to see the USAGE;  
```
usage: CNV-P_training_main.py [-h] [-m model] -s CNVcaller -fea featuresfile
                              -lab labelfile [-o outdir]
optional arguments:
  -h, --help            show this help message and exit
  -m model, --model model
                        which model you want to used
  -s CNVcaller, --soft CNVcaller
                        which CNV caller you used
  -fea featuresfile, --features featuresfile
                        file that provide traing features
  -lab labelfile, --labelfile labelfile
                        file that provide CNV label, true CNVs labeled as
                        1,false CNVs labeled as 0, The order should
                        corresponds to CNV_bed file(-bed/--CNV_bed) one to one
  -o outdir, --output outdir
                         output directory
```
#### 2.1 input parameters
**featuresfile**:  file that provide traing features, results from 'CNV-P_featureExtract_main.py'  
**labelfile**: one column, true CNVs labeled as 1,false CNVs labeled as 0  

#### 2.2 outputs:
**CNVcaller.model.train_model.m**: the classifier you trained  
**CNV-P_CNVcaller_model_Classifier.ROC.pdf, CNV-P_CNVcaller_model_Classifier.ROC.png**: the ROC of 10fold-cross_validation  

#### 2.3 running example
```
python script/CNV-P_training_main.py -s Lumpy -fea test-data/HG002.Lumpy.chr1.feature.txt -lab test-data/HG002.Lumpy.chr1.label.txt -o test-data/out/
```

Please help us improve CNV-P by reporting bugs or ideas on how to make things better.  

