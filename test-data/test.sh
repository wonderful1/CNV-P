python /data/lush-dev/wangtaifu/learning/CNV_P/CNV-P/src/CNV-P_predict_main.py -m RF -b test-data/HG002.test.bam -s Lumpy -bed test-data/HG002.Lumpy.fil.mer.bed -bas test-data/HG002.bam.bas -sam HG002 -o test-data/out/

python /data/lush-dev/wangtaifu/learning/CNV_P/CNV-P/src/CNV-P_featureExtract_main.py -b test-data/HG002.test.bam -bed test-data/HG002.Lumpy.fil.mer.bed -bas test-data/HG002.bam.bas -sam HG002 -o test-data/out/

python /data/lush-dev/wangtaifu/learning/CNV_P/CNV-P/src/CNV-P_training_main.py -s Lumpy -fea test-data/HG002.Lumpy.chr1.feature.txt -lab test-data/HG002.Lumpy.chr1.label.txt -o test-data/out/
