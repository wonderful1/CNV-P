# -*- coding: utf-8 -*-
#sklearn
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_blobs
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm
from sklearn.externals import joblib #save mode
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc ,precision_recall_curve

import sys
import math
import os
import matplotlib as mpl
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy import interp
import matplotlib.pyplot as plt
from itertools import cycle
import argparse

import Read_File as rf
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
    '-s',
    '--soft',
    type=str,
    required=True,
    metavar='CNVcaller',
    help='which CNV caller you used'
)

ap.add_argument(
    '-fea',
    '--features',
    type=str,
    required=True,
    metavar='featuresfile',
    help='file that provide traing features'
)

ap.add_argument(
    '-lab',
    '--labelfile',
    type=str,
    required=True,
    metavar='labelfile',
    help='file that provide CNV label, true CNVs labeled as 1,false CNVs labeled as 0, The order should corresponds to CNV_bed file(-bed/--CNV_bed) one to one'
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

def Make_Model(model):
    n_features = 16  # for RF & ETC
    if model == "RF":
        Classifier = RandomForestClassifier(n_estimators=1000, max_features=n_features, max_depth=None, min_samples_split=3, bootstrap=True)
    elif model == "SVM":
        # SVM-rbf
        Classifier = svm.SVC(kernel='rbf', C=0.8, gamma='auto', probability=True, random_state=1)
    elif model == "GBC":
        Classifier = GradientBoostingClassifier(n_estimators=1000, max_features=n_features, random_state=1, min_samples_split=3)
    else:
        print("no such model; model must be one of RF, GBC and SVM")
        exit(1)
    return Classifier

if __name__ == "__main__":
    print("running: check the parameter...")
    if args.model not in ["RF", "GBC", "SVM"]:
        print("ERROR: model(-m) should be one of RF, GBC and SVM. For other models, you need to change the code")
        exit(1)
    if not os.path.exists(args.output):
        print("ERROR: " + args.output + " not exists, please check the path")
        exit(1)
    if not os.path.exists(args.features):
        print("ERROR: " + args.features + " not exists, please check the path")
        exit(1)
    if not os.path.exists(args.labelfile):
        print("ERROR: " + args.labelfile + " not exists, please check the path")
        exit(1)

    print("Check finished, start running...")
    print("Step1: load features and labels...")
    outdir = args.output
    x, info = rf.Read_feat(args.features)
    y = rf.Read_label(args.labelfile)

    print("Step2: training a model; perform 10fold-cross_validation and draw ROC...")
    KF = KFold(n_splits=10)
    i = 0
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    for train_index, test_index in KF.split(x):
        Classifier = Make_Model(args.model)
        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        Classifier.fit(x_train, y_train)
        y1_pre_score = Classifier.predict_proba(x_test)
        n_class = 2
        lb_y_test = label_binarize(y_test, classes=[str(i) for i in range(n_class)])
        fpr, tpr, thresholds = roc_curve(lb_y_test, y1_pre_score[:, 1])
        tprs.append(interp(mean_fpr, fpr, tpr))
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        # plt.plot(fpr,tpr,lw=1 ,color="red",label='ROC fold %d(area=%0.2f)'%(i,roc_auc))
        plt.plot(fpr, tpr, lw=1, color="red", alpha=0.4)

        i += 1
    plt.plot([0, 1], [0, 1], 'k--', lw=2, alpha=.7)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[0] = 0
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(tprs, axis=0)
    plt.plot(mean_fpr, mean_tpr, color='b', label=r'Mean ROC (area=%0.2f)' % mean_auc, lw=2, alpha=.8)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    title="CNV-P_"+ args.soft + "_" + args.model
    plt.title(title)
    plt.legend(loc="lower right")
    plt.savefig(args.output + '/' + title + '_Classifier.ROC.png')
    plt.savefig(args.output + '/' + title + '_Classifier.ROC.pdf')

    print("Step3: save the final model...")
    Final_Classifier = Make_Model(args.model)
    Final_Classifier.fit(x, y)
    joblib.dump(Final_Classifier, args.output + "/" + args.soft + "." + args.model + ".train_model.m")

    print("All steps done: the results see " + args.output )