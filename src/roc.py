#!/usr/bin/env python
# ./roc.py inputfile
import numpy as np
import sys
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, f1_score, auc, average_precision_score
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot

def main():
    positive_gene_file=open('/home/tengxy/work/ennet/positive_set/positive_genes.txt','r')
    result_file=open(sys.argv[1],'r')
    base_name='.'.join(sys.argv[1].split('/')[-1].split('.')[:-1])
    positive_set=[]
    print("Inputfile:%s" % sys.argv[1])
    for line in positive_gene_file:
        line=line.strip()
        positive_set.append(line)
    
    y_true=[]
    y_score=[]
    
    for line in result_file:
        line=line.strip()
        li=line.split('\t')
        if li[0] in positive_set:
            y_true.append(1)
        else:
            y_true.append(0)
        y_score.append(float(li[6]))
    
    y_true=np.array(y_true)
    y_score=np.array(y_score)
    fpr, tpr, thresholds_1 = roc_curve(y_true, y_score, pos_label = 1)
    precision, recall, thresholds_2 = precision_recall_curve(y_true, y_score)
    ap=average_precision_score(y_true, y_score)
    print('ROC_AUC: %.3f' % roc_auc_score(y_true, y_score))
    print('PRC_AUC: %.3f ap=%.3f' % (auc(recall,precision),ap))    

    pyplot.plot([0, 1], [0, 1], linestyle='--')
    pyplot.plot(fpr, tpr, marker='.')
    pyplot.xlabel('False Positive Rate')
    pyplot.ylabel('True Positive Rate')
    pyplot.title('Receiver operating characteristic(ROC)')
    pyplot.legend(loc="lower right")
    pyplot.savefig(base_name+".roc.png")
    pyplot.close()

    pyplot.plot([0, 1], [0.5, 0.5], linestyle='--')
    pyplot.plot(recall, precision, marker='.')
    pyplot.xlabel('Recall')
    pyplot.ylabel('Precision')
    pyplot.title('Precision-Recall curve: AP={0:0.3f}'.format(ap))
    pyplot.savefig(base_name+".prc.png")
    pyplot.close()

    result_file.close()
    positive_gene_file.close()


if __name__=='__main__':
    main()
