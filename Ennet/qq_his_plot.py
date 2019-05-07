#!/usr/bin/env python

from scipy.stats import probplot
import matplotlib as mpl
import networkx as nx
mpl.use('Agg')
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import escore


def qq_his_plot(G, output_name): 
    '''
    Draw qq plot and histogram

    @parameter G - graph
    @parameter output_name - output_basename

    '''
    from operator import itemgetter
    emp_p_n = escore.get_value_from_graph(G, 'gene', 'emp_p_n')
    p = escore.get_value_from_graph(G, 'gene', 'emp_p')
    p = sorted(p.items(), key=itemgetter(1), reverse=False)
    p_length = len(p)
    quantile = p_length // 5
    gene_list = [p[i][0] for i in [0, quantile-1, 2*quantile-1, 3*quantile-1, 4*quantile-1, p_length-1]]
    random_size = len(emp_p_n[gene_list[0]])
    plt.figure(1)
    plt.subplot(321)
    probplot(emp_p_n[gene_list[0]], plot=plt)
    plt.title(gene_list[0])
    plt.subplot(322)
    probplot(emp_p_n[gene_list[1]], plot=plt)
    plt.title(gene_list[1])
    plt.subplot(323)
    probplot(emp_p_n[gene_list[2]], plot=plt)
    plt.title(gene_list[2])
    plt.subplot(324)
    probplot(emp_p_n[gene_list[3]], plot=plt)
    plt.title(gene_list[3])
    plt.subplot(325)
    probplot(emp_p_n[gene_list[4]], plot=plt)
    plt.title(gene_list[4])
    plt.subplot(326)
    probplot(emp_p_n[gene_list[5]], plot=plt)
    plt.title(gene_list[5])
    plt.tight_layout()
    plt.savefig(output_name+"_qqplot.png")
    plt.close()
    
    plt.figure(2)
    plt.subplot(321)
    sns.distplot(emp_p_n[gene_list[0]], hist=False, rug=True);
    plt.title(gene_list[0])
    plt.subplot(322)
    sns.distplot(emp_p_n[gene_list[1]], hist=False, rug=True);
    plt.title(gene_list[1])
    plt.subplot(323)
    sns.distplot(emp_p_n[gene_list[2]], hist=False, rug=True);
    plt.title(gene_list[2])
    plt.subplot(324)
    sns.distplot(emp_p_n[gene_list[3]], hist=False, rug=True);
    plt.title(gene_list[3])
    plt.subplot(325)
    sns.distplot(emp_p_n[gene_list[4]], hist=False, rug=True);
    plt.title(gene_list[4])
    plt.subplot(326)
    sns.distplot(emp_p_n[gene_list[5]], hist=False, rug=True);
    plt.title(gene_list[5])
    plt.tight_layout()
    plt.savefig(output_name+"_histogram.png")
    plt.close()
    

if __name__=="__main__":
    pkl_file = open("/home/niuyw/Project/RegulatorySNP_170808/ennet_190429/20903b5/brain_network.pickle",'r')
    G = pickle.load(pkl_file)
    qqplot(G,"brain")
    pkl_file.close()
