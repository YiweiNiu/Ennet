#!/usr/bin/env python

import networkx as nx
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pydotplus

def plot_hotspot_network(G, saveFile, layout = "graphviz",   #graphviz, circular, random, spectral, spring, shell
                          node_size = 1200,
                          width= 1 ,
                          font_weight = "bold",
                          with_labels = "True",
                          node_color = 'yellow',
                          node_alpha = 0.6,
                          edge_alpha = 0.4,
                          edge_width = 1,
                          edge_color = "k"):
    """
    Plot a undirected hotspot network.
    """
    node_options = {
        'node_size': node_size,
        'width': width,
        'font_weight': font_weight,
        'with_labels': with_labels,
        'node_color': node_color,
        'node_alpha': node_alpha
    }

    edge_options1 = {
        'edge_alpha': edge_alpha,
        'edge_width': edge_width,
        'edge_color': edge_color
    }

    if layout == "circular":
        pos=nx.circular_layout(G)
    elif layout == "random":
        pos=nx.random_layout(G)
    elif layout == "spectral":
        pos=nx.spectral_layout(G)
    elif layout == "graphviz":
        pos=nx.nx_pydot.graphviz_layout(G)
    elif layout == "spring":
        pos=nx.spring_layout(G, iterations=20, scale=2)
    elif layout == "shell":
        pos=nx.shell_layout(G)
    else:
        try:
            pos=nx.graphviz_layout(G)
        except:
            pos=nx.spring_layout(G,iterations=20)

    plt.rcParams['text.usetex'] = False
    if G.number_of_nodes() > 100:
        plt.figure(figsize=(30,30))
    elif G.number_of_nodes() < 100 and G.number_of_nodes() > 50:
        plt.figure(figsize=(20,20))
    else:
        plt.figure(figsize=(10,10))

    nx.draw_networkx_nodes(G, pos, **node_options)
    nx.draw_networkx_edges(G, pos, node_size=0, **edge_options1)
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=5, edge_color='m')
    nx.draw_networkx_labels(G, pos, fontsize=14)

    plt.axis('off')
    plt.savefig(saveFile, dpi=275)
    plt.show()

def plot_hotspot_network_known(G, saveFile, knownGenes, layout = "graphviz",   #graphviz, circular, random, spectral, spring, shell
                          node_size = 1200,
                          width= 1 ,
                          font_weight = "bold",
                          with_labels = "True",
                          node_alpha = 0.6,
                          edge_alpha = 0.4,
                          edge_width = 1,
                          edge_color = "k"):
    """
    Plot a undirected hotspot network with known genes.
    """
    node_options = {
        'node_size': node_size,
        'width': width,
        'font_weight': font_weight,
        'with_labels': with_labels,
        'node_alpha': node_alpha
    }

    node_body_options = {
        'node_size': 2300,
        'font_weight': font_weight,
        'with_labels': "False"
    }

    node_promoter_options = {
        'node_size': 2000,
        'font_weight': font_weight,
        'with_labels': "False"
    }

    edge_options1 = {
        'edge_alpha': edge_alpha,
        'edge_width': edge_width,
        'edge_color': edge_color
    }
    if layout == "circular":
        pos=nx.circular_layout(G)
    elif layout == "random":
        pos=nx.random_layout(G)
    elif layout == "spectral":
        pos=nx.spectral_layout(G)
    elif layout == "graphviz":
        pos=nx.nx_pydot.graphviz_layout(G)
    elif layout == "spring":
        pos=nx.spring_layout(G, iterations=20, scale=2)
    elif layout == "shell":
        pos=nx.shell_layout(G)
    else:
        try:
            pos=nx.graphviz_layout(G)
        except:
            pos=nx.spring_layout(G,iterations=20)

    with open(knownGenes) as known:
        arrs  = [l.split() for l in known]
    knownDict = dict((arr[0], arr[1]) for arr in arrs)

    nodeList_cancer_gene_known=[]
    nodeList_cancer_gene_candidate=[]
    nodeList_new=[]

    for node in G.nodes():
        if (node in knownDict) and knownDict[node] == "known":
            nodeList_cancer_gene_known.append(node)
        elif (node in knownDict) and knownDict[node] == "candidate":
            nodeList_cancer_gene_candidate.append(node)
        else:
            nodeList_new.append(node)

    plt.rcParams['text.usetex'] = False
    if G.number_of_nodes() > 100:
        plt.figure(figsize=(30,30))
    elif G.number_of_nodes() < 100 and G.number_of_nodes() > 50:
        plt.figure(figsize=(20,20))
    else:
        plt.figure(figsize=(10,10))

    if nodeList_cancer_gene_known:
        #nx.draw_networkx_nodes(G, pos, node_shape="8", nodelist=nodeList_cancer_gene_known[1:5], alpha=0.8, node_color='r', linewidths=0, **node_body_options)
        nx.draw_networkx_nodes(G, pos, nodelist=nodeList_cancer_gene_known, linewidths=0.5, alpha=0.8, node_color='red', **node_options)
    if nodeList_cancer_gene_candidate:
        nx.draw_networkx_nodes(G, pos, nodelist=nodeList_cancer_gene_candidate, node_color='yellow', **node_options)
    if nodeList_new:
        #nx.draw_networkx_nodes(G, pos, node_shape="8", nodelist=nodeList_new[1:5], alpha=0.8, node_color='r', linewidths=0, **node_body_options)
        nx.draw_networkx_nodes(G, pos, nodelist=nodeList_new, linewidths=0.5, alpha=0.8, node_color='w', **node_options)

    nx.draw_networkx_edges(G, pos, node_size=0, **edge_options1)
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=5, edge_color='m')
    nx.draw_networkx_labels(G, pos, fontsize=14)

    #x,y=pos["JAK2"]
    #print x,y
    #sizes = [60,30,10]
    #colors = ['red','yellowgreen','lightskyblue']
    #labels = ["a","b","c"]
    #explode = (0.05,0,0)
    #matplotlib.patches.Wedge((10,10), 5.8, 90, 360, width=10.5)
    #plt.text(10,20,s='some text', bbox=dict(facecolor='red', alpha=0.5),horizontalalignment='center')
    #plt.pie(sizes, explode=None, labels=None, colors=colors, radius=0.2, center=(100, 100))

    plt.axis('off')
    plt.savefig(saveFile+".png", dpi=275)
    plt.show()

def draw_ppi(G,known_cancer_gene_file,candidate_cancer_gene_file,network_out_file_name):
    '''
    '''
    cancer_gene_known_file=open(known_cancer_gene_file,"r")
    cancer_gene_candidate_file=open(candidate_cancer_gene_file,"r")
    cancer_gene_known=[]
    for line in cancer_gene_known_file:
        cancer_gene_known.append(line[0:-1])
    cancer_gene_candidate=[]
    for line in cancer_gene_candidate_file:
        cancer_gene_candidate.append(line[0:-1].upper())

    nodeList_cancer_gene_known=[]
    nodeList_cancer_gene_candidate=[]
    nodeList_new=[]
    for node in G.nodes():
        if node in cancer_gene_known:
            nodeList_cancer_gene_known.append(node)
        elif node in cancer_gene_candidate:
            nodeList_cancer_gene_candidate.append(node)
        else:
            nodeList_new.append(node)

    try:
        pos=nx.graphviz_layout(G)
    except:
        pos=nx.spring_layout(G,iterations=20)

    plt.rcParams['text.usetex'] = False
    if G.number_of_nodes() > 100:
        plt.figure(figsize=(30,30))
    elif G.number_of_nodes() < 100 and G.number_of_nodes() > 50:
        plt.figure(figsize=(20,20))
    else:
        plt.figure(figsize=(10,10))

    if nodeList_cancer_gene_known:
        nx.draw_networkx_nodes(G,pos,node_size=1200,nodelist=nodeList_cancer_gene_known,node_color='red',alpha=0.6)
    if nodeList_cancer_gene_candidate:
        nx.draw_networkx_nodes(G,pos,node_size=1200,nodelist=nodeList_cancer_gene_candidate,node_color='yellow',alpha=0.6)
    if nodeList_new:
        nx.draw_networkx_nodes(G,pos,node_size=1200,nodelist=nodeList_new,node_color='w',alpha=0.6)

    nx.draw_networkx_edges(G,pos,alpha=0.4,node_size=0,width=1,edge_color='k')
    nx.draw_networkx_edges(G,pos,alpha=0.3,width=5,edge_color='m')
    nx.draw_networkx_labels(G,pos,fontsize=14)

    plt.axis('off')
    plt.savefig(network_out_file_name, dpi=275)
    plt.show()

def draw_ppi(G,known_cancer_gene_file,candidate_cancer_gene_file,network_out_file_name):
    '''
    '''
    cancer_gene_known_file=open(known_cancer_gene_file,"r")
    cancer_gene_candidate_file=open(candidate_cancer_gene_file,"r")
    cancer_gene_known=[]
    for line in cancer_gene_known_file:
        cancer_gene_known.append(line[0:-1])
    cancer_gene_candidate=[]
    for line in cancer_gene_candidate_file:
        cancer_gene_candidate.append(line[0:-1].upper())

    nodeList_cancer_gene_known=[]
    nodeList_cancer_gene_candidate=[]
    nodeList_new=[]
    for node in G.nodes():
        if node in cancer_gene_known:
            nodeList_cancer_gene_known.append(node)
        elif node in cancer_gene_candidate:
            nodeList_cancer_gene_candidate.append(node)
        else:
            nodeList_new.append(node)

    try:
        pos=nx.graphviz_layout(G)
    except:
        pos=nx.spring_layout(G,iterations=20)

    plt.rcParams['text.usetex'] = False
    if G.number_of_nodes() > 100:
        plt.figure(figsize=(30,30))
    elif G.number_of_nodes() < 100 and G.number_of_nodes() > 50:
        plt.figure(figsize=(20,20))
    else:
        plt.figure(figsize=(10,10))

    if nodeList_cancer_gene_known:
        nx.draw_networkx_nodes(G,pos,node_size=1200,nodelist=nodeList_cancer_gene_known,node_color='red',alpha=0.6)
    if nodeList_cancer_gene_candidate:
        nx.draw_networkx_nodes(G,pos,node_size=1200,nodelist=nodeList_cancer_gene_candidate,node_color='yellow',alpha=0.6)
    if nodeList_new:
        nx.draw_networkx_nodes(G,pos,node_size=1200,nodelist=nodeList_new,node_color='w',alpha=0.6)

    nx.draw_networkx_edges(G,pos,alpha=0.4,node_size=0,width=1,edge_color='k')
    nx.draw_networkx_edges(G,pos,alpha=0.3,width=5,edge_color='m')
    nx.draw_networkx_labels(G,pos,fontsize=14)

    plt.axis('off')
    plt.savefig(network_out_file_name, dpi=275)
    plt.show()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, \
                                    usage='\n\npython draw.py [-h] [-n network file] [-e enhancer-gene pairs] [-m snp file] [-M method] [-l layout][-o output prefix]\n', \
                                    description='', epilog='Ennet: find sub-network with a substantial load of enhancer mutations.')
     parser.add_argument('-i', '--input',  required=True, help='Input graph.', metavar='', dest="input")
     parser.add_argument('-o', "output", required=True, help="Ouput prefix.", metavar='', dest="output")
    parser.add_argument('-l', '--layout',  default="graphviz", help='Choose layout style: graphviz, circular, random, spring, shell', metavar='', dest="layout")

    args = parser.parse_args()

    plot_hotspot_network_known(G, output, '/home/niuyw/Project/RegulatorySNP_170808/ennet/data/anno/netpath.tab', layout = layout)

