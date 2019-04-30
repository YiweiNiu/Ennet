#!/usr/bin/env python
import sys
import matplotlib as mpl
mpl.use('Agg')
import draw
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

cancer = sys.argv[1].strip()

G = nx.read_edgelist('%s_edges.txt' %(cancer), nodetype=str)

with open('%s_nodes.txt' %(cancer), 'r') as fin:
	for line in fin:
		line = line.strip()
		G.add_node(line)

draw.plot_hotspot_network_known(G, cancer, '/home/niuyw/Project/RegulatorySNP_170808/ennet/data/anno/netpath.tab', layout = "graphviz")

