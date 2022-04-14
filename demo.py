# %%
# Demo use of the ReactionNetwork class.
# python 3.8.10
# to draw hierarchy graph, pygraphviz is required.
from networkx.drawing.nx_agraph import graphviz_layout
import networkx as nx
import ReactionNetwork
import matplotlib.pyplot as plt
import numpy as np
import csv
import importlib
importlib.reload(ReactionNetwork)
# %%
# network is defied by a list of reactions
# read csv of reaction data

networkName = 'demo_MAPK'
# networkName = 'demo_network1'

with open(f'./{networkName}.csv', 'r') as f:
    reader = csv.reader(f)
    reaction_list = [[reac[0], reac[1].split(' '), reac[2].split(' ')] for reac in reader]
# %%
# construct a class
network = ReactionNetwork.ReactionNetwork(reaction_list)
#%%
#conserved quantities
print(network.cons_list)
print(network.cons_list_index)
# %%
# limitset
limitset_list = network.find_limitset_meansmat(N=100)
#%%
hoge=ReactionNetwork.compute_limitset(network)
print([a[0] for a in hoge])
#%%
network.make_hieredge(hoge)
# %%
# draw a hierarchy
#pygraphviz is required
hieredge_list = network.make_hieredge(limitset_list)

hier = []
pick = []
for edge in hieredge_list:
    s = ' '.join(edge[0])
    t = ' '.join(edge[1])
    hier.append([network.short_name(s), network.short_name(t)])
    pick.extend(list(edge))


hier_graph = nx.DiGraph(hier)
hier_agraph = nx.nx_agraph.to_agraph(hier_graph)
hier_agraph.layout(prog='dot')
hier_agraph.draw(f'./hier_{networkName}.png')
# %%
# A matrix
fig = plt.figure(figsize=(10, 10), facecolor='w')
plt.imshow(network.compute_amat())
plt.yticks(np.arange(network.R), network.reactionNames, size=20)
fig.autofmt_xdate(rotation=45)
plt.xticks(np.arange(network.M), network.cpd_list_noout, size=20)
plt.title(f'A matirx of {networkName}')
plt.colorbar()
plt.show()
# %%
# S matrix
fig = plt.figure(figsize=(10, 10), facecolor='w')
plt.imshow(network.compute_smat())
plt.xticks(np.arange(network.R), network.reactionNames, size=20)
plt.yticks(np.arange(network.M), network.cpd_list_noout, size=20)
plt.title(f'S matirx of {networkName}')
plt.colorbar()
plt.show()
# %%
