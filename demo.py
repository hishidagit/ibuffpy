# %%
# Demo use of the ReactionNetwork class.
# python 3.8.10
# to draw hierarchy graph, pygraphviz is required.
from networkx.drawing.nx_agraph import graphviz_layout
import networkx as nx
from SSA import ReactionNetwork
import matplotlib.pyplot as plt
import numpy as np
import importlib
import pandas as pd
importlib.reload(ReactionNetwork)
# %%
# network is defied by a list of reactions
# read csv of reaction data

networkName = 'demo_network1'
networkName = 'demo_MAPK'

# %%
#path to csv
reaction_path=f'./{networkName}.csv'
# construct a class
network = ReactionNetwork.from_csv(reaction_path)
#%%
#conserved quantities
print(network.cons_list)
print(network.cons_list_index)
#%%
# limitset
limitset_list = ReactionNetwork.compute_limitset(network)
# %%
# draw a hierarchy
#pygraphviz is required
hier_agraph = ReactionNetwork.make_hiergraph(limitset_list)
hier_agraph.draw(f'./hier_{networkName}.png')
#%%
# ignoring conserved quantities
limitset_list_noCQ = ReactionNetwork.compute_limitset(network, detectCQ=False)
hier_agraph_noCQ = ReactionNetwork.make_hiergraph(limitset_list_noCQ)
hier_agraph_noCQ.draw(f'./hier_{networkName}_noCQ.png')
# %%
# A matrix
fig = plt.figure(figsize=(10, 10), facecolor='w')
amat=network.compute_amat()
amat_binary=np.where(np.abs(amat)>0, 1, 0)
plt.imshow(amat_binary)
plt.yticks(np.arange(network.R), network.reactionNames, size=20)
fig.autofmt_xdate(rotation=45)
plt.xticks(np.arange(network.M), network.cpd_list_noout, size=20)
plt.title(f'A matirx of {networkName}')
plt.colorbar()
plt.show()
# %%
# S matrix
fig = plt.figure(figsize=(10, 10), facecolor='w')
smat=network.compute_smat()
smat_binary=np.where(np.abs(smat)>1.0e-10, 1, 0)
plt.imshow(smat_binary)
plt.xticks(np.arange(network.R), network.reactionNames, size=20)
plt.yticks(np.arange(network.M), network.cpd_list_noout, size=20)
plt.title(f'S matirx of {networkName}')
plt.colorbar()
plt.show()
# %%
rows=[]
pd.DataFrame(smat_binary)
# %%
df_smat = pd.DataFrame(data=network.compute_smat(), columns=['R{}'.format(i) for i in range(1, len (network.reaction_list)+1)]+['Coker{}'.format(i) for i in range(1, len (network.ns2)+1)], \
                       index=network.cpd_list_noout+['Ker{}'.format(i) for i in range(1, len (network.ns.T)+1)], dtype='float')
df_smat
# %%
