# Demo use of the ReactionNetwork class.
# demo_MAPK.csv
from SSApy import ReactionNetwork
import matplotlib.pyplot as plt
import numpy as np
import importlib
import pandas as pd

###########Example 1########################
reaction_path='demo_MAPK.csv'
# construct a class
network = ReactionNetwork.from_csv(reaction_path)#, ker_basis ="rref", cq_basis = "rref"
#%%
#conserved quantities
if len (network.ns2)>0:
    print("CQ", network.cons_list)
else:
    print ("No CQ")

# limitset
limitset_list = ReactionNetwork.compute_bs(network)#, cq_pert=False, cq_oc=False
# %%
# draw a hierarchy
#pygraphviz is required
hier_agraph = ReactionNetwork.make_hiergraph(limitset_list)
hier_agraph.draw(f'./result/hier_demo_MAPK.png')

###determine sign###
smat_sign=network.compute_smat_sign(N=100)
df_sign = pd.DataFrame (smat_sign[:network.M,], index = network.cpd_list_noout, columns = [x[0] for x in network.reac_cons_list])
df_sign.to_csv ("./result/demo_MAPK_sign.csv")
