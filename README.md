# SSApy

This package applies structure sensitivity analysis to chemical reaction networks to discover buffering structures within them. 


## Installation
To install SSApy, run the next command:

  pip install (path to the repository setup.py is located)

"-e" option enables you to edit the package by yourself.

## Usage
When you have installed SSApy, you can import ReactionNetwork module from the SSApy package;
```
from SSApy import ReactionNetwork
```
In SSApy, reaction networks are treated as ReactionNetwork.ReactionNetwork class with information of reaction name, metabolites, and stoichiometric matrices.
To read a network from a csv file, use ReactionNetwork.from_csv(" PATH_TO_CSV") where the csv file must contain reaction formulas such as "reaction_idx, substrate1 substrate2, product1 product2" (demo_network1.csv is an example). 
If the reaction equation contains more than one metabolite, separate them with a space. 
If the coefficient of a metabolite is n, it will appear n times on one side.
Outflow reaction or inflow reaction is represented as `['outflow_name',['out'],['metabolite_name]]` or `['outflow_name',['out'],['metabolite_name]]`.
Metabolite name 'out' is reserved for the node corresponding to outside of the network which does not appear in the stoichiometric matrix $\nu$.
```
network = ReactionNetwork.from_csv("PATH_TO_CSV")
network.info()
# list of reactions in the network
print(network.reaction_list)
# list of metabolites in the network except for 'out' node
print(network.cpd_list_nooout)
```
Structural sensitivity analysis computes the network sensitivity to parameter perturbation from A-matrix which reflects network structure.
`ReactionNetwork.compute_amat(network)` returns an A-matrix with random real values in nonzero entries of $\partial \bm r / \partial \bm x$, a basis of $\ker \nu$, and a basis of $\ker \nu^\top$.
ReactionNetwork class uses `scipy.linalg.null_space` to obtain bases by default, but you can also use `nullspace` function of `sympy` instead by `ker_basis='sympy'` option.
`ker_bais` or `coker_basis` accepts ndarray or list to give a self-made bases.
```
# numpy is used to calculate a basis by default.
network = ReactionNetwork.from_csv(path_to_csv, ker_basis='numpy', coker_basis='numpy')

# augumented matrix
amat = network.compute_amat()

# sensitivity matrix
smat = np.linalg.inv(amat)
```

You can find buffering structures by 
```
bs_list = ReactionNetwork.compute_bs(network)
``` 
and visualise them by 
```
graph = ReactionNetwork.make_hiergraph(bs_list) 
graph.draw(PATH_SAVE)
```