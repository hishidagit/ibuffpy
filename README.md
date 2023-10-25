# ibuffpy

This package applies structural sensitivity analysis to chemical reaction networks to discover buffering structures within them.


## Installation
To install ibuffpy, run the following command:

```
pip install ibuffpy
```
The pygraphviz package is required to visualize the hierarchy graph of buffering structures.
See https://pygraphviz.github.io/documentation/stable/install.html for the installation of pygraphviz.

## Usage
Once you have installed ibuffpy, you can import ```ReactionNetwork``` module from the ibuffpy package
```
from ibuffpy import ReactionNetwork
```
In ibuffpy, reaction networks are treated as instances of the ```ReactionNetwork.ReactionNetwork``` which contains information about reaction names, metabolites, and stoichiometric matrices.
To read a network from a csv file, use ```ReactionNetwork.from_csv(" PATH_TO_CSV")``` where the csv file must contain reaction formulas such as "reaction_idx, substrate1 substrate2, product1 product2" (an example can be found in demo_network1.csv).
If the reaction equation contains more than one metabolite, separate them with a space.
If the coefficient of a metabolite is $n$, it will appear n times on one side.
Outflow reactions or inflow reactions are represented as `['inflow_name',['metabolite_name'],['out']]` or `['outflow_name',['out'],['metabolite_name']]`.
The metabolite name 'out' is reserved for the node corresponding to outside of the network, which does not appear in the stoichiometric matrix.
```
network = ReactionNetwork.from_csv("PATH_TO_CSV")
network.info()
# list of reactions in the network
print(network.reaction_list)
# list of metabolites in the network except for 'out' node
print(network.cpd_list_nooout)
```
Structural sensitivity analysis computes the network sensitivity to parameter perturbation from the A-matrix, which reflects network structure.
`ReactionNetwork.compute_amat(network)` returns an A-matrix with random real values in nonzero entries of $\partial \boldsymbol r / \partial \boldsymbol x$, a basis of $\ker \nu$, and a basis of $\ker \nu^\top$.
`ReactionNetwork` class uses `scipy.linalg.null_space` to obtain bases by default, but you can also use `ker_basis='rref'` option.
`ker_bais` or `cq_basis` accepts ndarray or list to provide a self-made bases.
```
# numpy is used to calculate a basis by default.
network = ReactionNetwork.from_csv(path_to_csv, ker_basis='svd', cq_basis='rref')

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
