# Structural-Sensitivity-Analysis

This package applies structure sensitivity analysis to chemical reaction networks to discover buffering structures within them. 
In this package, reaction networks are treated as ReactionNetwork.ReactionNetwork class with information on reactions, metabolites and stoichiometric matrices.
To read a network from a csv file, use ReactionNetwork.from_csv(" PATH_TO_CSV")); the csv file must contain reaction formulas such as "reaction_idx, substrate1 substrate2, product1 product2" (demo_network1.csv is an example). 
If the reaction equation contains more than one metabolite, separate them with a space. 
If the coefficient of a metabolite is greater than 1, it will appear more than once on one side.

You can find buffering structures by limitset_list = ReactionNetwork.compute_limitset(network) and visualise them by graph = ReactionNetwork.make_hiergraph(limitset_list) and graph.draw(PATH_SAVE).
When the A-matrix of the network (network.compute_amat()) is singular, this method cannot be applied to it.

