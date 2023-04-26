"""ReactionNetwork class for SSA.
"""
# %%
import csv
import numpy as np
from scipy import linalg
import pandas as pd
# import func.compute_limitset_meansmat
# import func.compute_rref
# import func.make_hiergraph
# import func.compute_smat
from . import func
import cobra
from sklearn.decomposition import TruncatedSVD
from scipy.sparse import csr_matrix
import sympy
'''
format of reaction_list
= [('reaction0', [substrate0, substrate1], [product0, product1])
    ('reaction1', [substrate0], [product2], [activator0],[inhibitor0])]
'''
np.ndarray
class ReactionNetwork:
    """ReactionNetwork class for SSA."""
    def __init__(self, reaction_list_input, info=True, ker_basis="numpy", cq_basis="numpy"):
        """Initialize ReactionNetwork class.
        
        Parameters
        ----------
        reaction_list_input : list
            List of reactions. Each reaction contains at least 3 elements, that is name, list of substrate names, and list of product names.
            If the reaction is regulated, 4th element is a list of activator names, and 5th element is a list of inhibitor names.
        info : bool, optional
            If True, print information about the network when constructed, by default True
        ker_basis : str or numpy.ndarray, optional
            If "numpy", compute basis of nullspace of stoichiometric matrix using numpy, by default "numpy".
            "sympy" is also available, but it is slower than "numpy".
            It can also be given by yourself as a numpy array.
        cq_basis : str or numpy.ndarray, optional
            If "numpy", compute basis of CQ using numpy, by default "numpy".
            "sympy" is also available, but it is slower than "numpy".
            It can also be given by yourself as a numpy array.
        
        Examples
        --------
        >>> reaction_list = [('reaction0', [substrate0, substrate1], [product0, product1])
                            ('reaction1', [substrate0], [product2], [activator0],[inhibitor0])]
        >>> network = ReactionNetwork(reaction_list, info=True, ker_basis="numpy", cq_basis="numpy")
        """
        if info:
            print('constructed.')

        #drop duplicates
        nodup=[]
        for reac in reaction_list_input:
            if reac not in nodup:
                nodup.append(reac)
        reaction_list_input=nodup

        #regulation
        _reaction_list=[]
        for reac in reaction_list_input:
            if len(reac)==3:
                _reaction_list.append(reac)
            elif len(reac)==4:
                _reaction_list.append([reac[0], reac[1]+reac[3], reac[2]+reac[3]])
            elif len(reac)==5:
                _reaction_list.append([reac[0], reac[1]+reac[3]+reac[4], reac[2]+reac[3]+reac[4]])
            else:
                print('reaction_list format is not correct.')
                1/0

        self.reaction_list=_reaction_list
        self.reaction_list_noid = [[reac[1], reac[2]]
                                    for reac in self.reaction_list]
        self.reactionNames = [reac[0] for reac in self.reaction_list]
        self.reaction_list_reg=reaction_list_input

        self.cpd_list = []
        for rn in self.reaction_list_noid:
            for form in [rn[0], rn[1]]:
                for cpd in form:
                    if not cpd in self.cpd_list:
                        self.cpd_list.append(cpd)
        self.cpd_list = sorted(self.cpd_list)

        self.R = len(self.reaction_list_noid)
        self.M = len(self.cpd_list)
        if 'out' in self.cpd_list:
            self.M -= 1
            self.cpd_list_noout = self.cpd_list.copy()
            self.cpd_list_noout.remove('out')
            if info:
                print('out node')
        else:
            self.cpd_list_noout = self.cpd_list.copy()
            if info:
                print('no out node')
        if info:
            print('M = ', self.M)
            print('R = ', self.R)

        # default tolerance is set to 1e-10
        self.tol = 1e-10

        # stoichiometric matrix
        self.stoi = self.make_stoi()
        # self.ns = func.compute_rref.compute_rref(linalg.null_space(self.stoi).T).T
        # self.ns2 = func.compute_rref.compute_rref(linalg.null_space(self.stoi.T).T)

        #nullspace
        if not isinstance(ker_basis, (np.ndarray, list)):
            if ker_basis == "numpy":
                try :
                    ns=linalg.null_space(self.stoi)

                except np.linalg.LinAlgError:
                    # when network size is large, ns(cycles) are computed using truncated SVD
                    # dimension of nullspace of self.stoi
                    stoich = self.stoi
                    # compute truncated SVD of stoich.T * stoich
                    svd = TruncatedSVD(n_components=stoich.shape[1])
                    svd.fit(stoich.T @ stoich)
                    # vectors of right singular vectors correspoinding to small singular values
                    nullspace_mask = svd.singular_values_ < self.tol
                    nullspace = svd.components_[nullspace_mask].T
                    # check if random estimation is correct
                    if np.linalg.norm(stoich @ nullspace) > self.tol:
                        raise Exception('Error: nullspace estimation is not correct.')
                    else:
                        ns=nullspace

            elif ker_basis == "sympy":
                ns=sympy.Matrix(self.stoi).nullspace()
                if len (ns):
                    ns=np.concatenate([np.array(vec, dtype=float) for vec in ns], axis=1)

            else:
                raise ValueError ('ker_basis must be one of "numpy"", "sympy" or a numpy array/list')

            if len(ns)==0:
                self.ns=np.empty((self.R,0))
            else:
                self.ns=ns
        else:#basis of nullspace can be given by yourself
            self.ns = ker_basis


        #CQ
        if not isinstance(cq_basis, (np.ndarray, list) ):
            if cq_basis == "numpy":
                try:
                    ns2=linalg.null_space(self.stoi.T)

                except np.linalg.LinAlgError:
                    # dimension of nullspace of self.stoi
                    mrx = self.stoi.T
                    # compute truncated SVD of stoich.T * stoich
                    svd = TruncatedSVD(n_components=mrx.shape[1])
                    svd.fit(mrx.T @ mrx)
                    # vectors of right singular vectors correspoinding to small singular values
                    nullspace_mask = svd.singular_values_ < self.tol
                    nullspace = svd.components_[nullspace_mask].T
                    # check if random estimation is correct
                    if np.linalg.norm(mrx @ nullspace) > self.tol:
                        raise Exception('Error: nullspace estimation is not correct.')
                    else:
                        ns2=nullspace

            elif cq_basis == "sympy":
                ns2=sympy.Matrix(self.stoi.T).nullspace()
                if len (ns2):
                    ns2=np.concatenate([np.array(vec, dtype=float) for vec in ns2], axis=1)

            else:
                raise ValueError ('cq_basis must be one of "numpy"", "sympy" or a numpy array/list')



            if len(ns2)==0:
                self.ns2=np.empty((0,self.M))
            else:
                self.ns2=ns2.T
        else:#basis of conserved quantity can be given by yourself
            self.ns2 = cq_basis

        self.A = self.M+len(self.ns.T)
        self.graph = [self.cpd_list_noout, self.reaction_list]
        self.cons_list, self.cons_list_index=self.make_conslist()

        self.reac_cons_list=self.reaction_list+self.cons_list

        # # regularity of the A-matrix
        # if self.A < 100:
        #     self.regularity = np.linalg.matrix_rank(self.compute_amat()) == self.A
        # else:
        #     amat_sp = csr_matrix(self.compute_amat())
        #     self.regularity = np.linalg.matrix_rank(amat_sp) == self.A



    def info(self):
        """print information of the network"""
        print(f'M = {self.M}')
        print(f'R = {self.R}')
        if 'out' in self.cpd_list:
            print('outnode exists')
        else:
            print('no outnode')
        print('cyc = ', len(self.ns.T))
        print('cons = ', len(self.ns2))
        print('rank A = ', np.linalg.matrix_rank(self.compute_amat()))
        print('det A = ', np.linalg.det(self.compute_amat()))
        # print('regularity = ', self.regularity)

    def make_stoi(self):
        """make stoichiometric matrix"""
        stoi = np.zeros((self.M, self.R), dtype=float)
        for r,reac in enumerate(self.reaction_list):
            for m,cpd in enumerate(self.cpd_list_noout):
                cpd_sub = reac[1].count(cpd)  # count in substrate
                cpd_pro = reac[2].count(cpd)  # count in product
                stoi[m, r] = -cpd_sub+cpd_pro
        return stoi

    def make_conslist(self):
        """make list of conserved quantities"""
        cons_list=[]
        cons_list_index=[]
        ns2=self.ns2

        #metabolites are identified by its name
        for c in range(len(ns2)):
            m_list=[self.cpd_list_noout[m] for m in range(self.M) if np.abs(ns2[c,m])>self.tol]
            cons=['cons_'+str(c), m_list]
            cons_list.append(cons)

        #index of metabolites
        for c in range(len(ns2)):
            m_list=[m for m in range(self.M) if np.abs(ns2[c,m])>self.tol]
            cons=['cons_'+str(c), m_list]
            cons_list_index.append(cons)

        return cons_list, cons_list_index

    def compute_amat(self):
        """compute the A-matrix. 
        Random real value is assigned to each nonzero element in the left upper part of the A-matrxis.
        Cycles and conserved quantities are computed using numpy or sympy, or explicitly given.
        """
        R = self.R
        M = self.M
        A = self.A
        ns2 = self.ns2
        ns = self.ns
        cpd_list_noout = self.cpd_list_noout
        reaction_list_noid = self.reaction_list_noid

        amat = np.zeros((A, A), dtype=float)
        amat[R:A, :M] = -ns2
        amat[:R, M:A] = -ns

        # create rmat
        rmat = np.zeros((R, M))
        for r in range(R):
            for m in range(M):
                if cpd_list_noout[m] in reaction_list_noid[r][0]:  # subに含まれる
                    rmat[r, m] = np.random.rand()+0.1
        # create amat
        amat[:R, :M] = rmat

        return amat

    def compute_smat(self):
        """compute the sensitivity matrix"""
        return func.compute_smat.compute_smat(self)

    def compute_smat_mean(self, N=10, large_error=False):
        """compute the sensitivity matrix of mean"""
        return func.compute_smat.compute_smat_mean(self, N, large_error=large_error)

    def check_ocomp(self, subg):
        """check if output complete
        Args:
            subg (list): list of metabolites and reactions"""
        reaction_list = self.reaction_list
        R = self.R

        sub_m_list = subg[0]
        sub_r_list = subg[1]

        # check if output complete
        for cpd in sub_m_list:
            for rxn in reaction_list:
                if (cpd in rxn[1]) and (rxn[0] not in sub_r_list):
                    return False
        return True

    def compute_cyc(self, sub_r_index):
        """compute the number of cycles from the list of reaction index.
        
        Parameters
        ----------
        sub_r_index : list
            index list of reactions

        Returns
        -------
        num_cyc : int
            number of cycles composed of the reactions in sub_r_index
        """
        if len(sub_r_index) == 0:
            num_cyc = 0
        else:
            sub_rstoi = self.stoi[:, sub_r_index]
            num_cyc = len(linalg.null_space(sub_rstoi).T)
        return num_cyc

    def compute_cons(self, sub_m_index):
        """compute the number of conserved quantities from the list of metabolite index.
        Parameters
        ----------
        sub_m_index : list
            index list of metabolites
        
        Returns
        -------
        num_cons : int
            number of conserved quantities composed of the metabolites in sub_m_index
        """

        if len(self.ns2)==0:#no cons exists
            num_cons = 0

        else:
            nsub_m_index = [m for m in range (self.M)
                            if self.cpd_list_noout[m] not in sub_m_index]
            if not nsub_m_index:  # contains all metabolites
                num_cons = len(self.ns2)
            else:
                nsubstoi = self.stoi[nsub_m_index, :]
                # num_cons = len(self.ns2)-len(linalg.null_space(nsubstoi.T).T)
                # if network is large, null_space calculation does not converge
                num_cons=len(self.ns2)-(len(nsub_m_index)-np.linalg.matrix_rank(nsubstoi))
        return num_cons

    def index_subg(self, subg):
        """compute the index of subgraph
        
        Parameters
        ----------
        subg : list of list
            list of metabolites and reactions. reactions are identified by their names.
        
        Returns
        -------
        index : int
            index of subgraph
        """
        stoi = self.stoi
        cpd_list_noout = self.cpd_list_noout
        sub_m_list = subg[0]

        sub_r_index = []
        for n, reac in enumerate(self.reaction_list):
            for rname in subg[1]:
                if reac[0] == rname:
                    sub_r_index.append(n)
                    break

        # number of cycles in the subgraph
        num_cyc =self.compute_cyc(sub_r_index)

        # number of conserved quantities in the subgraph
        num_cons=self.compute_cons (sub_m_list)

        index = len(sub_m_list)+num_cyc-len(sub_r_index)-num_cons

        return index

    def short_name(self, name):
        """Insert one line break for every N characters.

        Parameters
        ----------
        name : str
            name of a box in hierarchy graph of buffering structures.
        
        Returns
        -------
        name : str"""
        l = len(name)
        N = 40  # One line break for every N characters.
        for i in range(l//N):
            name = name[:(i+1)*N]+'\n'+name[(i+1)*N:]
        return name

    def make_ocompSubg(self, subm_list):
        """make output complete subgraph from the list of metabolites.
        
        Parameters
        ----------
        subm_list : list
            list of metabolite names
        
        Returns
        -------
        subg : list of list
            list of metabolites and reactions. reactions are identified by their names.
        """
        reaction_list = self.reaction_list
        subr_list = []
        for reac in reaction_list:
            if not set(reac[1]).isdisjoint(set(subm_list)):
                subr_list.append(reac[0])
        return [subm_list, subr_list]

    def to_df(self):
        """Return pandas dataframe of reactions in the network.
        
        Returns
        -------
        df_reaction : pandas dataframe
        """
        maxlen_reaction=max([len(reac) for reac in self.reaction_list])
        if maxlen_reaction==3:
            df_reaction = pd.DataFrame(self.reaction_list,
             columns=['name','substrates', 'products'])
            df_reaction=df_reaction.set_index('name')
            return df_reaction

        else:
            print(maxlen_reaction,'inhibitor/activator is not supported in this function')

    def get_reacCons_by_id(self,_id):
        for rc in self.reac_cons_list:
            if rc[0]==_id:
                return rc

def compute_limitset(network,N=10,large_error=True,detectCQ=True):
    """compute the limit set of the network
    Parameters:
    ----------
    network : ReactionNetwork
        reaction network
    N : int, optional
        number of iterations, by default 10
    large_error : bool, optional
        if True, return error if computational error is large while S-matrix computation, by default True
    detectCQ : bool, optional
        if True, include conserved quantities in buffeing structures, by default True
    
    Returns
    -------
    limitset_list : list of list
        list of buffering structures
    """
    return func.compute_limitset_meansmat.compute_limitset_meansmat(network,N,large_error=large_error,detectCQ=detectCQ)

def make_hieredge(limitset_list):
    return func.make_hiergraph.make_hieredge(limitset_list)

def make_hiergraph(limitset_list):
    return func.make_hiergraph.make_hiergraph(limitset_list)
# %%
def from_csv(path,info=True,ker_basis="numpy",cq_basis="numpy"):
    """ load reaction network from csv file
    
    Parameters
    ----------
    path : str
        path to csv file
    info : bool, optional
        if True, print information of the network, by default True
    ker_basis : str, optional
        "numpy" or "sympy", by default "numpy"
    cq_basis : str, optional
        "numpy" or "sympy", by default "numpy"
    """
    reaction_list=[]
    with open(path, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            reaction=[]
            reaction.append(line[0]) #name
            for cpds in line[1:]:#substrate, product, activator, inhibitor
                _cpds=[]
                [_cpds.append(cpd) for cpd in cpds.split(' ') if cpd!='']
                reaction.append(_cpds)
            reaction_list.append(reaction)
    return ReactionNetwork(reaction_list,info=info,ker_basis=ker_basis,cq_basis=cq_basis)

def from_cobra(model,info=True):
    """convert cobra model to ReactionNetwork object
    
    Parameters
    ----------
    model : cobra model
        cobra model
    info : bool, optional
        if True, print information of the network, by default True
    
    Returns
    -------
    ReactionNetwork
    """
    # convert to SSA network
    reaction_list=[]
    for reac in model.reactions:
        # objective function
        # if 'BIOMASS' in reac.id:
        #     continue
        lhs=[]
        rhs=[]
        for reactant in reac.reactants:
            coef=abs(int(reac.metabolites[reactant]))
            lhs+=[reactant.id]*coef
        for product in reac.products:
            coef=abs(int(reac.metabolites[product]))
            rhs+=[product.id]*coef
        if len(lhs)==0:
            lhs=['out']
        if len(rhs)==0:
            rhs=['out']
        if not reac.reversibility:
            reaction_list.append([reac.id,lhs,rhs])
        elif lhs==['out'] or rhs==['out']: # reversible outflow
            reaction_list.append([reac.id,lhs,rhs])
            reaction_list.append([reac.id+'_rev',rhs,lhs])

        # reversible reaction has products as its regulator
        else:
            reaction_list.append([reac.id,lhs,rhs,rhs])

    network=ReactionNetwork(reaction_list,info=info)
    return network

def to_cobra(network,name=''):
    """convert to cobra model
    
    Parameters
    ----------
    network : ReactionNetwork
        ReactionNetwork object
    name : str, optional
        name of the model, by default ''
        
    Returns
    -------
    model : cobra.Model"""
    # convert from SSA network to cobra model
    model_name=name
    model=cobra.Model(model_name)

    for cpdname in network.cpd_list_noout:
        metab=cobra.Metabolite(cpdname)
        model.add_metabolites(metab)
    for reac in network.reaction_list_reg:
        #reaction name of cobra cannot contain white-space(" ").
        reacname=reac[0].replace(' ','_')
        cobra_reac=cobra.Reaction(reacname)
        metab_dict=dict()

        for cpd in reac[1]:#lhs
            if cpd!='out':
                metab_dict[model.metabolites.get_by_id(cpd)]=-reac[1].count(cpd)
            else:
                0
        for cpd in reac[2]:#rhs
            if cpd!='out':
                metab_dict[model.metabolites.get_by_id(cpd)]=reac[2].count(cpd)
            else:
                0
        cobra_reac.add_metabolites(metab_dict)

        #if reaction is reversible, reac[2] is in reac[1]
        if len(reac)>3 and reac[3]==reac[2]:
            if reac[2]!=reac[3]:
                print('reaction_list format is not for cobra model')
                raise(Exception)
            else:
                print(reac[0])
                cobra_reac.lower_bound=-1000

        model.add_reaction(cobra_reac)

    return model
#%%
class LargeErrorSmat(Exception):
    pass