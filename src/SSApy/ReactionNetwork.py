"""ReactionNetwork class for SSA.
"""
# %%
import csv
import numpy as np
from scipy import linalg
import pandas as pd
from .ftn import ftn_compute_bs_meansmat
#from .ftn import ftn_compute_rref
from .ftn import ftn_compute_nullspace
from .ftn import ftn_make_hiergraph
from .ftn import ftn_compute_smat
# from . import func
#import sympy
'''
format of reaction_list
= [('reaction0', [substrate0, substrate1], [product0, product1]),
    ('reaction1', [substrate0], [product2], [activator0],[inhibitor0])]
'''

class ReactionNetwork:
    """ReactionNetwork class for SSA.
        Parameters
        ----------
        reaction_list_input : list
            List of reactions. Each reaction contains at least 3 elements, that is name, list of substrate names, and list of product names.
            If the reaction is regulated, 4th element is a list of activator names, and 5th element is a list of inhibitor names.
        info : bool, optional
            If True, print information about the network when constructed, by default False
        ker_basis : str or numpy.ndarray, default = "svd"
            If "svd", compute the basis of nullspace of stoichiometric matrix using singular value decomposition (with "scipy.linalg.null_space").
            If "rref", compute the basis by putting a matrix into RREF.
            It can also be given by yourself as a numpy array or list. The Size should be M*K (M: The number of reactions, K: The number of cycles).
            "svd" is faster than "rref" and the calculation of the sensitivity is not affected by the choice of ker basis. Therefore we recommend to use "svd".


        cq_basis : str or numpy.ndarray, default = "rref"
            If "rref", compute the basis by putting a matrix into RREF.
            If "svd", compute the basis of nullspace of the stoichiometric matrix using the SVD (singular value decomposition) The scipy function "scipy.linalg.null_space" is used.
            It can also be given by yourself as a numpy array or list. The size should be N*L  ( N: The number of chemicals, L: The number of conserved quantities).
            When calculationg the basis of conserved quantities, "rref" is recommended. Using the SVD leads to an inappropriate choice of the basis (See Yamauchi et al, Apeendix A).

        Examples
        --------
        >>> reaction_list = [('reaction0', [substrate0, substrate1], [product0, product1])
                            ('reaction1', [substrate0], [product2], [activator0],[inhibitor0])]
        >>> network = ReactionNetwork(reaction_list, info=True, ker_basis="rref", cq_basis="rref")
    """

    def __init__(self, reaction_list_input, info=False, ker_basis="svd", cq_basis="rref"):
        """Initialize ReactionNetwork class.
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
        positive_regulation = {}#dict for Negative regulation
        negative_regulation = {}#dict for Negative regulation
        for reac in reaction_list_input:
            if len(reac)==3:
                _reaction_list.append(reac)
            elif len(reac)==4:
                _reaction_list.append([reac[0], reac[1]+reac[3], reac[2]+reac[3]])
                positive_regulation[reac[0]] = reac[3]#record positive regulators
            elif len(reac)==5:
                _reaction_list.append([reac[0], reac[1]+reac[3]+reac[4], reac[2]+reac[3]+reac[4]])
                negative_regulation[reac[0]] = reac[4]#record negative regulators
            else:
                print('reaction_list format is not correct.')
                1/0
        self.positive_regulation = positive_regulation
        self.negative_regulation = negative_regulation

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

        #calculate nullspace (cycle)
        if not isinstance(ker_basis, (np.ndarray, list)):
            if ker_basis == "svd":
                ns = ftn_compute_nullspace.cal_nullspace_svd (self.stoi, error=1.0e-10)

            elif ker_basis == "rref":
                ns = ftn_compute_nullspace.cal_nullspace_rref (self.stoi, error=1.0e-10)

            else:
                raise ValueError ('ker_basis must be one of "svd", "rref", or a numpy array/list')

            if len(ns)==0:
                self.ns=np.empty((self.R,0))
            else:
                self.ns=ns
        else:#basis of nullspace can be given by yourself
            self.ns = ker_basis##ker_basis must be a column vector

        #calculate coker vector (conserved quantity)
        if not isinstance(cq_basis, (np.ndarray, list) ):
            if cq_basis == "rref":
                ns2 = ftn_compute_nullspace.cal_nullspace_rref (self.stoi.T, error=1.0e-10)

            elif cq_basis == "svd":
                ns2 = ftn_compute_nullspace.cal_nullspace_svd (self.stoi.T, error=1.0e-10)

            else:
                raise ValueError ('cq_basis must be one of "rref", "svd", or a numpy array/list')

            if len(ns2)==0:
                self.ns2=np.empty((0,self.M))
            else:
                self.ns2=ns2.T#Convert ns2 into a row vector
        else:#basis of conserved quantity can be given by yourself
            self.ns2 = np.array (cq_basis).T#cq_basis sholud be a column vector→ns2 is a row vector


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
        print(f'A = {self.A}')
        if 'out' in self.cpd_list:
            print('outnode exists')
        else:
            print('no outnode')
        print('cyc = ', len(self.ns.T))
        print('cons = ', len(self.ns2))
        print('rank A = ', np.linalg.matrix_rank(self.compute_amat()))
        # print('det A = ', np.linalg.det(self.compute_amat()))
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
        Cycles and conserved quantities are computed using SVD or RREF, or explicitly given.
        """
        R = self.R
        M = self.M
        A = self.A
        ns2 = self.ns2
        ns = self.ns
        #random_min, random_max = self.random_min, self.random_max
        random_min, random_max = 0.1, 10
        cpd_list_noout = self.cpd_list_noout
        reaction_list_noid = self.reaction_list_noid
        reaction_list_reg =self.reaction_list_reg
        negative_regulation=self.negative_regulation

        amat = np.zeros((A, A), dtype=float)
        amat[R:A, :M] = -ns2
        amat[:R, M:A] = -ns

        # create rmat
        #rmat = np.zeros((R, M))
        #for r in range(R):
            #for m in range(M):
                #if cpd_list_noout[m] in reaction_list_noid[r][0]:  # subに含まれる
                    #rmat[r, m] = np.random.rand()+0.1
        # create rmat (New by Yamauchi) (distinguish positive regulation and negative regulation)
        rmat = np.zeros((R, M))
        for i,r in enumerate(reaction_list_reg):
            for j,m in enumerate(cpd_list_noout):
                if m in reaction_list_noid[i][0]:# subに含まれる
                    rmat[i, j] = np.random.rand()*random_max+random_min
                    if neg_reg:=negative_regulation.get(r[0]):
                        if m in neg_reg:#When m is a negaive regulator of the reaction
                            rmat[i, j] = -np.random.rand()*random_max-random_min
        # create amat
        amat[:R, :M] = rmat

        return amat

    def compute_smat(self):
        """compute the sensitivity matrix"""
        return ftn_compute_smat.compute_smat(self)

    def compute_smat_mean(self, N=10, large_error=False):
        """compute the sensitivity matrix of mean"""
        return ftn_compute_smat.compute_smat_mean(self, N, large_error=large_error)

    def compute_smat_sign(self, N=100):#add by Yamauchi
        """compute the signs of the sensitivity"""
        return ftn_compute_smat.compute_smat_sign(self, N)


    def check_ocomp(self, subg):
        """
        subg (list of list):[[cpd name list], [reaction name list]]
        """
        reaction_list = self.reaction_list_reg#8/22山内修正 (self.reaction_list)
        R = self.R

        sub_m_list = subg[0]
        sub_r_list = subg[1]

        # check if output complete
        for cpd in sub_m_list:
            for r in range(R):
                if (cpd in reaction_list[r][1]) and (reaction_list [r][0] not in sub_r_list):#8/22山内修正
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

    def to_csv(self,path):
        """save network to csv file

        parameters
        ----------
        network : ReactionNetwork
            reaction network
        path : str
            path to csv file
        """
        with open(path, 'w',newline='') as f:
            writer = csv.writer(f)
            for rxn in self.reaction_list:
                writer.writerow([rxn[0],' '.join(rxn[1]),' '.join(rxn[2])])

def compute_bs(network,N=10,large_error=True,detectCQ=True):
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
    bs_list : list of list
        list of buffering structures
    """
    return ftn_compute_bs_meansmat.compute_bs_meansmat(network,N,large_error=large_error,detectCQ=detectCQ)

def make_hieredge(bs_list):
    return ftn_make_hiergraph.make_hieredge(bs_list)

def make_hiergraph(bs_list):
    return ftn_make_hiergraph.make_hiergraph(bs_list)
# %%
def from_csv(path,info=True,ker_basis="svd",cq_basis="rref"):
    """ load reaction network from csv file

    Parameters
    ----------
    path : str
        path to csv file
    info : bool, optional
        if True, print information of the network, by default True
    ker_basis : str, optional
        "svd" or "rref", by default "svd"
    cq_basis : str, optional
        "svd" or "rref", by default "rref"
    """
    reaction_list=[]
    with open(path, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            if len(line)==0:
                continue
            reaction=[]
            reaction.append(line[0]) #name
            for cpds in line[1:]:#substrate, product, activator, inhibitor
                _cpds=[]
                [_cpds.append(cpd) for cpd in cpds.split(' ') if cpd!='']
                reaction.append(_cpds)
            reaction_list.append(reaction)
    return ReactionNetwork(reaction_list,info=info,ker_basis=ker_basis,cq_basis=cq_basis)


def from_pandas (df, sep = " ",info=True,ker_basis="svd",cq_basis="rref"):
    df=df.copy ()

    if len (df.columns)==3:
        df.columns = ["Reaction_index", "Substrate", "Product"]
        df["Activator"] = np.nan
        df["Inhibitor"] = np.nan
    elif len (df.columns)==4:
        df.columns = ["Reaction_index", "Substrate", "Product", "Activator"]
        df["Inhibitor"] = np.nan
    elif len (df.columns)==5:
        df.columns = ["Reaction_index", "Substrate", "Product", "Activator", "Inhibitor"]
    else:
        raise ValueError ("The number of columns in the dataframe must be between 3 and 5")


    if not df["Reaction_index"].dtypes =='str':
        df["Reaction_index"]= df["Reaction_index"].astype(str)

    df.loc[df["Substrate"].isnull(), "Substrate"]="out"
    df.loc[df["Product"].isnull(), "Product"]="out"

    df["Substrate"] = df["Substrate"].apply (lambda x:x.split (sep))
    df["Product"] = df["Product"].apply (lambda x:x.split (sep))
    df["Activator"] = df["Activator"].apply (lambda x:x.split (sep) if not pd.isnull(x) else "nan")
    df["Inhibitor"] = df["Inhibitor"].apply(lambda x:x.split (sep) if not pd.isnull(x) else "nan")

    list_df=  df.values.tolist()

    input_crn=[]
    for reac in list_df:
        if reac[4]!= 'nan':
            input_crn.append ([k if k !="nan" else [] for k in reac])
        else:
            input_crn.append ([x for x in reac if x != 'nan'] )
    return ReactionNetwork(input_crn,info=info,ker_basis=ker_basis,cq_basis=cq_basis)


def add_reactions(network,newrxns,info=True):
    """add reactions to network

    Parameters
    ----------
    network : ReactionNetwork
        reaction network
    newrxns : list of reactions
        reactions to be added
    info : bool, optional
        if True, print information of the network, by default True
    """
    if newrxns==[]:
        return network
    for rxn in newrxns:
        if rxn[0] in network.reactionNames:
            print('reaction name already exists')
            raise(Exception)
    reaction_list = network.reaction_list+newrxns
    return ReactionNetwork(reaction_list,info=info)

#%%
