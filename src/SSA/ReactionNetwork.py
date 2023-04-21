"""Applying Structural Sensitivity Analysis to a reaction network.
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

class ReactionNetwork:
    def __init__(self, reaction_list_input, info=True, ker_basis="numpy", cq_basis="numpy"):
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

        negative_regulation = {}#dict for Negative regulation
        for reac in reaction_list_input:
            if len(reac)==3:
                _reaction_list.append(reac)
            elif len(reac)==4:
                _reaction_list.append([reac[0], reac[1]+reac[3], reac[2]+reac[3]])
            elif len(reac)==5:
                _reaction_list.append([reac[0], reac[1]+reac[3]+reac[4], reac[2]+reac[3]+reac[4]])
                negative_regulation[reac[0]] = reac[4]#record negative regulators of each reaction
            else:
                print('reaction_list format is not correct.')
                1/0


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
        # print information
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
        #return the stoichiometric matrix
        stoi = np.zeros((self.M, self.R), dtype=float)
        for r,reac in enumerate(self.reaction_list):
            for m,cpd in enumerate(self.cpd_list_noout):
                cpd_sub = reac[1].count(cpd)  # count in substrate
                cpd_pro = reac[2].count(cpd)  # count in product
                stoi[m, r] = -cpd_sub+cpd_pro
        return stoi

    def make_conslist(self):
        #return the list of conserved quantities
        cons_list=[]
        cons_list_index=[]
        ns2=self.ns2

        #metabolites are identified by its name
        for c in range(len(ns2)):
            m_list=[self.cpd_list_noout[m] for m in range(self.M) if np.abs(ns2[c,m])>1.0e-10]
            cons=['cons_'+str(c), m_list]
            cons_list.append(cons)

        #index of metabolites
        for c in range(len(ns2)):
            m_list=[m for m in range(self.M) if np.abs(ns2[c,m])>1.0e-10]
            cons=['cons_'+str(c), m_list]
            cons_list_index.append(cons)

        return cons_list, cons_list_index

    def compute_amat(self):
        R = self.R
        M = self.M
        A = self.A
        ns2 = self.ns2
        ns = self.ns
        cpd_list_noout = self.cpd_list_noout
        reaction_list_noid = self.reaction_list_noid
        reaction_list_reg =self.reaction_list_reg
        negative_regulation=self.negative_regulation

        amat = np.zeros((A, A), dtype=float)
        amat[R:A, :M] = -ns2
        amat[:R, M:A] = -ns

        # create rmat
        rmat = np.zeros((R, M))
        for i,r in enumerate(reaction_list_reg):
            for j,m in enumerate(cpd_list_noout):
                if m in reaction_list_noid[i][0]:# subに含まれる
                    rmat[i, j] = np.random.rand()+0.1
                    if neg_reg:=negative_regulation.get(r[0]):
                        if m in neg_reg:#When m is a negaive regulator of the reaction
                            rmat[i, j] = -np.random.rand()-0.1



        # create amat
        amat[:R, :M] = rmat

        return amat

    def compute_smat(self):
        # compute the sensitivity matrix
        return func.compute_smat.compute_smat(self)

    def compute_smat_mean(self, N=10, large_error=False):
        # compute the sensitivity matrix of mean
        return func.compute_smat.compute_smat_mean(self, N, large_error=large_error)

    def check_ocomp(self, subg):
        reaction_list = self.reaction_list
        R = self.R

        sub_m_list = subg[0]
        sub_r_list = subg[1]

        # check if output complete
        for cpd in sub_m_list:
            for r in range(R):
                if (cpd in reaction_list[r][1]) and (r not in sub_r_list):
                    return False
        return True

    def compute_cyc(self, sub_r_list):
        # return the number of cycles from the list of reactions
        if len(sub_r_list) == 0:
            num_cyc = 0
        else:
            sub_rstoi = self.stoi[:, sub_r_list]
            num_cyc = len(linalg.null_space(sub_rstoi).T)
        return num_cyc

    def compute_cons(self, sub_m_list):
        # return the number of conserved quantities from the list of metabolites

        if len(self.ns2)==0:#no cons exists
            num_cons = 0

        else:
            nsub_m_index = [m for m in range (self.M)
                            if self.cpd_list_noout[m] not in sub_m_list]
            if not nsub_m_index:  # contains all metabolites
                num_cons = len(self.ns2)
            else:
                nsubstoi = self.stoi[nsub_m_index, :]
                # num_cons = len(self.ns2)-len(linalg.null_space(nsubstoi.T).T)
                # if network is large, null_space calculation does not converge
                num_cons=len(self.ns2)-(len(nsub_m_index)-np.linalg.matrix_rank(nsubstoi))
        return num_cons

    def index_subg(self, subg):
        # return the index of subgraph
        stoi = self.stoi
        cpd_list_noout = self.cpd_list_noout
        sub_m_list = subg[0]

        sub_r_index = []
        for n, reac in enumerate(self.reaction_list):
            for rname in subg[1]:
                if reac[0] == rname:
                    sub_r_index.append(n)
                    break

        ns2 = self.ns2

        # number of cycles in the subgraph
        num_cyc =self.compute_cyc(sub_r_index)

        # number of conserved quantities in the subgraph
        num_cons=self.compute_cons (sub_m_list)

        index = len(sub_m_list)+num_cyc-len(sub_r_index)-num_cons

        return index

    def short_name(self, name):
        l = len(name)
        N = 40  # N文字に一回改行
        for i in range(l//N):
            name = name[:(i+1)*N]+'\n'+name[(i+1)*N:]
        return name

    def make_ocompSubg(self, subm_list):
        reaction_list = self.reaction_list
        subr_list = []
        for reac in reaction_list:
            if not set(reac[1]).isdisjoint(set(subm_list)):
                subr_list.append(reac[0])
        return [subm_list, subr_list]

    def to_df(self):
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

def compute_limitset(network,N=10,large_error=True):
    return func.compute_limitset_meansmat.compute_limitset_meansmat(network,N,large_error=large_error)

def make_hieredge(limitset_list):
    return func.make_hiergraph.make_hieredge(limitset_list)

def make_hiergraph(limitset_list):
    return func.make_hiergraph.make_hiergraph(limitset_list)
# %%
def from_csv(path,info=True,ker_basis="numpy",cq_basis="numpy"):
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

#############################Convert panda.dataframe to the input of ReactionNetwork########################################
#########Parameters
#########df:pandas.dataframe (Input dataframe (Each row corresponds to one reaction))#########
#0th column = reaction id, 1th column = substrates of a reaction, 2nd column = products of a reaction, 3th column = regulator of a reaction)
#In principle, all regulators are considered to be positive regulators of a reaction.
#If you want to distinguish negative regulator of a reaction from positive regulators, 3 th columns = activator, 4th column = inhibitor
#########sep:{" ", ","} #########

#########ker_basis: {"numpy", "sympy", 2D numpy array}#########
# How to calculate the basis of nullspace
#########cq_basis:{"numpy", "sympy", 2D numpy array}#########
#How to calculate the basis of conserved quantities
def from_pandas (df, sep = " ",info=True,ker_basis="numpy",cq_basis="numpy"):
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

def from_cobra(model,info=True):
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
