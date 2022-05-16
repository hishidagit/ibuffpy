# %%
import csv
import numpy as np
from scipy import linalg
import networkx as nx
import func.compute_limitset_meansmat
import func.compute_rref
import func.make_hiergraph
import func.compute_smat
'''
format of reaction_list
= [('reaction0', [substrate0, substrate1], [product0, product1])
    ('reaction1', [substrate0], [product2], [regulator0, regulator1])]
'''

class ReactionNetwork:
    def __init__(self, reaction_list_input, info=True):
        if info:
            print('constructed.')

        #regulation 
        _reaction_list=[]
        for reac in reaction_list_input:
            if len(reac)==3:
                _reaction_list.append(reac)
            elif len(reac)==4:
                _reaction_list.append([reac[0], reac[1]+reac[3], reac[2]+reac[3]])
            else:
                print('reaction_list format is not correct.')
                1/0

        self.reaction_list=_reaction_list
        self.reaction_list_noid = [[reac[1], reac[2]]
                                    for reac in self.reaction_list]
        self.reactionNames = [reac[0] for reac in self.reaction_list]

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

        self.stoi = self.make_stoi()
        self.ns = func.compute_rref.compute_rref(linalg.null_space(self.stoi).T).T
        self.ns2 = func.compute_rref.compute_rref(linalg.null_space(self.stoi.T).T)
        self.A = self.M+len(self.ns.T)
        self.graph = [self.cpd_list_noout, self.reaction_list]
        self.cons_list, self.cons_list_index=self.make_conslist()

        self.reac_cons_list=self.reaction_list+self.cons_list



    def info(self):
        print(f'M = {self.M}')
        print(f'R = {self.R}')
        
        if 'out' in self.cpd_list:
            print('outあり')
        else:
            print('outなし')
        
        print('cyc = ', len(self.ns.T))
        print('cons = ', len(self.ns2))
        
        print('det A = ', np.linalg.matrix_rank(self.compute_amat()))

    def make_stoi(self):
        stoi = np.zeros((self.M, self.R), dtype=float)
        for r in range(self.R):
            for m in range(self.M):
                cpd = self.cpd_list_noout[m]
                cpd_sub = self.reaction_list_noid[r][0].count(cpd)  # subに出現する回数
                cpd_pro = self.reaction_list_noid[r][1].count(cpd)  # proに出現する回数
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

        amat = np.zeros((A, A), dtype=float)
        amat[R:A, :M] = -ns2
        amat[:R, M:A] = -ns

        # rmat作成
        rmat = np.zeros((R, M))
        for r in range(R):
            for m in range(M):
                if cpd_list_noout[m] in reaction_list_noid[r][0]:  # subに含まれる
                    rmat[r, m] = np.random.rand()
        # amat作成
        amat[:R, :M] = rmat

        return amat

    def compute_smat(self):
        return func.compute_smat.compute_smat(self)

    def compute_smat_mean(self, N=10):
        return func.compute_smat.compute_smat_mean(self, N)

    def check_ocomp(self, subg):
        # subgの反応は，reaction_listのindexで与える
        reaction_list = self.reaction_list
        R = self.R

        sub_m_list = subg[0]
        sub_r_list = subg[1]

        # cpdを選んで，そこから伸びる反応
        for cpd in sub_m_list:
            for r in range(R):
                if (cpd in reaction_list[r][1]) and (r not in sub_r_list):
                    # print('output completeでない')
                    return False
        return True

    def compute_cyc(self, sub_r_list):
        # 反応のリスト(index)を与えて，サイクル数を返す
        if len(sub_r_list) == 0:
            num_cyc = 0
        else:
            sub_rstoi = self.stoi[:, sub_r_list]
            num_cyc = len(linalg.null_space(sub_rstoi).T)
        return num_cyc

    def compute_cons(self, sub_m_list):
        # 物質のリストを与えて，保存量数を返す
        if len(self.ns2) == 0:
            num_cons = 0
        else:
            nsub_m_index = []
            for m in range(self.M):
                if self.cpd_list_noout[m] not in sub_m_list:
                    nsub_m_index.append(m)
            if len(nsub_m_index) == 0:  # すべての物質を含む
                num_cons = len(self.ns2)
            else:
                nsubstoi = self.stoi[nsub_m_index, :]
                num_cons = len(self.ns2)-len(linalg.null_space(nsubstoi.T).T)
        return num_cons

    def index_subg(self, subg):
        # 部分グラフの指数を返す関数
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

        # 部分グラフに含まれるサイクル数
        if len(sub_r_index) == 0:  # 反応が含まれない
            num_cyc = 0
        else:
            sub_rstoi = stoi[:, sub_r_index]
            num_cyc = len(linalg.null_space(sub_rstoi).T)

        # 保存量数
        if len(ns2) == 0:
            num_cons = 0
        else:
            nsub_m_index = []
            for m in range(self.M):
                if cpd_list_noout[m] not in sub_m_list:
                    nsub_m_index.append(m)
            if len(nsub_m_index) == 0:  # すべての物質を含む
                num_cons = len(ns2)
            else:
                nsubstoi = stoi[nsub_m_index, :]
                num_cons = len(ns2)-len(linalg.null_space(nsubstoi.T).T)

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

def compute_limitset(network):
    return func.compute_limitset_meansmat.compute_limitset_meansmat(network)

def make_hiergraph(limitset_list):
    return func.make_hiergraph.make_hiergraph(limitset_list)
# %%
def from_csv(path):
    reaction_list=[]
    with open(path, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            reaction=[]
            reaction.append(line[0])
            for elem in line[1:]:
                reaction.append(elem.split(' '))
            reaction_list.append(reaction)
    return ReactionNetwork(reaction_list)