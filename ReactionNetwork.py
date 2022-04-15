# %%

import numpy as np
from scipy import linalg
import networkx as nx
import func.compute_limitset_meansmat
import func.compute_rref
import func.make_hiergraph

class ReactionNetwork:

    def __init__(self, reaction_list, info=True):
        if info:
            print('constructed.')
        if len(reaction_list[0]) == 3:
            self.reaction_list = reaction_list
            self.reaction_list_noid = [[reac[1], reac[2]]
                                       for reac in reaction_list]
        else:
            print('reaction_list format is not correct.')
            1/0

        self.reactionNames = [reac[0] for reac in reaction_list]

        cpd_list = []
        for rn in self.reaction_list_noid:
            for form in [rn[0], rn[1]]:
                for cpd in form:
                    if not cpd in cpd_list:
                        cpd_list.append(cpd)
        cpd_list = sorted(cpd_list)
        self.cpd_list = cpd_list

        self.R = len(self.reaction_list_noid)
        self.M = len(self.cpd_list)
        if 'out' in self.cpd_list:
            self.M -= 1
            cpd_list_noout = cpd_list.copy()
            cpd_list_noout.remove('out')
            if info:
                print('outあり')
        else:
            cpd_list_noout = cpd_list.copy()
            if info:
                print('outなし')
        if info:
            print('M = ', self.M)
            print('R = ', self.R)

        self.cpd_list_noout = cpd_list_noout
        self.stoi = self.make_stoi()
        self.ns = linalg.null_space(self.stoi)
        self.ns2 = func.compute_rref.compute_rref(linalg.null_space(self.stoi.T).T)
        self.A = self.M+len(self.ns.T)
        self.graph = [cpd_list_noout, reaction_list]
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
        R = self.R
        M = self.M

        ns = self.ns
        ns2 = self.ns2

        A = R+len(ns2)
        if R+len(ns2) != M+len(ns.T):
            print('A行列が正方行列でない')
            1/0

        amat = self.compute_amat()
        # smat計算
        smat = np.linalg.inv(amat)

        return smat

    def compute_smat_mean(self, N=10):

        # smatをN回計算して，平均を取る
        smat_all = np.array([self.compute_smat() for i in range(N)])
        smat_mean = np.mean(smat_all, axis=0)

        # 閾値が1.0e-10なので，その周辺の値があると誤差の影響があるかもしれない
        np_mean_check = np.where(
            (smat_mean < 1.0e-8) & (smat_mean > 1.0e-10), 1, 0)
        if np.sum(np_mean_check) == 0.0:
            0
        else:
            print('large error')

        return smat_mean

    def find_limitset_meansmat(self, N=10):  # N;smatの計算回数

        M = self.M
        R = self.R
        cpd_list_noout = self.cpd_list_noout
        reaction_list = self.reaction_list
        A = R+len(self.ns2)  # smatのサイズ

        smat_mean = self.compute_smat_mean(N=N)

        # smat_meanをもとにしてlimitsetを出す
        
        #smatをbinaryに変換
        smat2 = np.zeros((A, A), dtype='int')
        for i in range(A):
            for j in range(A):
                if abs(smat_mean[i, j]) < 1.0e-10:
                    smat2[i, j] = 0
                else:
                    smat2[i, j] = 1

        # 限局集合を求めるのに使う関数
        # 注意:次の関数はsmat2, cpd_list_noout, raction_list3に依存する

        def choose_m(r_list, smat2, cpd_list_noout):  # 反応を与えて、それが影響する物質のリストを返す
            m_list = []
            for reac in r_list:
                r = reaction_list.index(reac)
                for m in range(M):
                    if smat2[m, r] == 1:
                        m_list.append(cpd_list_noout[m])
            m_list = list(set(m_list))
            return m_list

        def choose_r(m_list):  # 物質のリストを与えて、outputcompleteになるための反応のリストを返す
            r_list = []
            for reac in reaction_list:
                if set(reac[1]) & set(m_list) != set():
                    r_list.append(reac)
            return r_list

        # smat2から限局集合を出す
        # 各反応から始める。反応はindexで扱う
        limitsets_found = []
        for reac in reaction_list:
            eff_m = []
            eff_r = [reac]
            eff_r2 = []  # 新たに追加する反応のリスト

            while True:
                # eff_rが影響する物質
                # 注意:次の関数はsmat2, cpd_list_noout, raction_listに依存する
                eff_m = sorted(list( set(eff_m) | set(choose_m(eff_r, smat2, cpd_list_noout)) ))
                eff_r2 = [reac for reac in choose_r(
                    eff_m) if (reac not in eff_r)]
                if eff_r2 == []:  # 新しく追加する反応がない
                    break

                eff_r.extend(eff_r2)

            if (eff_m, eff_r) not in limitsets_found:
                limitsets_found.append((eff_m, eff_r))

        # limitset_listを並びかえ
        hoge = []
        i = 0
        while i < M+R+1:
            for lset in limitsets_found:
                if len(lset[0])+len(lset[1]) == i:
                    hoge.append(lset)
            i += 1
        limitset_list = hoge

        return limitset_list

    def make_hiermat(self, limitset_list):
        # return a matrix representing hierarchy

        limitset_list_all = []
        for lset in limitset_list:
            lset_all = lset[0]+[reac[0] for reac in lset[1]]  # reaction_idを追加
            limitset_list_all.append(lset_all)

        l = len(limitset_list_all)

        relmat = np.zeros((l, l), dtype=int)  # limitsetの包含関係を表す行列
        # 列<行のときにrelmatの要素を1にする
        for i, s in enumerate(limitset_list_all):
            for j, t in enumerate(limitset_list_all):
                if set(t) < set(s):
                    relmat[i, j] = 1

        relmat2 = []  # さらに下流
        for i in range(l):
            row = np.dot(relmat[i], relmat)
            row2 = np.array([1 if p > 0 else 0 for p in row])
            relmat2.append(row2)
        relmat2 = np.array(relmat2)

        hiermat = relmat-relmat2
        if np.min(hiermat) < 0:
            print('包含関係が成り立たない')
            1/0
        return hiermat

    def make_hieredge(self, limitset_list):
        # return a hierarchy graph as a list of edges

        limitset_list_all = []  # reactionとcompoundを一緒のリストにする
        for lset in limitset_list:
            lset_all = lset[0]+[reac[0] for reac in lset[1]]  # reaction_idを追加
            limitset_list_all.append(lset_all)

        l = len(limitset_list_all)
        hiermat = self.make_hiermat(limitset_list)

        # ヒエラルキーグラフのノードを作成
        node_list = []
        for i in range(l):
            elim = set()
            for j in range(l):
                if hiermat[i, j] == 1:  # jがiに含まれる
                    elim |= set(limitset_list_all[j])  # iのノードからjを除く
            node = sorted(list(set(limitset_list_all[i])-elim))

            node_list.append(node)

        # ヒエラルキーのエッジを作成
        edge_list = []
        for i, s in enumerate(node_list):  # source
            for j, t in enumerate(node_list):  # target
                if hiermat[i, j] == 1:
                    edge_list.append((s, t))

        return edge_list

    def make_hiergraph(self, hieredge_list):
        hier = []
        pick = []
        for edge in hieredge_list:
            s = ' '.join(edge[0])
            t = ' '.join(edge[1])
            hier.append([self.short_name(s), self.short_name(t)])
            pick.extend(list(edge))

        hier_graph = nx.DiGraph(hier)
        hier_agraph = nx.nx_agraph.to_agraph(hier_graph)
        hier_agraph.layout(prog='dot')
        return hier_agraph

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