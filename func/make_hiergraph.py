import numpy as np
import networkx as nx
def make_hiermat(limitset_list):
    # return a matrix representing the hierarchy of limitset_list

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

def make_hieredge(limitset_list):
    # return a hierarchy graph as a list of edges

    limitset_list_all = []  # reactionとcompoundを一緒のリストにする
    for lset in limitset_list:
        lset_all = lset[0]+[reac[0] for reac in lset[1]]  # reaction_idを追加
        limitset_list_all.append(lset_all)

    l = len(limitset_list_all)
    hiermat = make_hiermat(limitset_list)

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

def make_hiergraph(limitset_list):
    #return an Agraph of the hierarchy
    hieredge_list=make_hieredge(limitset_list)
    hier = []
    pick = []
    for edge in hieredge_list:
        s = ' '.join(edge[0])
        t = ' '.join(edge[1])
        hier.append([short_name(s), short_name(t)])
        pick.extend(list(edge))

    hier_graph = nx.DiGraph(hier)
    hier_agraph = nx.nx_agraph.to_agraph(hier_graph)
    hier_agraph.layout(prog='dot')
    return hier_agraph

def short_name(name):
    l = len(name)
    N = 40  # N文字に一回改行
    for i in range(l//N):
        name = name[:(i+1)*N]+'\n'+name[(i+1)*N:]
    return name