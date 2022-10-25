#%%
import numpy as np
from . import compute_smat
def compute_limitset_meansmat(network, N, large_error=True):  # N;smatの計算回数

    M = network.M
    R = network.R
    cpd_list_noout = network.cpd_list_noout
    reaction_list = network.reaction_list
    A = R+len(network.ns2)  # smatのサイズ

    smat_mean = compute_smat.compute_smat_mean(network,N=N,large_error=large_error)

    # smat_meanをもとにしてlimitsetを出す
    
    #smatをbinaryに変換
    smat2 = np.zeros((A, A), dtype='int')
    for i in range(A):
        for j in range(A):
            if abs(smat_mean[i, j]) < 1.0e-10:
                smat2[i, j] = 0
            else:
                smat2[i, j] = 1



    # smat2から限局集合を出す
    # 各反応から始める。反応はindexで扱う
    limitsets_found = []
    for perturbed in network.reac_cons_list:
        eff_m = []
        eff_rc = [perturbed]
        eff_rc2 = []  # 新たに追加する反応のリスト

        while True:
            # eff_rcが影響する物質
            eff_m = sorted(list( set(eff_m) | set(choose_m(eff_rc, smat2, network)) ))


            eff_rc2 = [reac for reac in choose_rc(eff_m, network) if (reac not in eff_rc)]
            if eff_rc2 == []:  # 新しく追加する反応がない
                break

            eff_rc.extend(eff_rc2)

        if (eff_m, eff_rc) not in limitsets_found:
            limitsets_found.append((eff_m, eff_rc))

    # limitset_listを並びかえ
    limitset_list=sorted(limitsets_found, key=lambda x: len (x[0])+len (x[1]))

    return limitset_list

# 限局集合を求めるのに使う関数
# 注意:次の関数はsmat2, cpd_list_noout, raction_list3に依存する

def choose_m(perturbed_list, smat2, network):
    # 反応を与えて、それが影響する物質のリストを返す
    m_list = []
    for perturbed in perturbed_list:
        #index of the reaction/cons.
        p=network.reac_cons_list.index(perturbed)
        for m in range(network.M):
            if smat2[m, p] == 1:
                m_list.append(network.cpd_list_noout[m])
    m_list = list(set(m_list))
    return m_list

def choose_rc(m_list, network):
    # 物質のリストを与えて、outputcompleteになるための反応と保存量のリストを返す
    rc_list = []
    for reacCons in network.reac_cons_list:
        #if one of metabolites is included in sub(reaction) or cons., 
        # then the reac/cons is included in the perturbed_list
        if set(reacCons[1]) & set(m_list) != set():
            rc_list.append(reacCons)
    return rc_list
# %%
