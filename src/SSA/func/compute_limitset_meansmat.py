#%%
import numpy as np
from . import compute_smat
def compute_limitset_meansmat(network, N, large_error=True):  # N;smatの計算回数

    M = network.M
    R = network.R
    cpd_list_noout = network.cpd_list_noout
    reaction_list = network.reaction_list
    A = R+len(network.ns2)  # sieze of smat

    smat_mean = compute_smat.compute_smat_mean(network,N=N,large_error=large_error)

    
    #conver smat to binary form
    smat2 = np.zeros((A, A), dtype='int')
    for i in range(A):
        for j in range(A):
            if abs(smat_mean[i, j]) < 1.0e-10:
                smat2[i, j] = 0
            else:
                smat2[i, j] = 1



    # compute limitsets from smat
    # start from each reaction.
    limitsets_found = []
    for perturbed in network.reac_cons_list:
        eff_m = []
        eff_rc = [perturbed]
        eff_rc2 = []  # new reactions to be added

        while True:
            # cpds affected perturbation on reacs in "eff_rc"
            eff_m=sorted(list(set(eff_m + choose_m(eff_rc,smat2,network))))
            eff_rc2 = [reac for reac in choose_rc(eff_m, network) if (reac not in eff_rc)]
            if eff_rc2 == []:  # no new reaction
                break
            else:
                eff_rc+=eff_rc2
        
        #sort reactions
        eff_rc_names=sorted([rc[0] for rc in eff_rc])
        eff_rc=[network.get_reacCons_by_id(_id) for _id in eff_rc_names]
        if (eff_m, eff_rc) not in limitsets_found:
            limitsets_found.append((eff_m, eff_rc))

    # sort limitset_list
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
    # given list of cpds, return list of reacs and cons for output-complete
    rc_list = []
    for reacCons in network.reac_cons_list:
        if not set(reacCons[1]).isdisjoint(set(m_list)):
            rc_list.append(reacCons)
    return rc_list
# %%
