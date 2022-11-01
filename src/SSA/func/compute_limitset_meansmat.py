#%%
import numpy as np
from tqdm import tqdm
from . import compute_smat
def compute_limitset_meansmat(network, N, large_error=True):  # N;smatの計算回数

    M = network.M
    R = network.R
    A = R+len(network.ns2)  # sieze of smat
    smat_mean = compute_smat.compute_smat_mean(network,N=N,large_error=large_error)

    #conver smat to binary form
    smat_bn = np.zeros((A, A), dtype='bool')
    for i in range(A):
        for j in range(A):
            if abs(smat_mean[i, j]) < 1.0e-10:
                smat_bn[i, j] = False
            else:
                smat_bn[i, j] = True

    dict_affect=dict()
    for rc in network.reac_cons_list:
        dict_affect[rc[0]]=choose_m([rc], smat_bn, network)
    
    dict_ocomp=dict()
    for cpd in network.cpd_list_noout:
        dict_ocomp[cpd]=choose_rc([cpd], network)

    # compute limitsets from smat
    # start from each reaction.
    limitsets_found = []
    for perturbed in tqdm(network.reac_cons_list):
        eff_rc = [perturbed]
        eff_m = []
        # cpds affected by perturbation on reacs in "eff_rc"
        for rc in eff_rc:
            eff_m+=dict_affect[rc[0]]
        eff_m=sorted(list(set(eff_m)))


        while True:
            eff_rc2 = []  # new reactions to be added
            for m in eff_m:
                eff_rc2+=dict_ocomp[m]
            eff_rc2=[rc for rc in eff_rc2 if rc not in eff_rc]
            
            for rc in eff_rc2:
                eff_m+=dict_affect[rc[0]]
            eff_m=sorted(list(set(eff_m)))

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
# 注意:次の関数はsmat_bn, cpd_list_noout, raction_list3に依存する

def choose_m(perturbed_list, smat_bn, network):
    # 反応を与えて、それが影響する物質のリストを返す
    m_list = []
    for perturbed in perturbed_list:
        #index of the reaction/cons.
        p=network.reac_cons_list.index(perturbed)
        for m in range(network.M):
            if smat_bn[m, p]:
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
