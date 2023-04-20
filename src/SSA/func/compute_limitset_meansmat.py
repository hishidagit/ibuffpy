#%%
import numpy as np
from . import compute_smat
def compute_limitset_meansmat(network, N, large_error=True):  # N;smatの計算回数
    """compute limitset from smat_mean
    
    Args: 
        network (ReactionNetwork): ReactionNetwork class
        N (int): number of smat_mean calculation
        large_error (bool, optional): if True, use large error. Defaults to True.
    
    Returns:
        list: list of limitset
    
    Examples:
        >>> limitset_list = ReactionNetwork.compute_limitset(network)

    """

    M = network.M
    R = network.R
    A = network.A  # sieze of smat
    smat_mean = compute_smat.compute_smat_mean(network,N=N,large_error=large_error)

    #conver smat to binary form
    smat_bn=np.where(np.abs(smat_mean)<1.0e-10,False,True)

    #
    cpdidx_all=np.arange(M)
    rcindex_all=np.arange(A)
    reacCons_names=[rc[0] for rc in network.reac_cons_list]
    mat_ocomp=np.array([[cpd in reaccons[1] for reaccons in network.reac_cons_list]
                for cpd in network.cpd_list_noout])
    mat_eff=smat_bn.T[:,:M]

    limitset_list=[]
    for perturbed in np.arange(A):
        eff_rc=np.array([perturbed])
        eff_m=cpdidx_all[np.any(mat_eff[eff_rc],axis=0)]
        eff_m=np.sort(np.unique(eff_m))

        while True:
            eff_rc_new=rcindex_all[np.any(mat_ocomp[eff_m],axis=0)]
            eff_rc_new=np.array([rc for rc in eff_rc_new if rc not in eff_rc])
            if len(eff_rc_new)>0:
                eff_m_new=cpdidx_all[np.any(mat_eff[eff_rc_new],axis=0)]
                eff_m=np.concatenate([eff_m,eff_m_new])
                eff_m=np.sort(np.unique(eff_m))
                eff_rc=np.concatenate([eff_rc,eff_rc_new])
                eff_rc=np.sort(np.unique(eff_rc))
            else:
                break
        lsetcpds=[network.cpd_list_noout[m] for m in eff_m]
        lsetrcs=[network.reac_cons_list[r][0] for r in eff_rc]# name only
        lsetcpds.sort()
        lsetrcs.sort()
        lset=[lsetcpds,lsetrcs]
        if lset not in limitset_list:
            limitset_list.append([lsetcpds,lsetrcs])

    # dict_affect=dict()
    # for rc in network.reac_cons_list:
    #     dict_affect[rc[0]]=choose_m([rc], smat_bn, network)
    
    # dict_ocomp=dict()
    # for cpd in network.cpd_list_noout:
    #     dict_ocomp[cpd]=choose_rc([cpd], network)

    # # compute limitsets from smat
    # # start from each reaction.
    # limitsets_found = []
    # for perturbed in tqdm(network.reac_cons_list):
    #     eff_rc = [perturbed]
    #     eff_m = []
    #     # cpds affected by perturbation on reacs in "eff_rc"
    #     for rc in eff_rc:
    #         eff_m+=dict_affect[rc[0]]
    #     eff_m=sorted(list(set(eff_m)))


    #     while True:
    #         eff_rc2 = []  # new reactions to be added
    #         for m in eff_m:
    #             eff_rc2+=dict_ocomp[m]
    #         eff_rc2=[rc for rc in eff_rc2 if rc not in eff_rc]
            
    #         for rc in eff_rc2:
    #             eff_m+=dict_affect[rc[0]]
    #         eff_m=sorted(list(set(eff_m)))

    #         if eff_rc2 == []:  # no new reaction
    #             break
    #         else:
    #             eff_rc+=eff_rc2
        
    #     #sort reactions
    #     eff_rc_names=sorted([rc[0] for rc in eff_rc])
    #     eff_rc=[network.get_reacCons_by_id(_id) for _id in eff_rc_names]
    #     if (eff_m, eff_rc) not in limitsets_found:
    #         limitsets_found.append((eff_m, eff_rc))

    # sort limitset_list
    limitset_list=sorted(limitset_list, key=lambda x: len (x[0])+len (x[1]))

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
