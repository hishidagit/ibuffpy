#%%
import numpy as np
from . import compute_smat
def compute_limitset_meansmat(network, N, large_error=True, detectCQ=True):  # N;smatの計算回数
    """compute limitsets from smat_mean.
    
    Args: 
        network (ReactionNetwork): ReactionNetwork class
        N (int): number of smat_mean calculation
        large_error (bool, optional): if True, use large error. Defaults to True.
        detectCQ (bool, optional): if True, detect conserved quantities. Defaults to True.
    
    Returns:
        list: list of limitset
    
    Examples:
        >>> limitset_list = ReactionNetwork.compute_limitset(network)
    """

    if not detectCQ:
        return compute_limitset_meansmat_noCQ(network, N, large_error)
    
    # compute mean of S-matrix
    smat_mean = compute_smat.compute_smat_mean(network,N=N,large_error=large_error)

    #conver smat to binary form by threshold
    smat_bn=np.where(np.abs(smat_mean)<network.tol,False,True)

    # reactions and compounds are treated as index to use numpy
    cpdidx_all=np.arange(network.M)
    rcindex_all=np.arange(network.A)
    # the matrix representing output complete and conserved quantities
    mat_ocomp=np.array([[cpd in reaccons[1] for reaccons in network.reac_cons_list]
                for cpd in network.cpd_list_noout])
    # effect of perturbation ignoring cycle part of S-matrix
    mat_eff=smat_bn.T[:,:network.M]

    # compute limitset, starting from each reaction
    limitset_list=[]
    for perturbed in np.arange(network.A):
        # reactions and CQs in a limitset
        eff_rc=np.array([perturbed])
        # compounds in a limitset, detected from smat as the first step
        eff_m=cpdidx_all[np.any(mat_eff[eff_rc],axis=0)]
        eff_m=np.sort(np.unique(eff_m))

        # select new reactions, CQs and compounds until there is no reaction or CQ to be added
        while True:
            # select new reactions to satisfy output-completeness and CQs constructed from some of selected compounds
            eff_rc_new=rcindex_all[np.any(mat_ocomp[eff_m],axis=0)]
            eff_rc_new=np.array([rc for rc in eff_rc_new if rc not in eff_rc])
            if len(eff_rc_new)>0:
                # select new compounds from S-matrix
                eff_m_new=cpdidx_all[np.any(mat_eff[eff_rc_new],axis=0)]
                eff_m=np.concatenate([eff_m,eff_m_new])
                eff_m=np.sort(np.unique(eff_m))
                eff_rc=np.concatenate([eff_rc,eff_rc_new])
                eff_rc=np.sort(np.unique(eff_rc))
            else:
                break
        
        # convert compounds and reactions from index to names
        lsetcpds=[network.cpd_list_noout[m] for m in eff_m]
        lsetrcs=[network.reac_cons_list[r][0] for r in eff_rc]# name only
        lsetcpds.sort()
        lsetrcs.sort()
        lset=[lsetcpds,lsetrcs]
        if lset not in limitset_list:
            limitset_list.append([lsetcpds,lsetrcs])

    # sort limitset by the number of compounds and reactions
    limitset_list=sorted(limitset_list, key=lambda x: len(x[0])+len(x[1]))

    return limitset_list

def compute_limitset_meansmat_noCQ(network, N, large_error=True):
    """compute limitsets from smat_mean without detecting conserved quantities."""
    
    smat_mean = compute_smat.compute_smat_mean(network,N=N,large_error=large_error)

    #conver smat to binary form by threshold
    smat_bn=np.where(np.abs(smat_mean)<network.tol,False,True)

    # reactions and compounds are treated as index to use numpy
    cpdidx_all=np.arange(network.M)
    rxnindex_all=np.arange(network.R) # only focus on reactions
    # the matrix representing output complete and conserved quantities
    mat_ocomp=np.array([[cpd in rxn[1] for rxn in network.reaction_list]
                for cpd in network.cpd_list_noout])
    # effect of perturbation ignoring cycle part of S-matrix
    mat_eff=smat_bn.T[:network.R,:network.M]

    # compute limitset, starting from each reaction
    limitset_list=[]
    for perturbed in np.arange(network.R):# only reactions
        # reactions and CQs in a limitset
        eff_r=np.array([perturbed])
        # compounds in a limitset, detected from smat as the first step
        eff_m=cpdidx_all[np.any(mat_eff[eff_r],axis=0)]
        eff_m=np.sort(np.unique(eff_m))

        # select new reactions and compounds until there is no reaction to be added
        while True:
            # select new reactions to satisfy output-completeness and CQs constructed from some of selected compounds
            eff_r_new=rxnindex_all[np.any(mat_ocomp[eff_m],axis=0)]
            eff_r_new=np.array([r for r in eff_r_new if r not in eff_r])
            if len(eff_r_new)>0:
                # select new compounds from S-matrix
                eff_m_new=cpdidx_all[np.any(mat_eff[eff_r_new],axis=0)]
                eff_m=np.concatenate([eff_m,eff_m_new])
                eff_m=np.sort(np.unique(eff_m))
                eff_r=np.concatenate([eff_r,eff_r_new])
                eff_r=np.sort(np.unique(eff_r))
            else:
                break
        
        # convert compounds and reactions from index to names
        lsetcpds=[network.cpd_list_noout[m] for m in eff_m]
        lsetrxns=[network.reac_cons_list[r][0] for r in eff_r]# name only
        lsetcpds.sort()
        lsetrxns.sort()
        lset=[lsetcpds,lsetrxns]
        if lset not in limitset_list:
            limitset_list.append([lsetcpds,lsetrxns])

    # sort limitset by the number of compounds and reactions
    limitset_list=sorted(limitset_list, key=lambda x: len(x[0])+len(x[1]))

    return limitset_list