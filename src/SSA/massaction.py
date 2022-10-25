#fluxを計算する
import numpy as np
from tqdm import tqdm

def make_reguMat(network):
    # reaction and its substrates
    reguMat=[]
    for cpd in network.cpd_list_noout:
        row=[network.reaction_list[r][1].count(cpd) for r in range(network.R)]
        reguMat.append(row)
    reguMat=np.array(reguMat)
    return reguMat

def activation(x):
    return x/(1+x)

def inhibition(x):
    return 1/(1+x)

def compute_flux(network,x,reguMat,params,activate=[],inhibit=[]):
    #activate/inhibit=[(index_substrate, index_reaction),...]
    M=network.M
    R=network.R
    consCoef=np.array([x[m]**reguMat[m] for m in range(M)])
    consProd=np.prod(consCoef, axis=0)
    flux=params*consProd

    for (m,r) in activate:
        flux[r]*=activation(x[m])
    for (m,r) in inhibit:
        flux[r]*=inhibition(x[m])
    
    return flux


def perturb(network, ini, steps, params, perturbed, dt=0.01):
    # numcalc of the network, petrubation to a specific reaction
    # steps=[N1, N2]
    # perturbed=[reactionName, perburbation]
    M = network.M
    R = network.R
    reguMat=make_reguMat(network)

    for i, reac in enumerate(network.reaction_list):
        if reac[0] == perturbed[0]:
            perturbed_index = i
    params2 = params.copy()
    params2[perturbed_index] += perturbed[1]
    ans = np.zeros((steps[1], M))
    ans[0] = ini

    for n in tqdm(range(1, steps[1])):
        if n<steps[0]+1:
            flux = compute_flux(network, ans[n-1], reguMat, params)
        else:
            flux=compute_flux(network, ans[n-1], reguMat, params2)
            
        ans[n] = ans[n-1]+np.dot(network.stoi, flux)*dt
        if np.max(ans[n]) > 1.0e20:
            print(n)
            print('overflow')
            return ans
    
    return ans