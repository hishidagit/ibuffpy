import numpy as np

def compute_smat(network):
    R = network.R
    M = network.M

    ns = network.ns
    ns2 = network.ns2

    A = R+len(ns2)
    if R+len(ns2) != M+len(ns.T):
        print('A行列が正方行列でない')
        1/0

    amat = network.compute_amat()
    # smat計算
    smat = np.linalg.inv(amat)

    return smat

def compute_smat_mean(network, N=10):

    # calculate smat for N times, get mean
    smat_all = np.array([compute_smat(network) for i in range(N)])
    smat_mean = np.mean(smat_all, axis=0)

    # check error size
    np_mean_check = np.where(
        (smat_mean < 1.0e-8) & (smat_mean > 1.0e-10), 1, 0)
    if np.sum(np_mean_check) == 0.0:
        0
    else:
        print('large error')

    return smat_mean