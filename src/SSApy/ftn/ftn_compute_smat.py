import numpy as np

def compute_smat(network):
    R = network.R
    M = network.M

    ns = network.ns
    ns2 = network.ns2

    A = R+len(ns2)
    if R+len(ns2) != M+len(ns.T):
        raise Exception('A matrix is not square.')

    amat = network.compute_amat()
    # computing S matrix
    # use np.linalg.solve instead of np.linalg.inv
    # smat = -np.linalg.inv(amat)
    smat = -np.linalg.solve(amat,np.eye(network.A))

    return smat

def compute_smat_mean(network, N, large_error=True):

    # calculate smat for N times, get mean

    # #when network is large, np.mean requires too much memory
    # smat_all = np.array([compute_smat(network) for i in range(N)])
    # smat_mean = np.mean(smat_all, axis=0)

    smat_sum=np.zeros((network.A, network.A))
    for n in range(N):
        smat_sum+=compute_smat(network)
    smat_mean=smat_sum/N

    # check error size
    np_mean_check = np.where(
        (np.abs(smat_mean) < 1.0e-8) & (np.abs(smat_mean) > 1.0e-10), 1, 0)
    if np.sum(np_mean_check) == 0.0:
        0
    else:
        if large_error:
            raise LargeErrorSmat('smat_mean have large error.')
        else:
            print('large error warning')

    return smat_mean

##############Determine signs of the sensitivity##############
def determine_sign (val):
    if val>10**(-10):
        return "+"
    if val<-10**(-10):
        return "-"
    else:
        return "0"

def fun_integrate (val1,val2):
    if val1==val2:
        return val1
    else:
        return "+/-"

def compute_smat_sign (network, N):
    determine_sign_vec = np.vectorize(determine_sign)
    fun_integrate_vec = np.vectorize(fun_integrate)
    smat= network.compute_smat()
    result_array = determine_sign_vec(smat)#Detrmine signs of Smat.
    for _ in range (N-1):
        smat=network.compute_smat()
        sign_array_now=determine_sign_vec(smat)#Detrmine signs of Smat.
        #For each entry of Smat, if check if the sign is equal to the sign of previous result
        result_array=np.array ([fun_integrate_vec (x,y) for (x, y) in zip(result_array, sign_array_now)])
    return result_array

class LargeErrorSmat(Exception):
    pass
