
import numpy as np
import sympy as sp
from string import ascii_uppercase

sp.init_printing()

'''
A symbolic toolbox for structural sensitivity analysis.
Users are suggested to use IPython when working with this toolbox.
One may consider Jupyter Notebook or Google Colab, in particular.
Besides, please keep in mind that basically the stoichiometric matrices are assumed to possess integer entries only.

As an example, a beginner may first give it a shot with the following codes:


nu = np.array([[  1, -1,  1, -1],
               [  0,  1, -1,  0],
               [  0,  1, -1,  0]])
X, R = vec_XandR(*nu.shape)
print("Given the ODE of a chemical reaction network")
display(sp.Eq(
    sp.symbols(r'\frac{d}{dt}'+sp.latex(X)),
    sp.symbols((sp.latex(sp.Matrix(nu)) + sp.latex(R)).replace(" ","")),
evaluate = False))

print("")
print("The the augmented matrix can given by")
A   = AugMat(nu, kineticstype = "<")
display(sp.Eq(sp.symbols(r'\bm{A}'), A, evaluate = False))
print("In a more genral kinetic model, the matrix may read")
A   = AugMat(nu, kineticstype = "!=")
display(sp.Eq(sp.symbols(r'\bm{A}'), A, evaluate = False))
print("and its determinants is obtained as")
display(sp.Eq(sp.symbols(r'\det{\bm{A}}'), sp.factor(A.det()), evaluate = False))

print("")
print("Meanwhile, the Jacobian matrix is")
J_f = nu @ A[:nu.shape[1], :nu.shape[0]]
display(sp.Eq(sp.symbols(r'J_f'), J_f, evaluate = False))
print("with its determinant being")
display(sp.Eq(sp.symbols(r'\det{J_f}'), sp.factor(J_f.det()), evaluate = False))

'''

def KerImg(mat):
    """
    To find bases of the image and null spaces for a given matrix by means of Gaussian elimination.
    Inputs: 
        mat
            2D array-like, of which entries are integers.
    Outputs:
        a list with 2 elements, in which
        [0] A matrix of which column vectors form a basis of the kernel. 
        [1] A matrix of which column vectors form a basis of the image.
    """
    M, N = mat.shape
    def _gaussian(ErMat, m, n):
        "Gaussian elimination on the m-th row for except for the first n's columns."
        headrow  = ErMat[N+m,n:]
        minval   = np.min(np.abs(headrow)[headrow!=0])
        k        = (np.arange(N-n)[np.abs(headrow)==minval])[0]
        if k != 0:
            ErMat[:, [n,n+k]] = ErMat[:, [n+k,n]] # swap
        for i in np.arange(1, N-n):
            if headrow[i] != 0:
                lcm = np.lcm(headrow[0], headrow[i])
                ErMat[:,n+i] = ((lcm/headrow[i]) * ErMat[:,n+i]).astype(int)
                ErMat[:,n+i] = ErMat[:,n+i] \
                        - ((lcm/headrow[0]) * ErMat[:,n]).astype(int)
        return ErMat.astype(int)
    ermat = np.concatenate([np.eye(N), mat.copy()], axis = 0).astype(int)
    # Computation for a basis of the kernel
    n = 0
    for n in range(N):
        mlist = (np.arange(M)[np.any(ermat[N:,n:], axis = 1)])
        if len(mlist) == 0:
            break
        ermat = _gaussian(ermat, mlist[0], n)
        n = n+1
    # Reduction of the basis of the image
    # this step is actually not necessary for the construction of the basis
    m = 0
    Q = np.sum(np.any(ermat[N:,:], axis = 0).astype(int))
    for q in range(Q):
        if (m >= M):
            break
        elif ermat[N+m, q] == 0:
            m = m + 1
            pass
        for i in range(Q):
            if (i == q):
                ermat[:,i] = np.sign(ermat[N+m,i]) * ermat[:,i]
                pass
            else:
                ermat[:,i] = ermat[:,i] \
                        - (ermat[N+m,i]/ermat[N+m,q]) * ermat[:,q]
    indice = ~np.any(ermat[N:,:], axis = 0)
    return ermat[:N, indice], ermat[N:, ~indice]

def AugMat(nu, kineticstype = "!="):
    """
    To construct a symbolic augmented matrix A given a stoichiometric matrix.
    Input:
        nu
            2D array-like, the stoichiometric matrix, of which entries are required to be integers.
        kineticstype
            string, which specifies the condition in which partial derivatives of r is nonzero.
            default to be "!=", which means (partial r_j)/(partial x_i) is not cancelled to be zero
            if nu_{ij} != 0.
    Output:
        spA
            the symbolic augmented matrix A.
    """
    col_stand = lambda mat: (mat / np.gcd.reduce(mat, axis = 0))
    C   = col_stand(KerImg(nu)[0])
    D   = col_stand(KerImg(nu.T)[0])
    dim = nu.shape[0]+C.shape[1]
    spA = sp.Matrix(np.zeros([dim, dim]).astype(int))
    if np.all(C.shape):
        spA[:C.shape[0],-C.shape[1]:] = sp.Matrix(C.astype(int))
    if np.all(D.shape):
        spA[-D.shape[1]:,:D.shape[0]] = sp.Matrix(D.T.astype(int))
    for i in range(nu.shape[0]):
        for j in range(nu.shape[1]):
            if eval("{} {} 0".format(nu[i,j], kineticstype)):
                spA[j,i] = sp.symbols('r_%s%s'%(j+1,ascii_uppercase[i]))
    return spA

def vec_XandR(M, N):
    """
    A quick construction of vectors representing the system state and reaction rates.
    Input:
        M
            A natural number, the number of the chemicals.
        N
            A natural number, the number of the reactions.
    Output:
        a list with 2 elements, in which
        [0] The system state X
        [1] The reaction rate R
    """
    X = [[sp.symbols(ascii_uppercase[i]) for i in range(M)]]
    R = [[sp.symbols('r_{}'.format(i)) for i in range(1, N+1)]]
    return sp.Matrix(X).T, sp.Matrix(R).T

## == END OF THE SCRIPT == ##
