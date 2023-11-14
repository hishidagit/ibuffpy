import numpy as np
from scipy import linalg
import scipy.sparse.linalg as sla

def rowreduce2(matrix, tol=1.0e-10):
    pivot = []
    mat=np.array(matrix, dtype = float)
    if len(mat)==0:
        return mat, pivot #empty matrix
    #rank=np.linalg.matrix_rank(mat)
    row, col= mat.shape
    #kはreduceが終わった行の数。
    for k in range(row):
        # the largest row index with pivot is k-1
        # focus on rows from k to the last row
        #i:row j:col
        i = k
        
        # if all elements from k to last row is zero, return the matrix
        if np.max(np.abs(mat[k:,:])) < tol:
            return mat, pivot
        
        # j is the first column index with non-zero element in rows from k to the last row
        sum_k_to_last_row = np.sum(np.abs(mat[k:,:]), axis=0)
        j = np.where(sum_k_to_last_row > tol)[0][0]

        pivot.append (j)#j列がpivotになることは確定。
        #j列目の、k行以下で、非ゼロ成分を最も上に持つ行をi行とする
        # the row index with nonzero element in the j-th column from k to the last row is set to be a pivot
        i = k + np.argmax(np.abs(mat[k:,j])) 

        #i行とk行を入れ替える
        # swap i-th row and k-th row
        mat[[k,i]]=mat[[i,k]]

        #ピボットを1に
        # make the pivot element 1
        mat[k]=mat[k]/mat[k,j]

        #j列のk行目以外の要素を0にする
        for r in range(row):
            if r==k:
                continue
            mat[r]=mat[r]-mat[k]*mat[r,j]
        # refine zero elements
        mat = np.where(np.abs(mat) < tol, 0, mat)
    return mat, pivot



def cal_nullspace_rref (matrix, tol=1.0e-10):
    mat = np.array(matrix, dtype = float)
    # calculate nullspace of matrix numerically and convert it to rref (this version may faster?)
    ker_svd = cal_nullspace_svd(mat, tol=tol)
    ker_rref = rowreduce2(ker_svd.T, tol=tol)[0].T
    if len(mat) ==0:
        return ker_rref
    if np.max(np.abs(mat @ ker_rref)) > tol:
        raise ValueError('Error: rref calculation of nullspace vectors have too large error.')
    return ker_rref

    # calculate nullspace of matrix by converting it to RREF
    n_row, n_col = mat.shape
    rref, pivot=rowreduce2(mat, tol=tol)
    non_pivot = [i for i in range (n_col) if i not in pivot]#pivotでない列のindex

    ker =[]
    for i in non_pivot:
        vec = np.zeros (n_col)
        vec [i] =1.
        for j, s in enumerate (pivot):
            vec[s] = -rref[j][i]
        ker.append (vec)
    ker=np.array (ker)

    #check dimension
    rank=np.linalg.matrix_rank(mat)
    dim_ker=len (ker)
    assert rank+dim_ker==n_col, f"The nullspace dimension might be wrong. rank={rank}, dim_ker={dim_ker}, col={n_col}"

    #check product
    if dim_ker>0:
        dot_pro = np.dot (mat, ker.T)
        min_pro= np.min (np.abs (dot_pro))
        assert min_pro < tol, "Nullspace bases might be wrong."


    return ker.T#Converting ker into a column vector

#calculate nullspace usig SVD
def cal_nullspace_svd (matrix, tol=1.0e-10):
    mat = np.array (matrix, dtype = float)

    # when network size is small, ns(cycles) are computed using null_space
    if np.max(mat.shape) < 1000:
        try :
            ns=linalg.null_space(mat)
            return ns
        except np.linalg.LinAlgError:
            pass

    
    # when network size is large, ns(cycles) are computed using truncated SVD
    # dimension of nullspace of self.stoi
    # from sklearn.decomposition import TruncatedSVD
    # compute SVD of mat.T * mat
    mat_hermitian = mat.T @ mat
    _,s,vh = np.linalg.svd(mat_hermitian,hermitian=True)
    # _,s,vh = sla.svds(mat.T@mat,which='SM', k=mat.shape[1]-1,solver='arpack')
    # _,s,vh = sla.svds(mat.T@mat,which='SM', k=mat.shape[1]-1,solver='lobpcg')
    nullspace_mask = np.abs(s) < tol
    nullspace = vh[nullspace_mask].T
    # svd = TruncatedSVD(n_components=mat.shape[1],algorithm='arpack')
    # svd.fit(mat.T @ mat)
    # # vectors of right singular vectors correspoinding to small singular values
    # nullspace_mask = svd.singular_values_ < tol
    # nullspace = svd.components_[nullspace_mask].T
    
    # check if nullspace estimation is correct
    error = np.max(np.abs(mat @ nullspace))

    # when error is large, scipy.linalg.null_space is used
    if error > tol:
        nullspace = linalg.null_space(mat)
        error=np.max(np.abs(mat @ nullspace))
    if error > tol:
        raise Exception('Error: nullspace estimation by svd have too large error:', error)
    else:
        ns=nullspace
    return ns#a column vector

#def cal_nullspace_sympy (matrix):
    #import sympy
    #from sympy.polys.matrices import DomainMatrix
    #mat = np.array (matrix, dtype = float)
    #dM = DomainMatrix.from_Matrix(sympy.Matrix(matrix))
    #ns = dM.to_field().nullspace()
    #ns = ns.__dict__['rep']
    #ns = np.array ([ns[k] for k in range (len(ns))], dtype=float).T
    #return ns#a column vector


############## Test codes ##########################################################
if __name__=="__main__":
    ########## Example 1 ########
    matrix_eg1 = np.array([[1, 2, 3, 0, 0], [4, 10, 0, 0, 1]])
    ker_eg1_rref=cal_nullspace_rref (matrix_eg1)#calculate by RREF
    print ("eg1, rref \n", ker_eg1_rref)
    ##[[-15.   -0.    1. ]
     ##[  6.   -0.   -0.5]
     ##[  1.    0.    0. ]
     ##[  0.    1.    0. ]
     ##[  0.    0.    1. ]]
    print ("eg1, rref-check \n", matrix_eg1 @ker_eg1_rref)#check the result
    ## [[0. 0. 0.]
      ##[0. 0. 0.]]

    ker_eg1_svd=cal_nullspace_svd (matrix_eg1)#calculate by SVD
    print ("eg1, svd \n", ker_eg1_svd)
    ## [[-0.92700459  0.          0.01900451]
      ##[ 0.36551773  0.         -0.10679985]
      ##[ 0.06532305  0.          0.06486506]
      ##[ 0.          1.          0.        ]
      ##[ 0.05284109  0.          0.99198042]]

    print ("eg1, svd-check \n", matrix_eg1 @ker_eg1_svd)#check the result
    ##[[-7.77156117e-16  0.00000000e+00  0.00000000e+00]
     ##[-8.60422844e-16  0.00000000e+00  2.22044605e-16]]


    ########## Example 2 ########
    matrix_eg2 = [[1, 1, 1, 2],[1, -1,-1, 1],[1, 3, 3, 3],[3, 1, 1, 5]]
    ker_eg2_rref=cal_nullspace_rref (matrix_eg2)#calculate by RREF
    print ("eg2, rref \n", ker_eg2_rref)
    ## [[-0.  -1.5]
      ##[-1.  -0.5]
      ##[ 1.   0. ]
      ##[ 0.   1. ]]
    print ("eg2, rref-check \n", np.array(matrix_eg2) @ker_eg2_rref)#check the result
     ##[[0. 0.]
     ##[0. 0.]
     ##[0. 0.]
     ##[0. 0.]]

    ker_eg2_svd=cal_nullspace_svd (matrix_eg2)#calculate by SVD
    print ("eg2, rref \n", ker_eg2_svd)
    ##[[ 0.52344428 -0.62663606]
    ##[ 0.62992346  0.3488767 ]
    ##[-0.45544203 -0.55775539]
    ##[-0.34896286  0.41775737]]
    print ("eg2, svd-check \n", np.array(matrix_eg2) @ker_eg2_svd)#check the result
     ##[[-2.22044605e-16  2.22044605e-16]
     ##[ 3.33066907e-16  6.66133815e-16]
     ##[-7.77156117e-16 -4.44089210e-16]
     ##[-4.44089210e-16  8.88178420e-16]]

    ########## Example 3 ########
    matrix_eg3 = np.array([[1, 2, 3], [4, 10, 1],[11, 10, 8]])
    ker_eg3_rref=cal_nullspace_rref (matrix_eg3)#calculate by RREF
    print ("eg3, rref \n", ker_eg3_rref)


    ker_eg3_svd=cal_nullspace_svd (matrix_eg3)#calculate by SVD
    print ("eg3, svd \n", ker_eg3_svd)
    #[]
