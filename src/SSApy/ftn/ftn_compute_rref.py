import numpy as np

def rowreduce(matrix,tol=1.0e-10):
    mat=np.array(matrix)
    if len(mat)==0:
        return mat #empty matrix
    rank=np.linalg.matrix_rank(mat)
    
    row=len(mat)
    col=len(mat[0])
    
    #i:row j:col
    i,j,k=0,0,0

    # error retry
    retry = 0
    
    #kはreduceが終わった行の数。k=rankになったら終了
    while True:
        
        #k-1行まで掃き出しが終わっている
        #k行以下に注目する
        i,j=k,k
        #k行以下で、非零成分が最初に現れる列をj
        while all([abs(elem)<tol for elem in mat[k:,j]]):
            j+=1
            
        #k行以下で、非ゼロ成分を最も左に持つ行の一つをi行とする
        while abs(mat[i,j])<tol:
                i=i+1

        #ピボットを1に
        mat[i]=mat[i]/mat[i,j]



        #他の列のj列を0にする
        for r in range(row):
            if r==i:
                continue
            mat[r]=mat[r]-mat[i]*mat[r,j]

        #i列とk列を入れ替える
        mat[[k,i]]=mat[[i,k]]
        
        k+=1
        if k==rank:
            return mat  
        
        # check numerical error
        check_error = mat.flatten()
        check_error = check_error[np.where(np.abs(check_error)>tol/10)]
        check_error = check_error[np.where(np.abs(check_error)<tol*10)]
        if len(check_error)>0:
            raise ValueError('numerical error')
        
        # zero overwrite
        mat = np.where(np.abs(np.abs(mat))<tol,0,mat)
