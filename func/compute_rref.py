import numpy as np

def compute_rref(matrix):
    import numpy as np
    mat=np.array(matrix)
    if len(mat)==0:
        return mat #empty matrix
    rank=np.linalg.matrix_rank(mat)
    
    row=len(mat)
    col=len(mat[0])
    
    #i:row j:col
    i,j,k=0,0,0
    
    #kはreduceが終わった行の数。k=rankになったら終了
    while True:
        
        #k-1行まで掃き出しが終わっている
        #k行以下に注目する
        i,j=k,k
        #k行以下で、非零成分が最初に現れる列をj
        while all([abs(elem)<1.0e-10 for elem in mat[k:,j]]):
            j+=1
            
        #k行以下で、非ゼロ成分を最も左に持つ行の一つをi行とする
        while abs(mat[i,j])<1.0e-10:
                i=i+1

        #ピボットを1に
        mat[i]=mat[i]/mat[i,j]



        #他の列のj列を0にする
        for r in range(row):
            if r==i:
                continue
            mat[r]=mat[r]-mat[i]*mat[r,j]

        #i列とk列を入れ替える
        exc1=mat[i].copy()
        exc2=mat[k].copy()
        mat[k]=exc1
        mat[i]=exc2
        
        k+=1
        if k==rank:
            return mat  