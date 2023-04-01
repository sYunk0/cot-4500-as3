import numpy as np


def EulerMethod(f,a,b,N,w0):
    stepSize = (b-a)/N

    w = w0
    for i in range(N):
        #print(i,w)
        x = a + stepSize*i
        y = w
        w = y + stepSize*f(x,y)

    return w

def Runge_KuttaMethod(f,a,b,N,w0):
    stepSize = (b-a)/N

    w = w0
    for i in range(N):
        #print(i,w)
        x = a + stepSize*i
        x2 = a + stepSize*(i+1)
        y = w

        k1 = stepSize*f(x,y)
        k2 = stepSize*f(x + stepSize/2, y + k1/2)
        k3 = stepSize*f(x + stepSize/2, y + k2/2)
        k4 = stepSize*f(x2, y + k3)

        w = y + (1/6) * (k1+ 2*k2 + 2*k3 + k4)

    return w

def GaussianElimination(A,b):

    m,n = A.shape
    C = np.concatenate([A,np.reshape(b,(len(b),1))],axis=1)
    #print(C)

    for k in range(m):
        E_k = C[k,:]
        a_kk = C[k,k]
        for j in range(k+1,m):
            E_j = C[j,:]
            a_jk = C[j,k]

            factor = a_jk / a_kk
            C[j,:] = E_j - factor*E_k
            #print("\n\tk:{}, j:{}".format(k,j))
            #print(C)
        

    #print(C)

    A_UpperDiagonal = C[:,:n]
    b_factored = C[:,n]

    #print(A_diagonal,b_factored)

    return A_UpperDiagonal,b_factored

def backwardSubstitution(A,b):
    """
    Takes an upper triangular matrix A and a vector b and partialy solves b = Ax by making A diagonal and changing b to reflect that change
    """
    m,n = A.shape
    C = np.concatenate([A,np.reshape(b,(len(b),1))],axis=1)

    for j in range(m-1,-1,-1):
        a_jj = C[j,j]
        E_j = C[j,:]
        
        if(a_jj != 0):
            for i in range(j-1,-1,-1):
                a_ij = C[i,j]
                E_i = C[i,:]

                if(a_ij == 0):
                    factor = 0
                else:
                    factor = a_ij / a_jj

                C[i,:] = E_i - factor*E_j

    A_diagonal = C[:,:n]
    b_factored = C[:,n]
    #print(C)
    return A_diagonal,b_factored

def reduceAB(A,b):

    solutionVector = np.zeros(b.shape)
    m,n = A.shape
    for i in range(m):
        solutionVector[i] = b[i] / A[i,i]
    
    return solutionVector






def Question1():
    def f(t,y):
        return t-(y**2)

    a = 0
    b = 2
    N = 10
    
    w0 = 1

    y = EulerMethod(f,a,b,N,w0)
    print(y)

def Question2():
    def f(t,y):
        return t-(y**2)

    a = 0
    b = 2
    N = 10
    
    w0 = 1

    y = Runge_KuttaMethod(f,a,b,N,w0)
    print(y)

def Question3():
    MatrixA = np.array([[ 2,-1, 1],
                        [ 1, 3, 1],
                        [-1, 5, 4]])
    vectorB = np.array([6,0,-3])

    #print(MatrixA,vectorB)
    A,b = GaussianElimination(MatrixA,vectorB)
    A,b = backwardSubstitution(A,b)

    solutions = reduceAB(A,b)
    print(solutions)

if(__name__ == "__main__"):
    Question3()