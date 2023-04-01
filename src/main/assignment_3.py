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














def Question1():
    def f(t,y):
        return t-(y**2)

    a = 0
    b = 2
    N = 10
    
    w0 = 1

    y = EulerMethod(f,a,b,N,w0)
    print(y)




if(__name__ == "__main__"):
    Question1()