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





if(__name__ == "__main__"):
    Question2()