import numpy as np 
import math 
import warnings
from scipy.special import gammaln


def delta(a,b,c):
    return np.sqrt(math.factorial(a + b - c) * math.factorial(a - b + c) * math.factorial(-a + b + c) / math.factorial(a + b + c + 1))
    
def wigner3j(j1, j2, j3, m1, m2, m3):
    
    j123 = np.array([j1,j2,j3])
    m123 = np.array([m1,m2,m3])

    
    if (j123<0).sum() != 0:
        warnings.warn("All j values must be > 0")
        return 0
    
    if any(j123 % 0.5) or any(m123 % 0.5):
        warnings.warn("Arguments must be int or half-int")
        return 0
    
    if any((j123 - m123) % 1):
        warnings.warn("J's and m's must match")
        return 0 
    
    
    if j3 > j1 + j2 or j3 < np.abs(j1 - j2) or m1 + m2 + m3 != 0 or any(np.abs(m123) > j123):
        return 0 
        

    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2

    
    tmin = min([0, t1, t2])
    tmax = min([t3, t4, t5])
    
    
    w3j = (-1)**(j1 - j2 - m3)
    w3j *= delta(j1,j2,j3)
    
    
    w3j *= np.sqrt(np.exp(gammaln(j1-m1 +1)) * np.exp(gammaln(j1+m1+1)) * np.exp(gammaln(j2-m2+1)) * np.exp(gammaln(j2+m2+1)) *np.exp(gammaln(j3-m3+1)) *np.exp(gammaln(j3+m3+1)))
        
    sum = 0
    
    
    for k in range(tmin, tmax + 1):
        X = np.exp(gammaln(k+1)) * np.exp(gammaln(j3 - j2 + k + m1+1)) * np.exp(gammaln(j3 - j1 + k - m2+1)) * np.exp(gammaln(j1 + j2 - j3 - k +1)) * np.exp(gammaln(j1 - k - m1 +1)) * np.exp(gammaln(j2 - k + m2 +1))
        sum += ((-1)**k)/X
    
    w3j *= sum

    
    return w3j


def cg(j1,m1,j2,m2,j3,m3):
    phase = (-1)**(-j1+j2-m3)
    return phase*np.sqrt(2*j3+1)*wigner3j(j1,j2,j3,m1,m2,-m3)


def wigner6j(j1,j2,j3,j4,j5,j6):
    
    w6j = 0
    t1 = j1+j2+j3
    t2 = j1+j5+j6
    t3 = j4+j2+j6
    t4 = j4+j5+j3

    
    inf = max([t1, t2, t3, t4])

    t5 = j1 + j2 + j4 + j5
    t6 = j2 + j3 + j5 + j6
    t7 = j3 + j1 + j6 + j4
    
    sup = min([t5, t6, t7])

 
    w6j = delta(j1,j2,j3) * delta(j1,j5,j6) * delta(j4,j2,j6) * delta(j4,j5,j3)
    sum = 0
    
    for k in range(inf, sup + 1):
        X = np.exp(gammaln(k-t1+1)) * np.exp(gammaln(k-t2+1)) * np.exp(gammaln(k-t3+1)) * np.exp(gammaln(k-t4+1)) * np.exp(gammaln(t5-k+1)) * np.exp(gammaln(t6-k+1)) * np.exp(gammaln(t7-k+1))  
        sum += (-1)**k * np.exp(gammaln(k+2)) / X
        
    w6j *= sum
    
    return w6j

