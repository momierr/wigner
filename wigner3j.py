import numpy as np 
import math 
import warnings
from scipy.special import gammaln


def delta(j123):
    j1 = j123[0]
    j2 = j123[1]
    j3 = j123[2]
    return np.sqrt(math.factorial(j1 + j2 - j3) * math.factorial(j1 - j2 + j3) * math.factorial(-j1 + j2 + j3) / math.factorial(j1 + j2 + j3 + 1))
    
def wigner3j(j123, m123):
    
    j1 = j123[0]
    j2 = j123[1]
    j3 = j123[2]
    
    m1 = m123[0]
    m2 = m123[1]
    m3 = m123[2]
    
    
    j123 = np.array(j123)
    m123 = np.array(m123)
    
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
    w3j *= delta(j123)
    
    print(j1-m1)
    print(j1+m1)
    print(j2-m2)
    print(j2+m2)
    print(j3-m3)
    print(j3+m3)
    
    w3j *= np.sqrt(np.exp(gammaln(j1-m1 +1)) * np.exp(gammaln(j1+m1+1)) * np.exp(gammaln(j2-m2+1)) * np.exp(gammaln(j2+m2+1)) *np.exp(gammaln(j3-m3+1)) *np.exp(gammaln(j3+m3+1)))
        
    sum = 0
    
    
    for k in range(tmin, tmax + 1):
        X = np.exp(gammaln(k+1)) * np.exp(gammaln(j3 - j2 + k + m1+1)) * np.exp(gammaln(j3 - j1 + k - m2+1)) * np.exp(gammaln(j1 + j2 - j3 - k +1)) * np.exp(gammaln(j1 - k - m1 +1)) * np.exp(gammaln(j2 - k + m2 +1))
        sum += ((-1)**k)/X
    
    w3j *= sum

    
    return w3j






    
