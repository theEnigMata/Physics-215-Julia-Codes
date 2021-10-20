"""@author: User"""

from numpy import arange, linspace, sqrt, sum, pi, sin
from scipy.special import ellipk
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt

plt.close('all')

def double_factorial(n):
    if n < 2:
        return 1
    else:
        return n*double_factorial(n-2)

def K(k):
    result = (pi/2) * sum([(double_factorial(2*n-1)/double_factorial(2*n))**2 * k**n for n in range(150)])
    
    return result

def main():
    k = arange(0, 1, 0.005)
    
    ellip = []
    
    for kdx in k:
        ellip.append(K(kdx))
        
    ellip1 = ellipk(k)
    
    fig, ax = plt.subplots()
    ax.plot(k, ellip)
    # ax.plot(k, ellip1)
    plt.show()
    
if __name__ == '__main__':
    main()