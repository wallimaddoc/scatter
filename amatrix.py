import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm

from scipy.special import jv , yv

import scipy as sp
from scipy import linalg as sla;


import numpy as np
from numpy import pi
from numpy import linalg as la

from math import cos, sin

## perhaps
from numpy.linalg import norm
from numpy.ma import conjugate
from ndg.httpsclient.utils import fetch_from_url
conj = np.frompyfunc(conjugate,1,1)

def jh (n,x):
    return jv(n,x)+1j*yv(n,x)

def diff(fun,x,h):
    return (fun(x)-fun(x-h))/h;
def besseldiff(bessel,n,x):
    return (bessel(n-1,x) - bessel(n+1,x))/2.0
    
def aentryh (D,C,x,m,n,h):
    result = n(x)*besseldiff(D, m, x) * C(m,x-h)-n(x-h)* D(m,x) * besseldiff(C, m, x-h)
    result = result * pi /2.0/1j*x;
    return result;   
def aentry (D,C,x,m,n,ndiff):
    result = besseldiff(D, m, x)*besseldiff(C, m, x)+(1-m*m/x/x/n(x)/n(x))*C(m,x)*D(m,x)
    result = result * pi /2.0/1j*x*x*ndiff(x)*n(x);
    return result; 

bj = np.frompyfunc(jv,2,1)
by = np.frompyfunc(yv,2,1)
bh = np.frompyfunc(jh,2,1)

na = np.frompyfunc(cos,1,1)
nadiff = np.frompyfunc(sin,1,1)
m=0.0;
R=1.0
N=10000.0
h = R/N
r = np.append(np.arange(0, R-h/2, R/(N-1)),R)
pos = 50;
x = r[pos];
x=1;

v = aentryh(bj, bj, x, m, na,h)
w = aentry(bj, bj, x, m, na,nadiff)
print v/h
print w;
print (v/h/w)






