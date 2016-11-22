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
conj = np.frompyfunc(conjugate,1,1)

x=2.0-0.3j;y=1.0+4j;
testMatrix = np.array([[0,x],[y,0.0]]);
tm = np.array([[0,-x],[-y,0.0]]);
tm =  sla.expm(tm)
v = expMatrix = sla.expm(testMatrix)
t11 = np.cos(1j*np.sqrt(x*y));
t12 = -1j*np.sqrt(x/y)*np.sin(1j*np.sqrt(x*y));
t21 = -1j*np.sqrt(y/x)*np.sin(1j*np.sqrt(x*y));

print v;
t2 = np.array([[t11,t12],[t21,t11]]);
print t2
print np.dot(v,tm);
w =1;
