import numpy as np
from scipy import integrate
from scipy.special import jv , yv
import matplotlib.pyplot as plot

def jh (n,x):
    return jv(n,x)+1j*yv(n,x)
bj = np.frompyfunc(jv,2,1)
by = np.frompyfunc(yv,2,1)
bh = np.frompyfunc(jh,2,1)


def solvr(Y, t):
#    return [Y[1], - c*c*np.power(t,2*n)* (1.0 - m*m/c/c*(n+1.0)*(n+1.0)/np.power(t,2.0*(n+1.0)))* Y[0]- 1.0/t* Y[1]];

    return [Y[1],  (m*m*(n+1.0)*(n+1.0)/np.power(t,2.0)-c*c*np.power(t,2*n))* Y[0]- 1.0/t* Y[1]];
def ex (t):
    return jv(m,c/(n+1.0)*np.power(t,(n+1)));
def ex_diff (t):
    return (jv(m-1,c/(n+1.0)*np.power(t,(n+1)))-jv(m+1,c/(n+1.0)*np.power(t,(n+1.0))))*c*np.power(t,n)/2.0;


ex_sol = np.frompyfunc(ex,1,1)
ex_sol_diff = np.frompyfunc(ex_diff,1,1)

n=1.5;
m=2;
c=2.13;
a=0.1;
a_t = np.arange(a, 10.0, 0.01)
asol = integrate.odeint(solvr, [ex_sol(a), ex_sol_diff(a)], a_t);
sol = [0]*len(asol);
for y in range(len(asol)):
    sol[y]=1*asol[y][0];
w = ex_sol(a_t);
plot.plot(a_t,sol,'r+',a_t,w,'g--');
plot.show()
print(sol)
print w;
