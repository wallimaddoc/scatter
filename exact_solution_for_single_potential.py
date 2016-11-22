import numpy as np
from scipy.special import gamma

from scipy import integrate
from scipy.special import jv , yv
import matplotlib.pyplot as plot

def jh (n,x):
    return jv(n,x)+complex(0,1)*yv(n,x)

def besseldiff(bessel,n,x,ord=1):
    if (ord==1):
        return (bessel(n-1,x) - bessel(n+1,x))/2.0
    else: 
        return besseldiff(bessel,n-1,x,ord-1)/2.0 - besseldiff(bessel,n+1,x,ord-1)/2.0

bj = np.frompyfunc(jv,2,1)
by = np.frompyfunc(yv,2,1)
bh = np.frompyfunc(jh,2,1)

bdiff = np.frompyfunc(besseldiff,4,1)

# Todo Complex faktor works not with integrate.odeint. We have to use integrate.integrate.odeint
# which has an other syntax
fak= 1j*np.pi/2.0;
fak= np.pi/2.0;

def solvr(Y, r):
#    return [Y[1], - c*c*np.power(t,2*n)* (1.0 - m*m/c/c*(n+1.0)*(n+1.0)/np.power(t,2.0*(n+1.0)))* Y[0]- 1.0/t* Y[1]];

    return [Y[1],  (m*m/np.power(r,2.0)-k*k*np.power(r,2*nu))* Y[0]- 1.0/r* Y[1]];
def ex (r):
    return fak*jv(m/(nu+1),k/(nu+1.0)*np.power(r,(nu+1)));
def ex_diff (r):
    return (jv(m/(nu+1)-1,k/(nu+1.0)*np.power(r,(nu+1)))-jv(m/(nu+1)+1,k/(nu+1.0)*np.power(r,(nu+1.0))))*k*np.power(r,nu)/2.0*fak;

def index (r):
    return k*np.power(r,nu);
def index_r (r):
    return k*np.power(r,nu)*r;
def my_real (r):
    return r.real;
def my_imag (r):
    return r.imag;


ex_sol = np.frompyfunc(ex,1,1)
ex_sol_diff = np.frompyfunc(ex_diff,1,1)
index_sol = np.frompyfunc(index,1,1)
index_sol_r = np.frompyfunc(index_r,1,1)

my_imag_function = np.frompyfunc(my_imag,1,1)
my_real_function = np.frompyfunc(my_real,1,1)


nu=2.0/2.0;
m=1.0;
k=np.sqrt(1.0+0.00*complex(0,1.0));
a=1.1;
b=a+5.52;
a_t = np.arange(a, b, 0.01)
asol = integrate.odeint(solvr, [ex_sol(a), ex_sol_diff(a)], a_t);
sol = [0]*len(asol);
sol_s = [0]*len(asol); 
one = [1]*len(asol);
minus_one = [-1]*len(asol); 
null_line =  [0]*len(asol); 
 
for y in range(len(asol)):
    sol[y]=asol[y][0];
    sol_s[y]=asol[y][1];

y = ex_sol(a_t);
ys = ex_sol_diff(a_t);
nr = index_sol(a_t);
nr_r = index_sol_r(a_t);

#plot.figure(1);
#fig1 = plot.plot(a_t,sol,'r+',a_t,y,'g--');
#plot.figure(2);
#fig2 = plot.plot(a_t,sol_s,'r+',a_t,ys,'g--');
#plot.figure(3);
#fig2 = plot.plot(a_t,nr*a_t,'r+',a_t,nr_r,'g--');

plot.figure(4);
fig2 = plot.plot(a_t,my_real_function(-a_t*ys*bj(m,nr*a_t)+nr*a_t*y*bdiff(bj,m,nr*a_t,1)),'r-+');
fig2 = plot.plot(a_t,my_imag_function(-a_t*ys*bj(m,nr*a_t)+nr*a_t*y*bdiff(bj,m,nr*a_t,1)),'y-+');

fig2 = plot.plot(a_t,my_imag_function(a_t*ys*bh(m,nr*a_t)-nr*a_t*y*bdiff(bh,m,nr*a_t,1)),'b--');
fig2 = plot.plot(a_t,my_real_function(a_t*ys*bh(m,nr*a_t)-nr*a_t*y*bdiff(bh,m,nr*a_t,1)),'g--');
fig2 = plot.plot(a_t,one);
fig2 = plot.plot(a_t,minus_one);
fig2 = plot.plot(a_t,null_line);

plot.figure(5);

# Near field # nu = 1;  k = m = 1
#a=0.1;b=a+14.52;nu=1.0;k=1.0;m=1.0;step = 0.01 # nu = 1;  k = m = 1
# Far field # nu = 1;  k = m = 1
#a=25.1;b=a+0.52;nu=1.0;k=1.0;m=1.0;step = 0.00001  # nu = 1;  k = m = 1

# Near field # nu = 2;  k = m = 1
#a=0.1;b=a+4.52;nu=2.0;k=1.0;m=1.0;step = 0.01 # nu = 2;  k = m = 1
# Far field # nu = 2;  k = m = 1
#a=25.1;b=a+0.052;nu=2.0;k=1.0;m=1.0;step = 0.00001  # nu = 2;  k = m = 1

# Near field # nu = 3;  k = m = 1
#a=0.1;b=a+4.52;nu=3.0;k=1.0;m=1.0;step = 0.01 # nu = 3;  k = m = 1
#Far field # # nu = 3;  k = m = 1
a=25.1;b=a+0.00052;nu=3.0;k=1.0;m=1.0;step = 0.00001  # nu = 3;  k = m = 1



a_t = np.arange(a, b, step)
y = ex_sol(a_t);
ys = ex_sol_diff(a_t);
one = [1]*len(a_t);
minus_one = [-1]*len(a_t); 
null_line =  [0]*len(a_t);
nr = index_sol(a_t);

# Far field # nu = 1;  k = m = 1
#far_field = -1.25*bj(m,1.0/(nu+1)*k*nr*a_t)*np.power(a_t,1.0) # nu = 1;  k = m = 1

# Far field # nu = 2;  k = m = 1
#far_field = -1.775*bj(m,1.99995/(nu+1)*k*nr*a_t)*np.power(a_t,1.5); ## nu = 2;  k = m = 1

#Far field # # nu = 3;  k = m = 1
far_field = -11.00*bj(m,2.999806/(nu+1)*k*nr*a_t)*np.power(a_t,1.5); # nu = 3;  k = m = 1



plot.figure(5);
fig2 = plot.plot(a_t,my_real_function(-a_t*ys*bj(m,nr*a_t)+nr*a_t*y*bdiff(bj,m,nr*a_t,1)),'r--');
fig2 = plot.plot(a_t,far_field,'b--');
fig2 = plot.plot(a_t,null_line);

plot.figure(6);
fig2 = plot.plot(a_t,my_real_function(-a_t*ys*bj(m,nr*a_t)+nr*a_t*y*bdiff(bj,m,nr*a_t,1)-far_field),'r--');
fig2 = plot.plot(a_t,null_line);


plot.show();

