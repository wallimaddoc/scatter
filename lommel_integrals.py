import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm


from scipy.special import jv , yv
import scipy.integrate as integrate

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
from gtk.keysyms import integral
conj = np.frompyfunc(conjugate,1,1)

def jh (n,x):
    return jv(n,x)+1j*yv(n,x)

def diff(fun,x,h):
    return (fun(x)-fun(x-h))/h;


def besseldiff(bessel,n,x,ord=1):
    if (ord==1):
        return (bessel(n-1,x) - bessel(n+1,x))/2.0
    else: 
        return besseldiff(bessel,n-1,x,ord-1)/2.0 - besseldiff(bessel,n+1,x,ord-1)/2.0


bj = np.frompyfunc(jv,2,1)
by = np.frompyfunc(yv,2,1)
bh = np.frompyfunc(jh,2,1)

m=8.43;
R=1.0
a=0.0000001;
b=2;
t=34.3;
k=t;
mu=0;



def integrand_plus_real (x,C,D,m,t,k,mu):
    return np.real(np.power(x,mu+1)*C(m+mu,t*x)*D(m,k*x));
def integrand_plus_imag (x,C,D,m,t,k,mu):
    return np.imag(np.power(x,mu+1)*C(m+mu,t*x)*D(m,k*x));
def integrand_minus_real (x,C,D,m,t,k,mu):
    return np.real(np.power(x,mu+1)*C(m-mu,t*x)*D(m,k*x));
def integrand_minus_imag (x,C,D,m,t,k,mu):
    return np.imag(np.power(x,mu+1)*C(m-mu,t*x)*D(m,k*x));

def primitive_plus (x,C,D,m,t,k,mu):
    return np.power(x,mu+2)/2.0/(mu+1)*(C(m+mu,t*x)*D(m,k*x)-C(m+mu+1,t*x)*D(m-1,k*x));
def primitive_minus (x,C,D,m,t,k,mu):
    return np.power(x,mu+2)/2.0/(mu+1)*(C(m-mu,t*x)*D(m,k*x)-C(m-mu-1,t*x)*D(m+1,k*x));

def primitive_strich (x,C,D,m,t,k,mu):
    return x*x/2.0*C(m,t*x)*D(m,k*x)-(m*m/t/t/2.0)*C(m,t*x)*D(m,k*x)+x*x/2.0*besseldiff(C,m,t*x)*besseldiff(D,m,t*x);

int1 = bj;
int2 = bh;

def check_plus (): 
    real_plus = integrate.quad(integrand_plus_real,a,b,args=(int1,int2,m,t,k,mu));
    imag_plus = integrate.quad(integrand_plus_imag,a,b,args=(int1,int2,m,t,k,mu));
    stammfkt_num_plus = real_plus[0]+1j*imag_plus[0];
    stammfkt_plus = primitive_plus(b,int1,int2,m,t,k,mu)-primitive_plus(a,int1,int2,m,t,k,mu);
    stammfkt_2 = primitive_strich(b,int1,int2,m,t,k,mu)-primitive_strich(a,int1,int2,m,t,k,mu);

    print stammfkt_2;
    print stammfkt_num_plus;
    print stammfkt_plus;

def check_minus ():
    real_minus = integrate.quad(integrand_minus_real,a,b,args=(int1,int2,m,t,k,mu));
    imag_minus = integrate.quad(integrand_minus_imag,a,b,args=(int1,int2,m,t,k,mu));
    stammfkt_num_minus = real_minus[0]+1j*imag_minus[0];
    stammfkt_minus = primitive_minus(b,int1,int2,m,t,k,mu)-primitive_minus(a,int1,int2,m,t,k,mu);
    
    print stammfkt_num_minus;
    print stammfkt_minus;
    
def check_limis ():
    int1 = bj;
    int2 = bj;

    r= 55100.0;
    #stammfkt_minus = np.linalg.norm(primitive_minus(r,int1,int2,m,t,k,mu));
    #stammfkt_plus = np.linalg.norm(primitive_plus(r,int1,int2,m,t,k,mu));
    c= -1/pi/t/t*np.exp(2j*t*r)*np.exp(-1j*pi*m);
    #c= 1j*1/pi*m/t/t;

    print primitive_strich(r,int1,int2,m,t,k,mu)/c;

    #print np.linalg.norm(primitive_strich(r,int1,int2,m,t,k,mu));
    #r= 12899393;
    #c= -1/pi/t/t*np.exp(2j*t*r)*np.exp(-1j*pi*m);
    #c= 1j*1/pi*m/t/t;


    #print primitive_strich(r,int1,int2,m,t,k,mu)/c;

    #print np.linalg.norm(primitive_strich(r,int1,int2,m,t,k,mu));

    #print stammfkt_num;


    #print stammfkt_minus;
    #print stammfkt_plus;

def check_CStrich_D ():
    
    # Investigate
    # integral ( s^2 C'(m,tx) D(m,tx)
    
    def integrand_CStrich_real (x,C,D,m,t,k):
        return np.real(np.power(x,2)*besseldiff(C,m,t*x)*D(m,k*x));
    def integrand_CStrich_imag (x,C,D,m,t,k):
        return np.imag(np.power(x,2)*besseldiff(C,m,t*x)*D(m,k*x));

    CStrich_real = integrate.quad(integrand_CStrich_real,a,b,args=(int1,int2,m,t,k));
    CStrich_imag = integrate.quad(integrand_CStrich_imag,a,b,args=(int1,int2,m,t,k));
    
    def primitive_mu_1_0_minnus_mu_1_1 (x,C,D,m,t,k):
        mu=1;
        result1=0;
        #result1 = C(m-mu,t*x)*D(m,k*x);
        #result1-= C(m+mu,t*x)*D(m,k*x);
        #result1 = 2 * besseldiff(C,m,t*x)*D(m,k*x);
        
        #result2 = -C(m-mu-1,t*x)*D(m+1,k*x);
        #result2 = -(-C(m,t*x) + 2 * (m-1)/t/x * C(m-1,t*x) ) * D(m+1,k*x);
        #result2 = -(-C(m,t*x) + 2 * (m-1)/t/x * C(m-1,t*x) ) * (m/t/x*D(m,t*x)- besseldiff(D,m,t*x));
        #result2 = -(-C(m,t*x) + 2 * (m-1)/t/x * (m/t/x*C(m,t*x) + besseldiff(C,m,t*x)) ) * (m/t/x*D(m,t*x)- besseldiff(D,m,t*x));
        #result2 = (C(m,t*x) - 2 * (m-1)/t/x * (m/t/x*C(m,t*x) + besseldiff(C,m,t*x)) ) * (m/t/x*D(m,t*x)- besseldiff(D,m,t*x));
        #result2 = ((1-2 * (m-1)*m/t/t/x/x)*C(m,t*x) - 2 * (m-1)/t/x *besseldiff(C,m,t*x) ) * (m/t/x*D(m,t*x)- besseldiff(D,m,t*x));
       
        result2 =0;
        #  result2 = (1-2 * (m-1)*m/t/t/x/x)*C(m,t*x) * m/t/x*D(m,t*x);
        #  result2-= (1-2 * (m-1)*m/t/t/x/x)*C(m,t*x) * besseldiff(D,m,t*x);
        #  result2-= 2 * (m-1)/t/x *besseldiff(C,m,t*x) * m/t/x*D(m,t*x);
        #  result2+= 2 * (m-1)/t/x *besseldiff(C,m,t*x) * besseldiff(D,m,t*x);

        #result2+= C(m+mu+1,t*x)*D(m-1,k*x);
        #result2+= (-C(m,t*x) + 2 * (m+1)/t/x * C(m+1,t*x) ) * D(m-1,k*x);
        #result2+= (-C(m,t*x) + 2 * (m+1)/t/x * C(m+1,t*x) ) * (m/t/x*D(m,t*x) + besseldiff(D,m,t*x));
        #result3 = (-C(m,t*x) + 2 * (m+1)/t/x * (m/t/x*C(m,t*x)- besseldiff(C,m,t*x))) * (m/t/x*D(m,t*x) + besseldiff(D,m,t*x));

        result3 =0;
        
        # result3 = (-1 + 2 * (m+1)*m/t/t/x/x ) * C(m,t*x)  * m/t/x*D(m,t*x);
        # result3+= (-1 + 2 * (m+1)*m/t/t/x/x ) * C(m,t*x)  * besseldiff(D,m,t*x);
        # result3-= 2 * (m+1)/t/x * besseldiff(C,m,t*x) * m/t/x*D(m,t*x);
        # result3-= 2 * (m+1)/t/x * besseldiff(C,m,t*x) * besseldiff(D,m,t*x);
        
        #result4 = (2 * (m-1)/t/x - 2 * (m+1)/t/x ) *besseldiff(C,m,t*x) * besseldiff(D,m,t*x);
        result4 = -2/t/x  *besseldiff(C,m,t*x) * besseldiff(D,m,t*x);
        result4+= ( 1 - 2 * m*m/t/t/x/x ) * ( besseldiff(C,m,t*x) * D(m,t*x) - C(m,t*x) * besseldiff(D,m,t*x));
        result4+= 2 * m*m/t/t/t/x/x/x   * C(m,t*x) * D(m,t*x);
          
        result= np.power(x,mu+2)/2.0/(mu+1)*(result1+result2+result3+result4);
        return result;
        #return np.power(x,mu+2)/2.0/(mu+1)*(C(m-mu,t*x)*D(m,k*x)-C(m-mu-1,t*x)*D(m+1,k*x))-np.power(x,mu+2)/2.0/(mu+1)*(C(m+mu,t*x)*D(m,k*x)-C(m+mu+1,t*x)*D(m-1,k*x));

    
    
    print (CStrich_real[0]+1j*CStrich_imag[0]);
    print (primitive_minus(b,int1,int2,m,t,k,1)-primitive_minus(a,int1,int2,m,t,k,1)-primitive_plus(b,int1,int2,m,t,k,1)+primitive_plus(a,int1,int2,m,t,k,1))/2.0;
    print primitive_mu_1_0_minnus_mu_1_1(b,int1,int2,m,t,k)-primitive_mu_1_0_minnus_mu_1_1(a,int1,int2,m,t,k)


def check_CStrichStrich_D ():
    
    # Investigate
    # integral ( s^2 C'(m,tx) D(m,tx)
    
    def integrand_CStrichStrich_real (x,C,D,m,t,k):
        return np.real(np.power(x,3)*besseldiff(C,m,t*x,2)*D(m,k*x));
    def integrand_CStrichStrich_imag (x,C,D,m,t,k):
        return np.imag(np.power(x,3)*besseldiff(C,m,t*x,2)*D(m,k*x));

    CStrich_real = integrate.quad(integrand_CStrichStrich_real,a,b,args=(int1,int2,m,t,k));
    CStrich_imag = integrate.quad(integrand_CStrichStrich_real,a,b,args=(int1,int2,m,t,k));
    
    def primitive_mu_2_0 (x,C,D,m,t,k):
        mu=1;
        result1=0;
        #result1 = C(m-mu,t*x)*D(m,k*x);
        #result1-= C(m+mu,t*x)*D(m,k*x);
        #result1 = 2 * besseldiff(C,m,t*x)*D(m,k*x);
        
        #result2 = -C(m-mu-1,t*x)*D(m+1,k*x);
        #result2 = -(-C(m,t*x) + 2 * (m-1)/t/x * C(m-1,t*x) ) * D(m+1,k*x);
        #result2 = -(-C(m,t*x) + 2 * (m-1)/t/x * C(m-1,t*x) ) * (m/t/x*D(m,t*x)- besseldiff(D,m,t*x));
        #result2 = -(-C(m,t*x) + 2 * (m-1)/t/x * (m/t/x*C(m,t*x) + besseldiff(C,m,t*x)) ) * (m/t/x*D(m,t*x)- besseldiff(D,m,t*x));
        #result2 = (C(m,t*x) - 2 * (m-1)/t/x * (m/t/x*C(m,t*x) + besseldiff(C,m,t*x)) ) * (m/t/x*D(m,t*x)- besseldiff(D,m,t*x));
        #result2 = ((1-2 * (m-1)*m/t/t/x/x)*C(m,t*x) - 2 * (m-1)/t/x *besseldiff(C,m,t*x) ) * (m/t/x*D(m,t*x)- besseldiff(D,m,t*x));
       
        result2 =0;
        #  result2 = (1-2 * (m-1)*m/t/t/x/x)*C(m,t*x) * m/t/x*D(m,t*x);
        #  result2-= (1-2 * (m-1)*m/t/t/x/x)*C(m,t*x) * besseldiff(D,m,t*x);
        #  result2-= 2 * (m-1)/t/x *besseldiff(C,m,t*x) * m/t/x*D(m,t*x);
        #  result2+= 2 * (m-1)/t/x *besseldiff(C,m,t*x) * besseldiff(D,m,t*x);

        #result2+= C(m+mu+1,t*x)*D(m-1,k*x);
        #result2+= (-C(m,t*x) + 2 * (m+1)/t/x * C(m+1,t*x) ) * D(m-1,k*x);
        #result2+= (-C(m,t*x) + 2 * (m+1)/t/x * C(m+1,t*x) ) * (m/t/x*D(m,t*x) + besseldiff(D,m,t*x));
        #result3 = (-C(m,t*x) + 2 * (m+1)/t/x * (m/t/x*C(m,t*x)- besseldiff(C,m,t*x))) * (m/t/x*D(m,t*x) + besseldiff(D,m,t*x));

        result3 =0;
        
        # result3 = (-1 + 2 * (m+1)*m/t/t/x/x ) * C(m,t*x)  * m/t/x*D(m,t*x);
        # result3+= (-1 + 2 * (m+1)*m/t/t/x/x ) * C(m,t*x)  * besseldiff(D,m,t*x);
        # result3-= 2 * (m+1)/t/x * besseldiff(C,m,t*x) * m/t/x*D(m,t*x);
        # result3-= 2 * (m+1)/t/x * besseldiff(C,m,t*x) * besseldiff(D,m,t*x);
        
        #result4 = (2 * (m-1)/t/x - 2 * (m+1)/t/x ) *besseldiff(C,m,t*x) * besseldiff(D,m,t*x);
        result4 = -2/t/x  *besseldiff(C,m,t*x) * besseldiff(D,m,t*x);
        result4+= ( 1 - 2 * m*m/t/t/x/x ) * ( besseldiff(C,m,t*x) * D(m,t*x) - C(m,t*x) * besseldiff(D,m,t*x));
        result4+= 2 * m*m/t/t/t/x/x/x   * C(m,t*x) * D(m,t*x);
          
        result= np.power(x,mu+2)/2.0/(mu+1)*(result1+result2+result3+result4);
        return result;
        #return np.power(x,mu+2)/2.0/(mu+1)*(C(m-mu,t*x)*D(m,k*x)-C(m-mu-1,t*x)*D(m+1,k*x))-np.power(x,mu+2)/2.0/(mu+1)*(C(m+mu,t*x)*D(m,k*x)-C(m+mu+1,t*x)*D(m-1,k*x));

    
    
    print (CStrich_real[0]+1j*CStrich_imag[0]);
    print primitive_minus(b,int1,int2,m,t,k,2)-primitive_minus(a,int1,int2,m,t,k,2);
    print primitive_mu_2_0(b,int1,int2,m,t,k)-primitive_mu_2_0(a,int1,int2,m,t,k)

#check_limis();

#check_CStrich_D();
check_CStrichStrich_D();




N=10000.0
h = R/N
r = np.append(np.arange(0, R-h/2, R/(N-1)),R)







