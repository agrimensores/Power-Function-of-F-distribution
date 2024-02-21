import math

def powernp(alfa,power,h,f):

##
#powernp :  Gives the non-centrality parameter for a power function
#           with (alfa, power, h and f)
#--------------------------------------------------------------------
#OUTPUT  :
#           np      = Lower Bound of the Non-Centrality Parameter
#INPUT   :
#           alfa    = Significance level
#           power   = Power of the test
#           h       = First degrees of freedom (numerator)
#           f       = Second degrees of freedom (denominator)
#
#NOTE      : Take f=1e8 for computing the noncentrality parameter of
#            chi-square distribution
#--------------------------------------------------------------------
#Written by Cuneyt Aydin (2021-2023)
#           Yildiz Technical University-Civil Engineering Faculty
#           Geodesy Division, Istanbul/Turkey
#           Contact: caydin78@gmail.com;caydin@yildiz.edu.tr
#           https://avesis.yildiz.edu.tr/caydin
#--------------------------------------------------------------------
#References:
#[1] Aydin C and Gunes O (2023), Power function of F-distribution:
#Revisiting its computation and solution for geodetic studies
#--------------------------------------------------------------------

    epsL=1e-20 #convergence error

    if f>1e12: #this is for saving time. To see whether it works even for f>10^12, remove this if-end
        f=1e12

    if power<alfa:
        np='power should be bigger than or equal to alfa'
        return np

    if alfa<1e-5 or power>0.9999:
        np='alfa should be bigger than 10^-5, power should be smaller than 0.9999'
        return np

    if h<=0 or f<=0:
        np='h and f should be bigger than zero'
        return np

    if (math.floor(h) - h < 0) or (math.floor(f) - f < 0):
        np='h and f should be integers'
        return np


    F=f_inv(1-alfa,h,f)
    a=h/2
    b=f/2
    x=h*F/(f+h*F)
    
    done1=0
    say=0

    delta=LNPapp_c(F,alfa,power,h,f)/2
    if (delta*2)>1e6:
        np='Solution is out of capacity'
        return np

    if delta<1:
        delta=math.log(1-alfa)-math.log(1-power)
        if delta==0:
            delta=1e-24
            

    while done1==0:
        done=0
        i=0
        logT=math.log(1-alfa)
        m=1-alfa
        sumP=math.exp(logT-delta)

        while done==0:
            i=i+1
            c=betaa(x,a+i-1,b)
            m=m-c
            if m>0:
                logT=logT+math.log(delta)+math.log(m)-math.log(m+c)-math.log(i) #log of each term T(i+1)
                P=math.exp(logT-delta); #each term P(i+1);
                sumP=sumP+P

                if i>delta:
                    if P<epsL:
                        done=1
            if m<0:
                done=1
            
        say=say+1

        sumP=math.exp(math.log(sumP))

        k_k=math.log(sumP/(1-power))
        D_D=delta

        if say==1:
            k_1=k_k
            D_1=D_D
    
        dk=k_1-k_k
        k_1=k_k
        dD=D_D-D_1
        D_1=D_D

        if say<2:
            aa=1
        else:
            aa=dD/dk
    
        delta=delta+aa*k_k
    
        if (delta*2)>1e6:
            np='Solution is out of capacity'
            return np
    
        #for convergence problems redefine initial value with a better algorithm-------------------
        if delta<0:
            delta=LNPapp(alfa,power,h,f,1e-5)/2
            say=0
    
        if say==15:
            delta=LNPapp(alfa,power,h,f,1e-5)/2
            say=0
    
        #------------------------------------------------------------------------------------------

        if abs(sumP+power-1)<1e-5:
            done1=1    

    np=delta*2
    return np
    

#------------------------------------------------------------------------------
def betaa(x, a, b):
    d = (a) * math.log(x) + b * math.log(1 - x) + gln(a + b) - gln(a + 1) - gln(b)
    c = math.exp(d)

    return c
#------------------------------------------------------------------------------

def gln(s):
    if s > 100:
        g = (s - 0.5) * math.log(s - 1) - (s - 1) + 0.5 * math.log(2 * math.pi) + 1 / (12 * (s - 1)) - 1 / (
                360 * (s - 1) ** 3) + 1 / (1260 * (s - 1) ** 5)
    else:
        h = float(s) - int(s)
        if h == 0.5:
            s = s - 0.5
            g = math.log(math.sqrt(math.pi)) - s * math.log(2)
            for i in range(1, int(s)+1):
                g = g + math.log(2 * s - 2 * i + 1)
        if h == 0:
            g = 0
            for i in range(1, int(s)):
                g = g + math.log(i)

    return g
#------------------------------------------------------------------------------

def f_inv(cl, h, f):
#Upper alpha percentage value of F-distribution
#CL=1-alpha
#h=first degrees of freedom
#f=second degrees of freedom
#Note=The program is stable for CL<0.99 and f>1. For CL getting close to 1
#the percentage value go to infitinity, and so big errors may occur. However,
#these errors do not affect on the power of the test which is the main aim 
#of the program "powerf"
    fa = 0
    fb = 1e10
    eps2 = 1e-12
    if cl > 0.99:
        eps2 = 1e-5
    fo = (fa + fb) / 2
    done = 0
    while done == 0:
        c = cdf0(fo, h, f)
        if c < cl:
            fa = fo
        else:
            fb = fo
        fo = (fa + fb) / 2
        if abs(fa - fb) < eps2:
            done = 1
            ff = fo

    return ff
#------------------------------------------------------------------------------
def cdf0(fo, h, f):
# cdf0: CDF of central F-distribution
#      Fo=value
#      h=first degrees of freedom
#      f=second degrees of freedom
    a = h / 2
    b = f / 2
    x = h * fo / (f + h * fo)
    if x * (a + b + 2) < (a + 1):
        p = icbeta(x, a, b)
    else:
        p = 1 - icbeta(1 - x, b, a)

    return p

#------------------------------------------------------------------------------

def icbeta(x, p, q):
#icbeta:  Incomplete beta function ratio while x*(p+q-2)<(p-1).
#         Otherwise, pp=1-icbeta(1-x,q,p). 
#
#         This function is created from Phien (1990) for computing the upper 
#         percentage value of F-distribution. It gives the most efficient 
#         solution among many algorithms for incomplete beta function.
#
#Created: by C. Aydin (on April-2021)
#         Yildiz Technical University, Geodesy Division
#         Istanbul, Turkey
#         caydin78@gmail.com;caydin@yildiz.edu.tr
#         https://avesis.yildiz.edu.tr/caydin
#
#References:
#     [1] Phien, HN, (1990). A note on the computation of the Incomplete
#         beta function, Adv. Eng. Software, 12/1, pp.39-44.
#     [2] Abramowitz, M and Stegun, IA, (1972). Handbook of mathematical
#         functions, Dover Publications, New York.

    eps1 = 1e-12
    an = 1
    bn = 1
    az = 1
    qab = p + q
    qap = p + 1
    qan = p - 1
    bz = 1 - qab * x / qap
    done = 0
    n = 0
    bab = gln(p) + gln(q) - gln(p + q)
    c = math.exp(p * math.log(x) + q * math.log(1 - x) - bab) / p  # constant
    while done == 0:
        n = n + 1
        d = n * (q - n) * x / ((qan + 2 * n) * (p + 2 * n))
        ap = az + d * an
        bp = bz + d * bn
        d = -(p + n) * (qab + n) * x
        d = d / ((p + 2 * n) * (qap + 2 * n))
        app = ap + d * az
        bpp = bp + d * bz
        aold = az
        an = ap / bpp
        bn = bp / bpp
        az = app / bpp
        bz = 1
        if abs(az - aold) < (eps1 * abs(az)):
            done = 1
            pp = az
    pp = pp * c

    return pp

#------------------------------------------------------------------------------

def LNPapp_c(c_value,alfa,power,h,f):

# Approximation for non-centrality parameter computation using
# normal-distribution approximation

    eps2=1e-5
    la=0
    lb=1

    while lb-la>eps2:

        lo=(la+lb)/2
        powerd=ffcdf(c_value,h,f,lo)-(1-power)
        if abs(powerd)<eps2:
            break
        if powerd<0:
            lb=lo
        else:
            la=lo
            if lo>lb*0.9:
                lb=lb*2
        
    delta=(la+lb)/2
    return(delta)

#------------------------------------------------------------------------------

def ffcdf(x,h,f,L):
#Approximate solution for CDF of non-central F-distribution 
#--------------------------------------------------------------------------
#References:
#[1] Severo NC and Zelen M (1960), Normal approximation to the chi-square
#and non-central F probability functions, Biometrika, 47(3-4),411-416.
#[2] Tiku ML, (1965),Laguerre series forms of non-central chi-square and 
#F distributions, Biometrika,52 (3-4), 415-427.
#[3] Vazquez-Leal H, Castaneda-Sheissa R, Filobello-Nino U, 
#Sarmiento-Reyes A, Sanchez Orea H, (2012), High accurate simple 
#approximation of normal distribution integral, Mathematical Problems in 
#Engineering, vol. 2012, Article ID 124029, 22 pages,
#https://doi.org/10.1155/2012/124029.
#--------------------------------------------------------------------------

#approximation of Severo and Zelen (1960), Biometrica
    x1=(h*x/(h+L))**(1/3)
    x2=1-2/(9*f)
    x4=2*(h+2*L)/(9*(h+L)**2)
    x3=1-x4
    x5=(1-x2)*((h*x/(h+L))**(2/3))
    xx=(x1*x2-x3)/math.sqrt(x4+x5)

#approximation of Tiku 1965, Biometrika
#     x1=math.sqrt(2*f-1);
#     x2=h*x/f;
#     x3=2*(h+L);
#     x4=(h+2*L)/(h+L);
#     xx=(x1*math.sqrt(x2)-math.sqrt(x3-x4))/math.sqrt(x2+x4);

#approximation of cdf of ....mathematical problems in engineering (high
#accurate simple approximation of normal distribution
    b1=math.exp(-358*xx/23+111*math.atan(37*xx/294))
    bbo=1/(b1+1)
    c=bbo
    return(c)

#------------------------------------------------------------------------------

def LNPapp(alfa,power,h,f,eps2):
    
# Second (but better!) approximation for non-centrality parameter using
# powerf original function
# This is a bisection algorithm; however, the upper_limit of the 
# searching interval is enlarged during iterative procedure. This is
# therefore a rough estimation of the non-centrality parameter. For an
# exact solution, the upper_limit should be decided before iteration.

    la=0
    lb=1
    while lb-la>eps2:
        lo=(la+lb)/2
        powerd=(1-powerf(alfa,h,f,lo))-(1-power)
        if abs(powerd)<eps2:
            break
    
        if powerd<0:
            lb=lo
        else:
            la=lo
            if lo>lb*0.9:
                lb=lb*2
    delta=(la+lb)/2
    return(delta)
#------------------------------------------------------------------------------

def powerf(alfa, h, f, np):
#
#powerf :yields the power of the test for statistics following
#        F-distribution (complementary of CDF of non-central F-distribution)
#------------------------------------------------------------------------------
#OUTPUT :
#       p    = Power of the test
#INPUT  :
#       alfa = Significance level
#       h    = First degrees of freedom
#       f    = Second degrees of freedom
#       NP   = Non-centrality parameter
#
#------------------------------------------------------------------------------
#Written by Cuneyt Aydin (2020-2021, 2023)
#       Yildiz Technical University-Civil Engineering Faculty
#       Geodesy Division, Istanbul/Turkey
#       Contact: caydin78@gmail.com;caydin@yildiz.edu.tr
#       https://avesis.yildiz.edu.tr/caydin
#------------------------------------------------------------------------------
#References:
#[1] Aydin C and Gunes O (2023), Power function of F-distribution:
#Revisiting its computation and solution for geodetic studies.
#------------------------------------------------------------------------------
    if f>1e12: #this is for saving time. To see whether it works even for f>10^10, remove this line
        f=1e12

    done = 0
    i = 0
    epsL = 1e-20
    ff=f_inv(1-alfa,h,f) #upper alfa percentage value of F-distribution
    if np<1e-20:
        np=1e-24
        
    delta = np / 2
    a = h / 2
    b = f / 2
    x = h * ff / (f + h * ff)
    t = 1 - alfa
    logt = math.log(t)
    m = 1 - alfa  # the first term
    sump = math.exp(logt-delta)
    while done == 0:
        i = i + 1
        c = betaa(x, a + i - 1, b)
        m = m - c

        if m > 0:
            logt = logt + math.log(delta) + math.log(m) - math.log(m + c) - math.log(i)
            pp = math.exp(logt - delta)
            sump = sump + pp
            if i > delta:
                if pp < epsL:
                    done = 1
        if m <= 0:
            done = 1
        if sump < 1e-20:
            sump=1e-20
    p = math.exp(math.log(sump))  # CDF
    p = 1 - p  # power

    return p
#------------------------------------------------------------------------------



